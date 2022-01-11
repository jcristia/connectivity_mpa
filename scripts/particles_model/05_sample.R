
# ...

library(tidyverse)
library(sdmTMB)
library(arrow)
library(patchwork)
library(sf)
library(glue)



########################################
########################################
# Inputs

# feather files proportion label
#labels <- list('0001', '0003', '0005', '0007', '0009', '001', '002', '004', '006', '008', '01')
labels <- list('001')
# feather files directory
base_dir <- '~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all'
# coastline shapefile
land <- 'C:/Users/jcristia/Documents/GIS/MSc_Projects/MPA_connectivity/spatial/Coastline/landmask_FINAL.shp'



########################################
########################################
# Functions


### Create mesh and barrier mesh ###
bmesh <- function(data, land, res) {
  land_sf <- st_read(land, quiet=TRUE)
  mesh <- make_mesh(data, c('origin_x_km', 'origin_y_km'), cutoff = res)
  bm <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=FALSE, range_fraction=0.2)
  return(bm)
}


### Predict ###
predict_model <- function(model, newdata) {
  predictions <- predict(model, newdata)
  predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
  predictions <- st_set_crs(predictions, value = "epsg:3005")
  return(predictions)
}


### Residuals ###
resid_plots <- function(model, predictions, filename) {
  predictions$residuals <- residuals(model)
  ggplot(predictions, aes(sample=residuals)) + stat_qq_line() + stat_qq()
  ggsave(paste('qqnorm_', filename, sep=''))
  ggplot(predictions, aes(residuals)) + geom_density()
  ggsave(paste('pdf_', filename, sep=''))
  ggplot(predictions, aes(x = est, y = residuals)) + geom_point()
  ggsave(paste('resid_', filename, sep=''))
}

### R-squared ###
rsq <- function (x, y) cor(x, y) ^ 2


### Average log-likelihood ###
ll <- function(model, newdata, response_var) {

  # Get our observed and predicted data
  object <- model
  cv_data <- newdata
  response <- response_var
  predicted_obj <- predict(object, newdata = cv_data, return_tmb_obj = TRUE)
  cv_data$cv_predicted <- object$family$linkinv(predicted_obj$data$est)
  withheld_y <- predicted_obj$data[[response]]
  withheld_mu <- cv_data$cv_predicted

  # Function that establishes the gamma distribution from the modeled data and then
  # gets the probabilities of drawing the observed data from it (i.e. log likelihood)
  ll_gamma <- function(object, withheld_y, withheld_mu) {
    .shape <- exp(object$model$par[["ln_phi"]])
    stats::dgamma(x = withheld_y, shape = .shape, scale = withheld_mu / .shape, log = TRUE)
  }
  ll_binomial <- function(object, withheld_y, withheld_mu) {
    stats::dbinom(x = withheld_y, size = 1, prob = withheld_mu, log = TRUE)
  }
  family <- object$family$family
  if (family=='Gamma'){
    cv_data$cv_loglik <- ll_gamma(object, withheld_y, withheld_mu)
  } else if (family=='binomial') {
    cv_data$cv_loglik <- ll_binomial(object, withheld_y, withheld_mu)
  }

  # This gets the sum of the probabilities. I'm not sure why it is done the way
  # it is done. Doing sum(exp(x)) should get us the same answer. We just want to
  # transform back to non-log probs and then add them. Oh well, just know that
  # this gets us the sum of the probs.
  log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }
  # To get the average in ln space, we sum the probs and subtract the length.
  # Remember, in log space adding=mutliplying and subtracting=dividing.
  m1_elpd = log_sum_exp(cv_data$cv_loglik) - log(length(cv_data$cv_loglik))
  return(m1_elpd)
}



########################################
########################################
# Run

for (label in labels){

  print(paste('processing', label))

  file_t = glue('{base_dir}/sample_{label}_training.feather')
  file_h = glue('{base_dir}/sample_{label}_holdout.feather')

  df_t <- read_feather(file_t) # training data
  df_h <- read_feather(file_h) # holdout data

  # km required for mesh
  df_t$origin_x_km <- df_t$origin_x/1000
  df_t$origin_y_km <- df_t$origin_y/1000
  df_h$origin_x_km <- df_h$origin_x/1000
  df_h$origin_y_km <- df_h$origin_y/1000


  ###### Binomial model #######
  # For this model I ask: what are the chances of
  # a particles leaving its home area as predicted by
  # exposure, space, time?

  # create a 1/0 field to indicate if a particles settles at its origin
  df_t <- df_t %>% mutate(settle_home = case_when(distance==0 ~ 0, distance>0 ~ 1))
  df_h <- df_h %>% mutate(settle_home = case_when(distance==0 ~ 0, distance>0 ~ 1))

  # create barrier mesh
  bmesh_b <- bmesh(df_t, land, 5)

  # run sdmTMB model
  print('Running settle-at-home model')
  start_time <- Sys.time()
  m1 <- sdmTMB(
    settle_home ~ 1 + log(total_exposure), data = df_t,
    spde = bmesh_b, family = binomial(link='logit'),
    time = 'month',
    include_spatial = TRUE,
    fields = "IID"
  )
  i <- 1
  while(i<11){ # max out at 10 optimizations
    if(max(m1$gradients)>0.01){
      print("Running extra optimization")
      m1 <- sdmTMB::run_extra_optimization(m1, nlminb_loops = 1L, newton_loops = 1L)
    }
    i <- i+1
  }
  end_time <- Sys.time()
  runtime_b <- end_time - start_time
  print(paste('Binomial model runtime:', runtime_b))
  saveRDS(m1, glue("m_{label}_binom.rds"))

  # predict
  predictions_b <- predict_model(m1, df_t)

  # log-likelihood
  ll_b <- ll(m1, m1$data, 'settle_home')

  # Holdout data log-likelihood
  ll_b_ho <- ll(m1, df_h, 'settle_home')


  ###### Gamma model #######
  # For this model I exclude 0 distances and ask:
  # If a particle leaves its origin, how far does it travel,
  # as predicted by pld, exposure, space, and time?

  # dataset without 0 distances for gamma model
  df_tg <- filter(df_t, distance > 0)

  # create barrier mesh
  bmesh_g <- bmesh(df_tg, land, 5)

  # run sdmTMB model
  print('Running distance model')
  start_time <- Sys.time()
  m2 <- sdmTMB(
    distance ~ 1 + log(pld) + log(total_exposure),
    data = df_tg,
    spde = bmesh_g,
    family = Gamma(link='log'),
    time = 'month',
    include_spatial = TRUE,
    fields = "IID"
  )
  i <- 1
  while(i<11){ # max out at 10 optimizations
    if(max(m2$gradients)>0.01){
      print("Running extra optimization")
      m2 <- sdmTMB::run_extra_optimization(m2, nlminb_loops = 1L, newton_loops = 1L)
    }
    i <- i+1
  }
  end_time <- Sys.time()
  runtime_g <- end_time - start_time
  print(paste('Distance model runtime:', runtime_g))
  saveRDS(m2, glue("m_{label}_gamma.rds"))

  # predictions, residuals, r2
  predictions_g <- predict_model(m2, df_tg)
  resid_plots(m2, predictions_g, glue('g_{label}.jpg'))
  r2_g <- rsq(log(predictions_g$distance), predictions_g$est)

  # predicted vs.observed
  ggplot(predictions_g, aes(x=distance, y=exp(est)))+
    geom_abline(intercept=0, slope=1, color='blue') +
    geom_point() +
    scale_y_log10()+scale_x_log10()
  ggsave(glue('predobs_g_{label}.jpg'))

  # log-likelihood
  ll_g <- ll(m2, m2$data, 'distance')

  # Holdout data log-likelihood
  df_hog <- filter(df_h, distance > 0)
  ll_g_ho <- ll(m2, df_hog, 'distance')


  ### Compile df and write to csv ###

  m1_df <- tidy(m1)
  m1_df_rand <- tidy(m1, effects=c('ran_pars'))
  m2_df <- tidy(m2)
  m2_df_rand <- tidy(m2, effects=c('ran_pars'))

  # general fields
  samp_frac = strtoi(label) * (10^-(nchar(label)))
  # binomial model fields
  length_bi_tr = nrow(df_t)
  length_bi_ho = nrow(df_h)
  runtime_bi = runtime_b[[1]]
  intercept_bi = m1_df$estimate[m1_df$term=='(Intercept)']
  cof_exposure_bi = m1_df$estimate[m1_df$term=='log(total_exposure)']
  se_exposure_bi = m1_df$std.error[m1_df$term=='log(total_exposure)']
  mrp_bi = m1_df_rand$estimate[m1_df_rand$term=='range'][1]
  spat_sd_bi = m1_df_rand$estimate[m1_df_rand$term=='sigma_O']
  spatemp_sd_bi = m1_df_rand$estimate[m1_df_rand$term=='sigma_E']
  loglik_b_tr = ll_b
  loglik_b_ho = ll_b_ho
  # gamma model fields
  length_ga_tr = nrow(df_tg)
  length_ga_ho = nrow(df_hog)
  runtime_ga = runtime_g[[1]]
  intercept_ga = m2_df$estimate[m2_df$term=='(Intercept)']
  cof_pld_ga = m2_df$estimate[m2_df$term=='log(pld)']
  se_pld_ga = m2_df$std.error[m2_df$term=='log(pld)']
  cof_exposure_ga = m2_df$estimate[m2_df$term=='log(total_exposure)']
  se_exposure_ga = m2_df$std.error[m2_df$term=='log(total_exposure)']
  mrp_ga = m2_df_rand$estimate[m2_df_rand$term=='range'][1]
  disp_ga = m2_df_rand$estimate[m2_df_rand$term=='phi']
  spat_sd_ga = m2_df_rand$estimate[m2_df_rand$term=='sigma_O']
  spatemp_sd_ga = m2_df_rand$estimate[m2_df_rand$term=='sigma_E']
  loglik_ga_tr = ll_g
  loglik_ga_ho = ll_g_ho
  r2_ga = r2_g

  # df
  df_m1 <- data.frame(samp_frac, length_bi_tr, length_bi_ho, runtime_bi,
                      intercept_bi, cof_exposure_bi, se_exposure_bi, mrp_bi,
                      spat_sd_bi, spatemp_sd_bi, loglik_b_tr, loglik_b_ho,
                      length_ga_tr, length_ga_ho, runtime_ga, intercept_ga,
                      cof_pld_ga, se_pld_ga, cof_exposure_ga, se_exposure_ga,
                      mrp_ga, disp_ga, spat_sd_ga, spatemp_sd_ga, loglik_ga_tr,
                      loglik_ga_ho, r2_ga)

  # write to csv
  tbl_res <- 'tbl_res.csv'
  if (!file.exists(tbl_res)) {
    write.csv(df_m1, tbl_res, row.names=FALSE)
  } else {
    dfr <- read.csv(tbl_res)
    dfr <- rbind(dfr, df_m1)
    write.csv(dfr, tbl_res, row.names=FALSE)
  }

}


