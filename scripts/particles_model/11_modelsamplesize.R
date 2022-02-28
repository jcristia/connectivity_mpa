
# Increase dataset sample size and evaluate with average cv log likelihood to
# determine minimal sample size needed.
# Increase mesh size to gauge sensitivity.


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
#labels <- list('00001', '00005', '0001', '0005', '001', '005', '01', '02', '03', '04, '05', '06', '07')
labels <- list('07')
# feather files directory
base_dir <- '~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all'
# coastline shapefile
land <- 'C:/Users/jcristia/Documents/GIS/MSc_Projects/MPA_connectivity/spatial/Coastline/landmask_FINAL.shp'


########################################
########################################
# Average log-likelihood function

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

  ########################################
  ########################################
  # Data

  file_t = glue('{base_dir}/sample_{label}_training.feather')
  file_h = glue('{base_dir}/sample_{label}_holdout.feather')

  df_t <- read_feather(file_t) # training data
  df_h <- read_feather(file_h) # holdout data

  # km required for mesh
  df_t$origin_x_km <- df_t$origin_x/1000
  df_t$origin_y_km <- df_t$origin_y/1000
  df_h$origin_x_km <- df_h$origin_x/1000
  df_h$origin_y_km <- df_h$origin_y/1000

  # datasets without 0 distances for gamma model
  df_tg <- filter(df_t, (distance > 2500 ))
  df_hog <- filter(df_h, (distance > 2500 ))


  ########################################
  ########################################
  # Model

  # Create mesh and barrier mesh
  land_sf <- st_read(land, quiet=TRUE)
  mesh <- make_mesh(df_tg, c('origin_x_km', 'origin_y_km'), cutoff = 5)
  #plot(mesh)
  bm <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                         proj_scaling = 1000, plot=FALSE, range_fraction=0.2)

  # sdmTMB model
  print('Running distance model')
  start_time <- Sys.time()
  df_tg$mpa_part_id_orig <- factor(df_tg$mpa_part_id_orig)
  m2 <- sdmTMB(
    distance ~ 1 + log(pld) + log(total_exposure) + (1|mpa_part_id_orig),
    data = df_tg,
    spde = bm,
    family = Gamma(link='log'),
    time = 'month',
    include_spatial = TRUE,
    fields = "IID"
  )
  i <- 1
  while(i<30){ # max out at 30 optimizations
    if(max(m2$gradients)>0.01){
      print("Running extra optimization")
      m2 <- sdmTMB::run_extra_optimization(m2, nlminb_loops = 1L, newton_loops = 1L)
      print(max(m2$gradients))
    }
    i <- i+1
  }
  end_time <- Sys.time()
  runtime_g <- (end_time - start_time)[[1]]
  print(paste('Distance model runtime:', runtime_g))
  saveRDS(m2, glue("m_{label}_gamma.rds"))
  m2 <- readRDS(glue("m_{label}_gamma.rds"))


  ########################################
  ########################################
  # Evaluate model

  # Predict
  #predictions <- predict(m2)
  predictions <- predict(m2, m2$data %>% filter(pld==10))   # filtered for just 1 PLD
  predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
  predictions <- st_set_crs(predictions, value = "epsg:3005")

  # R-squared
  rsq <- function (x, y) cor(x, y) ^ 2
  r2_g <- rsq(log(predictions$distance), predictions$est)

  # average log likelihood
  ll_g <- ll(m2, m2$data, 'distance')
  df_hog$mpa_part_id_orig <- factor(df_hog$mpa_part_id_orig)
  ll_g_ho <- ll(m2, df_hog, 'distance')


  ########################################
  ########################################
  ### Compile df and write to csv ###

  m2_df <- tidy(m2)
  m2_df_rand <- tidy(m2, effects=c('ran_pars'))

  # general fields
  samp_frac = strtoi(label) * (10^-(nchar(label)))
  # gamma model fields
  length_ga_tr = nrow(df_tg)
  length_ga_ho = nrow(df_hog)
  runtime_ga = runtime_g
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
  df_m1 <- data.frame(samp_frac,
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
###########################################

# interpreting:

# regardless of where the metrics "settle", anything below the 0.02 fraction
# results in some MPAs only having 1 data point per time per pld (84 * 0.02).
# Going below this level results in MPAs that have more points than 84 released
# not being selected because the fraction is too small.
# Therefore, I need to at least use this level, but probably higher.
# 0.03 looks ok, it dips lowest there and then back up for 0.04.


###########################################



########################################
########################################
# Test sensitivity of grid resolution

# apparently a model can IMPROVE with increasing grid resolution. I suspect that
# this is only if you started with a very small resolution, but I will test it
# to be sure.

label <- '03'
file_t = glue('{base_dir}/sample_{label}_training.feather')
file_h = glue('{base_dir}/sample_{label}_holdout.feather')
df_t <- read_feather(file_t) # training data
df_h <- read_feather(file_h) # holdout data
df_t$origin_x_km <- df_t$origin_x/1000
df_t$origin_y_km <- df_t$origin_y/1000
df_h$origin_x_km <- df_h$origin_x/1000
df_h$origin_y_km <- df_h$origin_y/1000
# datasets without 0 distances for gamma model
df_tg <- filter(df_t, (distance > 2500 ))
df_hog <- filter(df_h, (distance > 2500 ))

# Model
cutoffs <- list(5, 7, 10, 20)
for (cutoff in cutoffs){

  print(paste('processing cutoff', toString(cutoff)))

  # Create mesh and barrier mesh
  land_sf <- st_read(land, quiet=TRUE)
  mesh <- make_mesh(df_tg, c('origin_x_km', 'origin_y_km'), cutoff = cutoff)
  bm <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                         proj_scaling = 1000, plot=FALSE, range_fraction=0.2)

  start_time <- Sys.time()
  df_tg$mpa_part_id_orig <- factor(df_tg$mpa_part_id_orig)
  m2 <- sdmTMB(
    distance ~ 1 + log(pld) + log(total_exposure) + (1|mpa_part_id_orig),
    data = df_tg,
    spde = bm,
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
  runtime_g <- (end_time - start_time)[[1]]
  print(paste('Distance model runtime:', runtime_g))
  cutoff_str <- toString(cutoff)
  saveRDS(m2, glue("m_{label}_gamma_grid{cutoff_str}.rds"))
  m2 <- readRDS(glue("m_{label}_gamma_grid{cutoff_str}.rds"))

  # Predict
  #predictions <- predict(m2)
  predictions <- predict(m2, m2$data %>% filter(pld==10))   # filtered for just 1 PLD
  predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
  predictions <- st_set_crs(predictions, value = "epsg:3005")

  # R-squared
  rsq <- function (x, y) cor(x, y) ^ 2
  r2_g <- rsq(log(predictions$distance), predictions$est)

  # average log likelihood
  ll_g <- ll(m2, m2$data, 'distance')
  df_hog$mpa_part_id_orig <- factor(df_hog$mpa_part_id_orig)
  ll_g_ho <- ll(m2, df_hog, 'distance')

  ### Compile df and write to csv ###

  m2_df <- tidy(m2)
  m2_df_rand <- tidy(m2, effects=c('ran_pars'))

  # general fields
  samp_frac = strtoi(label) * (10^-(nchar(label)))
  grid_res = cutoff
  # gamma model fields
  length_ga_tr = nrow(df_tg)
  length_ga_ho = nrow(df_hog)
  runtime_ga = runtime_g
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
  df_m1 <- data.frame(samp_frac, cutoff,
                      length_ga_tr, length_ga_ho, runtime_ga, intercept_ga,
                      cof_pld_ga, se_pld_ga, cof_exposure_ga, se_exposure_ga,
                      mrp_ga, disp_ga, spat_sd_ga, spatemp_sd_ga, loglik_ga_tr,
                      loglik_ga_ho, r2_ga)

  # write to csv
  tbl_res <- 'tbl_res_gridresolution.csv'
  if (!file.exists(tbl_res)) {
    write.csv(df_m1, tbl_res, row.names=FALSE)
  } else {
    dfr <- read.csv(tbl_res)
    dfr <- rbind(dfr, df_m1)
    write.csv(dfr, tbl_res, row.names=FALSE)
  }

}
###########################################

# interpreting:
# https://peerj.com/articles/12783/
# Too fine of a spatial resolution could result in an underrepresentation of perceived uncertainty unless one can visualize uncertainty in the prediction at each location, a feat that is highly computationally demanding in predictive process SDMs. Too coarse of a resolution may increase uncertainty in derived products like total biomass estimates as this involves integrating over an increasing extent of modeled or unmodeled habitat covariates. A general rule of thumb for assessing an initial resolution is that it should typically be no finer than the scale of the sampling unit of the survey design. A rationale for making predictions at coarser scales may include trying to link the response of an SDM to environmental covariates predicted from global climate models, which often have a minimum resolution of ~10 km (Porfirio et al., 2014), to determine support for climate drivers of past changes or project future distribution changes.

# 5 km makes sense for me (just slightly bigger than scale of NEP36 model)
# I don't expect too many significant changes to ocean currents over this scale.
# It should capture most of the variation in how things travel and deviate.
# Also, I had trouble getting things to run at a finer scale anyways.


###########################################








