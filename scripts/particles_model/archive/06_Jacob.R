
# test script to show Jacob

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
label <- '0005'
# feather files directory
base_dir <- '~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all'
# coastline shapefile
land <- 'C:/Users/jcristia/Documents/GIS/MSc_Projects/MPA_connectivity/spatial/Coastline/landmask_FINAL.shp'



########################################
########################################

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


### data ###
file_t = glue('{base_dir}/sample_{label}_training.feather')
file_h = glue('{base_dir}/sample_{label}_holdout.feather')
df_t <- read_feather(file_t) # training data
df_h <- read_feather(file_h) # holdout data
df_t$origin_x_km <- df_t$origin_x/1000 # km required for mesh
df_t$origin_y_km <- df_t$origin_y/1000
df_h$origin_x_km <- df_h$origin_x/1000
df_h$origin_y_km <- df_h$origin_y/1000

#### exploratory plots ####
ggplot(df_t, aes(log(distance))) +
  geom_histogram()
ggplot(df_t, aes(log(total_exposure))) +
  geom_histogram()
ggplot(df_t, aes(total_exposure, distance)) +
  geom_point() +
  scale_y_log10()+scale_x_log10()
ggplot(df_t, aes(pld, log(distance))) +
  geom_jitter()

# datasets without 0 distances for gamma model
df_tg <- filter(df_t, (distance > 2500 ))
df_hog <- filter(df_h, (distance > 2500 ))

#### exploratory plots ####
ggplot(df_tg, aes(log(distance))) +
  geom_histogram()
ggplot(df_tg, aes(log(total_exposure))) +
  geom_histogram()
ggplot(df_tg, aes(total_exposure, distance)) +
  geom_point() +
  scale_y_log10()+scale_x_log10()


########################################
########################################



### Create mesh and barrier mesh ###
land_sf <- st_read(land, quiet=TRUE)
mesh <- make_mesh(df_tg, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bm <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                      proj_scaling = 1000, plot=TRUE, range_fraction=0.2)


# run sdmTMB model
print('Running distance model')
start_time <- Sys.time
m2 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)


# TO TRY: (based on notes from Jacob meeting)
# see notes in Evernote for more info
m2 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + (1|mpa_part_id_orig) + mpa_area,
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
# first try without space and time and just mpa_part_id (and pld, no exposure)
# then with area
# then also with mpa_area^2
# then with space and time and just mpa_part_id
# then again with area
# then also with mpa_area^2

# then also try with 3+ ocean groupings (Salish sea, outer coast, ...)
# first see how exposure and distance values cluster based on those groupings (overlay the distributions)

# From Patrick: try with splines (gam), which is now implemented.
# Sean has a brief example of this.
m2 <- sdmTMB(
  log(distance) ~ 1 + pld + s(total_exposure),
  data = df_tg,
  spde = bm,
  #family = Gamma(link='log'),
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

m2 <- readRDS(glue("m_{label}_gamma.rds"))


########################################
########################################


### Predict ###
predictions <- predict(m2, m2$data)
predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
predictions <- st_set_crs(predictions, value = "epsg:3005")


### Residuals ###
filename <- glue('g_{label}.jpg')
predictions$residuals <- residuals(m2)
ggplot(predictions, aes(sample=residuals)) + stat_qq_line() + stat_qq()
ggsave(paste('qqnorm_', filename, sep=''))
ggplot(predictions, aes(residuals)) + geom_density()
ggsave(paste('pdf_', filename, sep=''))
ggplot(predictions, aes(x = log(distance), y = residuals)) + geom_point()
ggsave(paste('resid_', filename, sep=''))


### R-squared ###
rsq <- function (x, y) cor(x, y) ^ 2
r2_g <- rsq(log(predictions$distance), predictions$est)


# predicted vs.observed
ggplot(predictions, aes(x=distance, y=exp(est)))+
  geom_abline(intercept=0, slope=1, color='blue') +
  geom_point() +
  scale_y_log10()+scale_x_log10()
ggsave(glue('predobs_g_{label}.jpg'))


# log-likelihood
ll_g <- ll(m2, m2$data, 'distance')


# Holdout data log-likelihood
ll_g_ho <- ll(m2, df_hog, 'distance')


m2
ll_g
ll_g_ho
r2_g


##########################################
##########################################
# fixed/random effects plots


predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(est)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Overall predicted")+
  facet_wrap(~month)

predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(est_non_rf)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Fixed effects component")+
  facet_wrap(~month)

# spatial random effects that represent consistent deviations in space through
# time that are not accounted for by our fixed effects. In other words, these
# deviations represent consistent biotic and abiotic factors that are affecting
# distance but are not accounted for in the model
predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(omega_s)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial component")+
  facet_wrap(~month)

predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(epsilon_st)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatiotemporal component")+
  facet_wrap(~month)




