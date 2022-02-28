# Building on the last script...
# Figure out how to calculate log likelihood without the sdmTMB cross validation function
# Compare additional model configurations (log pld, sqrt exposure) on just distance
# values greater than 0.

library(tidyverse)
library(sdmTMB)
library(arrow)
library(patchwork)
library(sf)
library(glue)


########################################
# Read
# All sampling and filtering done in Python

label <- '0003'

file_t = glue('~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all/sample_{label}_training.feather')
file_h = glue('~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all/sample_{label}_holdout.feather')
df_t <- read_feather(file_t) # training data
df_h <- read_feather(file_h) # holdout data


#########################################
# Create mesh and barrier mesh

df_t$origin_x_km <- df_t$origin_x/1000
df_t$origin_y_km <- df_t$origin_y/1000

land <- 'C:/Users/jcristia/Documents/GIS/MSc_Projects/MPA_connectivity/spatial/Coastline/landmask_FINAL.shp'
land_sf <- st_read(land)

mesh <- make_mesh(df_t, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bmesh <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=FALSE, range_fraction=0.2)


##########################################
# Model

start_time <- Sys.time()
m1 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_t,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
end_time <- Sys.time()
runtime <- end_time - start_time
if(max(m1$gradients)>0.01){
  m1 <- sdmTMB::run_extra_optimization(m1, nlminb_loops = 1L, newton_loops = 1L)
}

saveRDS(m1, glue("m1_{label}.rds"))
m1 <- readRDS(glue("m1_{label}.rds"))


#########################################
# Model analysis

summary(m1)

predictions1 <- predict(m1)
predictions1 <- st_as_sf(predictions1, coords = c("origin_x", "origin_y"))
predictions1 <- st_set_crs(predictions1, value = "epsg:3005")

# residuals
predictions1$residuals <- residuals(m1)
qqnorm(predictions1$residuals, ylim=c(-5,5));abline(a = 0, b = 1)
plot(density(predictions1$residuals))
ggplot(predictions1, aes(x = est, y = residuals)) + geom_point()

# predicted vs. observed
ggplot(predictions1, aes(distance)) +
  geom_histogram()
ggplot(predictions1, aes(exp(est))) +
  geom_histogram() # Notice how the 1m distances skew the distribution in comparison to the predicted data

ggplot(predictions1, aes(x=distance, y=exp(est)))+
  geom_point()
ggplot(predictions1, aes(x=distance, y=exp(est)))+
  geom_abline(intercept=0, slope=1, color='blue') +
  geom_point() +
  scale_y_log10()+scale_x_log10()

# Calculate r-squared
rsq <- function (x, y) cor(x, y) ^ 2
#rsq <- function(x, y) summary(lm(y~x))$r.squared
r2_m1 <- rsq(log(predictions1$distance), predictions1$est)

# Patrick thinks that maybe we should log pld in the model
ggplot(predictions1, aes(x=log(pld), y=log(distance)))+
  geom_jitter(width=0.1)
ggplot(predictions1, aes(x=pld, y=log(distance)))+
  geom_jitter(width=0.1)
# Similar to the Salish Sea, we see a slight trend up between 1 and 3 days, but not much after that

aic1 <- AIC(m1)



#########################################
# Figure out Log likelihood

# We are no longer relying on the sdmTMB_cv function to get log likelihood.
# Since we have so much holdout data there is no reason to do cross validation.
# So now we need to manually calculate log-likelihood instead of getting it
# from the cv function.
# https://github.com/pbs-assess/sdmTMB/blob/master/R/cross-val.R

# Get our observed and predicted data
object <- m1
cv_data <- m1$data
response <- 'distance'
predicted_obj <- predict(object, newdata = cv_data, return_tmb_obj = TRUE)
cv_data$cv_predicted <- object$family$linkinv(predicted_obj$data$est)
withheld_y <- predicted_obj$data[[response]]
withheld_mu <- cv_data$cv_predicted

# Function that established the gamma distribution from the modeled data and then
# get the probabilities of drawing the observed data from it (i.e. log likelihood)
ll_gamma <- function(object, withheld_y, withheld_mu) {
  .shape <- exp(object$model$par[["ln_phi"]])
  stats::dgamma(x = withheld_y, shape = .shape, scale = withheld_mu / .shape, log = TRUE)
}
cv_data$cv_loglik <- ll_gamma(object, withheld_y, withheld_mu)


# Average log likelihood probabilities

# This gets the sum of the probabilities.
# I'm not sure why it is done the way it is done.
# Doing sum(exp(x)) should get us the same answer. We just want to transform back
# to non-log probs and then add them.
# Oh well, just know that this gets us the sum of the probs.
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}
# To get the average in ln space, we sum the probs and subtract the length.
# Remember, in log space adding=mutliplying and subtracting=dividing.
m1_elpd = log_sum_exp(cv_data$cv_loglik) - log(length(cv_data$cv_loglik))





#########################################
# Compare additional model setups
# -log pld
# -sqrt exposure

# Since we will be doing a hurdle model that does not include distances of 0,
# do my comparisons with those removed.

df_samp <- filter(df_t, distance > 352)

mesh2 <- make_mesh(df_samp, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bmesh2 <- add_barrier_mesh(spde_obj = mesh2, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=FALSE, range_fraction=0.2)


# same as m1
m2 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_samp,
  spde = bmesh2, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
m2 <- sdmTMB::run_extra_optimization(m2, nlminb_loops = 1L, newton_loops = 1L)

# log pld
m3 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure), data = df_samp,
  spde = bmesh2, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m3$gradients)>0.01){
  m3 <- sdmTMB::run_extra_optimization(m3, nlminb_loops = 1L, newton_loops = 1L)
}

# sqrt exposure
m4 <- sdmTMB(
  distance ~ 1 + pld + sqrt(total_exposure), data = df_samp,
  spde = bmesh2, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m4$gradients)>0.01){
  m4 <- sdmTMB::run_extra_optimization(m4, nlminb_loops = 1L, newton_loops = 1L)
} # can't get this one to optimize, moving on

#log10 exposure
m5 <- sdmTMB(
  distance ~ 1 + log(pld) + log10(total_exposure), data = df_samp,
  spde = bmesh2, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m5$gradients)>0.01){
  m5 <- sdmTMB::run_extra_optimization(m5, nlminb_loops = 1L, newton_loops = 1L)
}


AIC(m2, m3, m4, m5)
# logging pld and exposure is the best one
# log 10 makes no difference

predictions2 <- predict(m2)
predictions3 <- predict(m3)
predictions4 <- predict(m4)
predictions5 <- predict(m5)

# residuals
predictions2$residuals <- residuals(m2)
predictions3$residuals <- residuals(m3)
predictions4$residuals <- residuals(m4)
predictions5$residuals <- residuals(m5)

qqnorm(predictions2$residuals, ylim=c(-5,5));abline(a = 0, b = 1)
qqnorm(predictions3$residuals, ylim=c(-5,5));abline(a = 0, b = 1)
qqnorm(predictions4$residuals, ylim=c(-5,5));abline(a = 0, b = 1)
qqnorm(predictions5$residuals, ylim=c(-5,5));abline(a = 0, b = 1)
# they are all still very similar with the skew at high values
# hopefully more data takes care of this
# Otherwise, using a gamma dist and log link may be wrong(???)

# look at a few other things:
ggplot(predictions2, aes(x = est, y = residuals)) + geom_point()
ggplot(predictions3, aes(x = est, y = residuals)) + geom_point()
plot(density(predictions2$residuals))
plot(density(predictions3$residuals))
test <- sample(predictions3$residuals, 5000)
shapiro.test(test) # test of normality. Above 0.05 is normal. Ours is definitely not.

ggplot(predictions3, aes(est))+
  geom_histogram()
ggplot(predictions3, aes(log(distance)))+
  geom_histogram()



# Try lognormal distribution
m6 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure), data = df_samp,
  spde = bmesh2, family = lognormal(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m6$gradients)>0.01){
  m6 <- sdmTMB::run_extra_optimization(m6, nlminb_loops = 1L, newton_loops = 1L)
}
predictions6 <- predict(m6)
predictions6$residuals <- residuals(m6) # I get an error on this step.
qqnorm(predictions6$residuals, ylim=c(-5,5));abline(a = 0, b = 1)

# try one more time to get tweedie to run
m7 <- sdmTMB(
  distance ~ 1, data = df_samp,
  spde = bmesh2, family = tweedie(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
) # NOPE

# try to get the version in the vignette to run
library(sdmTMB)
pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 12)
plot(pcod_spde)

m8 <- sdmTMB(density ~ 1, data = pcod,
             spde = pcod_spde, family = tweedie(link = "log"),
             spatial_trend = TRUE, time = "year",
             spatial_only = TRUE)
# runs super fast

# try with a basic mesh (no barrier) and less data
label <- '0001'
file_t = glue('~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all/sample_{label}_training.feather')
df_t <- read_feather(file_t) # training data
df_t$origin_x_km <- df_t$origin_x/1000
df_t$origin_y_km <- df_t$origin_y/1000
df_samp <- filter(df_t, distance > 352)
mesh <- make_mesh(df_samp, c('origin_x_km', 'origin_y_km'), cutoff = 15)
m9 <- sdmTMB(distance ~ 1, data = df_samp,
             spde = mesh, family = tweedie(link = "log"),
             spatial_trend = TRUE, time = "month",
             spatial_only = TRUE)
# still doesn't work


# review my distributions
ggplot(df_samp, aes(log(distance))) +
  geom_histogram()
ggplot(df_samp, aes(distance)) +
  geom_histogram()
ggplot(df_samp, aes(log(total_exposure))) +
  geom_histogram()
# It looks like log(total_exposure) doesn't actually follow a gamma distribution

# try a weird boxcox thing
bc_model <- lm(distance~total_exposure, data=df_samp)
library(MASS)
boxcox(bc_model, plotit=TRUE)
ggplot(df_samp, aes(((distance^-0.35)-1)/-0.35)) +
  geom_histogram()

