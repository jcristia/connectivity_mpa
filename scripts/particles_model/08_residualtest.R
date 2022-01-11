
# following this to check residuals in a new way for glmms:
# https://github.com/pbs-assess/sdmTMB/blob/sim2/vignettes/residual-checking.Rmd


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
  geom_histogram() #+ facet_grid(rows=vars(group_id))
ggplot(df_t, aes(log(total_exposure))) +
  geom_histogram() +
  facet_grid(rows=vars(group_id)) # group by ocean grouping
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
plot(mesh)
bm <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                      proj_scaling = 1000, plot=TRUE, range_fraction=0.2)


# run sdmTMB model (ORIGINAL)
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

# Printing out model summary:
m2
# Matern range scale: it is in km and tells you the distance that two points can be
# considered the same
# For my coefficient estimates, I also get a standard error (se). I can tell if
# my coefficients are significant by subtracting se*2 from them. If this doesn't
# overlap with zero, then they are significant relationships.


########################################
########################################

# For my maps - do them for 1 PLD. He would do them for the median PLD. He had
# me change the code to filter for a PLD, so some of my plots are updated for
# this which cleaned them up nicely. PLD is independent of exposure (they don't
# interact, they are just additive), so PLD just raises and lowers my total
# values proportionally anyways. Therefore, just pick 1 pld value.


### Predict ###
predictions <- predict(m2)
#predictions <- predict(m2, m2$data %>% filter(pld==3))   # filtered for just 1 PLD
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

ggplot(predictions, aes(x = log(pld), y = residuals)) + geom_jitter() # residual issues at higher PLDs
ggplot(predictions, aes(x = log(total_exposure), y = residuals)) + geom_point() # this kinda looks ok






### New residual stuff ###
# https://pbs-assess.github.io/sdmTMB/articles/residual-checking.html

stan_fit <- tmbstan::tmbstan(m2$tmb_obj, iter = 50, warmup = 49, chains = 1) # SLOOOOW, get warnings
stan_eta <- predict(m2, tmbstan_model = stan_fit)
stan_mu <- as.numeric(m2$family$linkinv(stan_eta))
mcmc_res <- sdmTMB:::qres_gamma(m2, y = df_tg$distance, mu = stan_mu)
qqnorm(mcmc_res);qqline(mcmc_res)


s_gamma <- simulate(m2, nsim = 500)
dim(s_gamma)

pred_fixed <- m2$family$linkinv(predict(m2)$est_non_rf)
r_gamma <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = df_tg$distance,
  fittedPredictedResponse = pred_fixed
)
plot(r_gamma)
DHARMa::testResiduals(r_gamma)
DHARMa::testSpatialAutocorrelation(r_gamma, x = df_tg$origin_x_km, y = df_tg$origin_y_km)
# spat auto doesn't work because of duplicate coordinates. I couldn't figure out
# how to fix this. Moving on.

# So in interpreting all of this:
# The qqplot doesn't look terrible, but...
# The dispersion of expected data differs significantly from that of the observed,
# There are more outliers than expected.

# However, is it a problem?
# The qqplot actually looks pretty good.
# It's not going to change my conclusion - that we are missing some variables to help us predict.

# See "General remarks on interpreting residual patterns and tests"
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

# "There are many situations in statistics where it is common practice to work with “wrong models”. For example, many statistical models used shrinkage estimators, which purposefully bias parameter estimates to certain values. Random effects are a special case of this. If DHARMa residuals for these estimators are calculated, they will often show a slight pattern in the residuals even if the model is correctly specified, and tests for this can get significant for large sample sizes. "

#"Significance in hypothesis tests depends on at least 2 ingredients: strength of the signal, and the number of data points. If you have a lot of data points, residual diagnostics will nearly inevitably become significant, because having a perfectly fitting model is very unlikely. That, however, doesn’t necessarily mean that you need to change your model. "

# "Important conclusion: DHARMa only flags a difference between the observed and expected data - the user has to decide whether this difference is actually a problem for the analysis!"


# SO KEEP MOVING ON













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

# wouldn't we need to split this out by PLD?
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
  geom_sf(aes(color = (est_non_rf)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Fixed effects component")

# NOTE: in these next plots, Patrick had me remove the exp() as well as the trans
# log10 for the color.
# Make a version of each again so that I can understand why.

# spatial random effects that represent consistent deviations in space through
# time that are not accounted for by our fixed effects. In other words, these
# deviations represent consistent biotic and abiotic factors that are affecting
# distance but are not accounted for in the model
predictions %>%
  ggplot()+
  geom_sf(aes(color = (omega_s)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial component")
  #facet_wrap(~month)
# For the spatial map. Don't do by month. They are all the same anyways. It
# didn't make sense to do it this way to begin with.

predictions %>%
  ggplot()+
  geom_sf(aes(color = (epsilon_st)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatiotemporal component")+
  facet_wrap(~month)
# Interpreting: spatiotemporal: once you account for fixed effects and spatial
# random effects - this is what is left In January we see outer coast,
# disappears in May, then see it in the Salish Sea in August


###########################################
# Distance vs. Exposure, 1 plot per PLD, with prediction line overlaid

# ggplot(predictions, aes(x=log(total_exposure), y=log(distance))) +
#   geom_point() +
#   facet_grid(rows=vars(pld[1]))
#
# ggplot(predictions, aes(x=log(total_exposure), y=est)) +
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)+
#   facet_grid(rows=vars(pld))

ggplot(predictions) +
  geom_point(aes(x=log(total_exposure), y=log(distance))) +
  geom_smooth(aes(x=log(total_exposure), y=est), method='loess', se=TRUE) +
  facet_grid(rows=vars(pld))

ggplot(predictions) +
  geom_point(aes(x=log(total_exposure), y=est)) +
  geom_smooth(aes(x=log(total_exposure), y=est), method='loess', se=TRUE) +
  facet_grid(rows=vars(pld))

#Would this be a good way to show model results?

# He said instead, I should make a gradient of new exposure values and get
# predictions for these:
# Make a dataframe, could do it for 1 or a few PLDs then
# do: predict(m, gradient, se_fit = TRUE, re_form = NA), where gradient is my
# new dataframe this is spatial agnostic which is what I want - so no need to
# put in a spatial location (since those are random effects anyways - we're not
# trying to say that it is a predictor)

# Then, plot as a ribbon plot:
# geom_ribbon(aes(ymin = exp(est - est_se*2), ymax = exp(est + est_se*2)), color = NA, alpha = 0.35)
# notice... the se*2 is the standard error times 2.

# Lastly, he said I could also do a similar plot for just ONE exposure
# so I would look at a gradient of PLD
# this would require choosing a few common exposure values
# basically, I am taking slices through the model


###########################################
