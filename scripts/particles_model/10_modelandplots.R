
# Summarize model selection
# Predictions
# Plots


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
label <- '03'
# feather files directory
base_dir <- '~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all'
# coastline shapefile
land <- 'C:/Users/jcristia/Documents/GIS/MSc_Projects/MPA_connectivity/spatial/Coastline/landmask_FINAL.shp'



########################################
########################################
# Data

file_t = glue('{base_dir}/sample_{label}_training.feather')
file_h = glue('{base_dir}/sample_{label}_holdout.feather')
df_t <- read_feather(file_t) # training data
df_h <- read_feather(file_h) # holdout data
df_t$origin_x_km <- df_t$origin_x/1000 # km required for mesh
df_t$origin_y_km <- df_t$origin_y/1000
df_h$origin_x_km <- df_h$origin_x/1000
df_h$origin_y_km <- df_h$origin_y/1000

# exploratory plots
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

# exploratory plots
ggplot(df_tg, aes(log(distance))) +
  geom_histogram()
ggplot(df_tg, aes(log(total_exposure))) +
  geom_histogram()
ggplot(df_tg, aes(total_exposure, distance)) +
  geom_point() +
  scale_y_log10()+scale_x_log10()



########################################
########################################
# Model

# Create mesh and barrier mesh
land_sf <- st_read(land, quiet=TRUE)
mesh <- make_mesh(df_tg, c('origin_x_km', 'origin_y_km'), cutoff = 5)
#plot(mesh)
bm <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                       proj_scaling = 1000, plot=FALSE, range_fraction=0.2)


# Summary of model selection from previous scripts:
# Given the structure of my question, I already knew that my model will include
# PLD and total_exposure. Therefore, varying the model only involved:
# - log(), sqrt(), ()**2, or no adjustment on pld and exposure
# - with and without space and time as sdmtmb random effects
# - with and without mpa_part_id_orig as a random effect

# See the 'archive' folder for the history of this testing. Model selection is
# not a focus of this chapter since I'm not trying to select the 'best' model to
# predict dispersal distance. Instead, I want to answer "Is there a
# relationship between dispersal distance and pld/exposure?".
# Therefore, I'm not including the history of model selection in this final
# script.


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

# model summary
m2
tidy(m2, conf.int = TRUE)
# Matern range scale: it is in km and tells you the distance that two points can
# be considered the same.
# Standard error (se): +/- se*2 from coefficients. Significant relationship if
# it doesn't overlap with zero.

# Get the spatial and spatiotemporal SD values and the
# confidence intervals.
m2$sd_report
tidy(m2, 'ran_pars', conf.int=TRUE)


########################################
########################################
# Evaluate model

### Predict ###
predictions <- predict(m2)
predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
predictions <- st_set_crs(predictions, value = "epsg:3005")

### Residuals - traditional ###
# need to run this with full model predictions
filename <- glue('g_{label}.jpg')
predictions$residuals <- residuals(m2)
ggplot(predictions, aes(sample=residuals)) + stat_qq_line() + stat_qq()
ggsave(paste('qqnorm_', filename, sep=''))
ggplot(predictions, aes(residuals)) + geom_density()
ggsave(paste('pdf_', filename, sep=''))
ggplot(predictions, aes(x = log(distance), y = residuals)) + geom_point()
ggsave(paste('resid_', filename, sep=''))

ggplot(predictions, aes(x = log(pld), y = residuals)) + geom_jitter() # residual issues at higher PLDs
ggsave(paste('residpld_', filename, sep=''))
ggplot(predictions, aes(x = log(total_exposure), y = residuals)) + geom_point() # this kinda looks ok,
# and it confirms that all my variability is happening when there are higher
# exposure values, which kinda makes sense: when low exposure, things only go
# a short distance, but when high exposure things go short and long distances
# and it is hard to predict.
ggsave(paste('residexposure_', filename, sep=''))

### Residuals - GLMMs ###
# https://pbs-assess.github.io/sdmTMB/articles/residual-checking.html
s_gamma <- simulate(m2, nsim = 500)
dim(s_gamma)
pred_fixed <- m2$family$linkinv(predict(m2)$est_non_rf)
r_gamma <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = df_tg$distance,
  fittedPredictedResponse = pred_fixed)
plot(r_gamma, quantreg=FALSE)
DHARMa::testResiduals(r_gamma)
DHARMa::plotQQunif(r_gamma, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)


# Interpreting:
# The qqplot doesn't look terrible, but...
# The dispersion of expected data differs significantly from that of the observed,
# there are more outliers than expected.
# However, is it a problem?
# The qqplot actually looks pretty good.
# Either way, it's not going to change my conclusion - that the model isn't
# adequate to predict, but we do see a relationship.

# See "General remarks on interpreting residual patterns and tests"
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# "There are many situations in statistics where it is common practice to work with “wrong models”. For example, many statistical models used shrinkage estimators, which purposefully bias parameter estimates to certain values. Random effects are a special case of this. If DHARMa residuals for these estimators are calculated, they will often show a slight pattern in the residuals even if the model is correctly specified, and tests for this can get significant for large sample sizes. "
#"Significance in hypothesis tests depends on at least 2 ingredients: strength of the signal, and the number of data points. If you have a lot of data points, residual diagnostics will nearly inevitably become significant, because having a perfectly fitting model is very unlikely. That, however, doesn’t necessarily mean that you need to change your model. "
# "Important conclusion: DHARMa only flags a difference between the observed and expected data - the user has to decide whether this difference is actually a problem for the analysis!"

### R-squared ###
rsq <- function (x, y) cor(x, y) ^ 2
r2_g <- rsq(log(predictions$distance), predictions$est)

# predicted vs.observed
ggplot(predictions, aes(x=distance, y=exp(est)))+
  geom_abline(intercept=0, slope=1, color='blue') +
  geom_point() +
  scale_y_log10()+scale_x_log10()
ggsave(glue('predvsobs_g_{label}.jpg'))



########################################
########################################
# log-likelihood

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


ll_g <- ll(m2, m2$data, 'distance')
ll_g_ho <- ll(m2, df_hog, 'distance')



##########################################
##########################################
# Maps of predictions (fixed, spatial, and spatiotemporal effects)


# Plot predictions for just 1 PLD (median values).
predictions <- predict(m2, m2$data %>% filter(pld==10))   # filtered for just 1 PLD
predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
predictions <- st_set_crs(predictions, value = "epsg:3005")
# PLD is independent of exposure (they don't interact, they are just additive),
# so PLD just raises and lowers my total values proportionally anyways.
# Therefore, just show results for 1 pld value.
# This also makes the points on the map easier to see and understand.

# Fixed effects + all random effects
predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(est)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Prediction (fixed effects + all random effects)",
          subtitle = 'Dispersal distance traveled from origin location') +
  facet_wrap(~month)

# Fixed effects (i.e. just the effect of exposure on distance)
# No month included, since that is a random effect.
predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(est_non_rf)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Prediction (fixed effects only")

# Spatial random effects only
# "spatial random effects that represent consistent deviations in space through
# time that are not accounted for by our fixed effects. In other words, these
# deviations represent consistent biotic and abiotic factors that are affecting
# distance but are not accounted for in the model"
# omega_s - indicates a relative magnitude and direction in terms of the
# deviation from the fixed effect prediction?
predictions %>%
  ggplot()+
  geom_sf(aes(color = (omega_s)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial random effects only")
# so it seems like on the outer coast we get better predictions with exposure,
# but in the Salish Sea there is a lot of variation unaccounted for.

# Spatiotemporal random effects
# "random effects that represent deviation from the fixed effect predictions and
# the spatial random effect deviations. These represent biotic and abiotic
# factors that are changing through time and are not accounted for in the model."
predictions %>%
  ggplot()+
  geom_sf(aes(color = (epsilon_st)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatiotemporal random effects only")+
  facet_wrap(~month)
# Interpreting: In January we see outer coast have higher deviations, this
# disappears in May, then we see it higher in the Salish Sea in August.


# Export model predictions so that I can make these maps pretty in ArcGIS
pred_out <- predictions %>% st_drop_geometry() # you get errors in arc if this column is kept in
pred_out['origin_x'] <- pred_out$origin_x_km * 1000
pred_out['origin_y'] <- pred_out$origin_y_km * 1000
pred_out['est_km'] <- exp(pred_out$est) / 1000
pred_out['est_non_rf_km'] <- exp(pred_out$est_non_rf) / 1000
# to get omega and epsilon in km, you can't just exp() the value.
# I'm calculating these assuming that:
# omega_s contribution is applied after fixed effects
# then epsilon_st is the remaining difference.
# I'm not sure about this and need to look into it further.
# Patrick doesn't think that this is right. Sean does plot just exp(est_non_rf)
# in one of his examples, but Patrick thinks that this was just for illustrative
# purposes and that you can't convert this to any meaningful values.
#pred_out['omegas_km'] <- exp(pred_out$est_non_rf + pred_out$omega_s)/1000 - pred_out$est_non_rf_km
#pred_out['epsilonst_km'] <- pred_out$est_km - pred_out$omegas_km - pred_out$est_non_rf_km
write.csv(pred_out, 'predictions.csv')



##########################################
##########################################
# Plots of model predictions
# Distance vs. Exposure
# Predict over a gradient of exposure values and then plot to show the predicted
# relationship between distance and exposure.

# build dataframe
# Must include: exposure, pld, mpa_part_id_orig
# Must also include a time value.
min_exposure = log(min(df_tg$total_exposure))
max_exposure = log(max(df_tg$total_exposure))
gradient <- seq(min_exposure, max_exposure,length.out=25)
gradient <- exp(gradient)
mpa_unique <- sample(unique(df_tg$mpa_part_id_orig),5) # this doesn't end up
# mattering for the way I make the plot
pld_vals <- c(1,3,7,10,21,30,40,60)
month <- 5L
# expand.grid creates a dataframe of all unique combinations
predictor_df <- expand.grid(
  mpa_part_id_orig=mpa_unique,
  pld=pld_vals,
  total_exposure=gradient,
  month=month
)

# Predict
# re_form=NA excludes spatial/spatiotemporal random effects, and therefore, we
# don't need to include any coordinates.
# re_form_iid=NA gives us population level predictions (so it excludes
# mpa_origin_id from the predictions).
p1 <- predict(m2, predictor_df, se_fit=TRUE, re_form=NA, re_form_iid=NA)

# test it without re_form_iid. This will give us individual results, but how
# would I plot this? I suppose by excluding mpa id I am missing out on showing
# a lot of the variation between mpas.
#p2 <- predict(m2, predictor_df, se_fit=TRUE, re_form=NA)

# print out p1 and p2 to see the difference.


# Plot distance vs. exposure by PLD
p1$pld <- as.factor(p1$pld)
ggplot(p1, aes(log(total_exposure), exp(est),
               ymin = exp(est - 1.96 * est_se),
               ymax = exp(est + 1.96 * est_se),
               group = pld
)) +
  geom_line(aes(colour = pld), lwd = 1) +
  geom_ribbon(aes(fill = pld), alpha = 0.1) +
  scale_colour_viridis_d(name='PLD (days)', guide = guide_legend(reverse=TRUE)) +
  scale_fill_viridis_d(name='PLD (days)', guide = guide_legend(reverse=TRUE)) +
  coord_cartesian(expand = F) +
  labs(x = "Exposure at origin site  log(m)", y = "Predicted distance (m)")+
  theme_classic() +
  theme(aspect.ratio=1/1.5)
ggsave('fig07_preddistanceVexposure.jpg')





###########################################
