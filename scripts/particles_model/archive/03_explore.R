# Building on the last script...
# Do cross validation and log likelihood
# R2
# Compile table to compare model runs
# Look at effect of exposure in just a basic GLM compared to sdmTMB glmm
# Try additional model setups (tweedie, hurdle, set 0 distance values differently)

library(tidyverse)
library(sdmTMB)
library(arrow)
library(patchwork)
library(sf)
library(glue)


########################################
# Read
# All sampling and filtering done in Python

label <- '0001'

file_t = glue('~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all/sample_{label}_training.feather')
file_h = glue('~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_df/output_all/sample_{label}_holdout.feather')
df_t <- read_feather(file_t) # training data
df_h <- read_feather(file_h) # holdout data


#########################################
# Create mesh and barrier mesh

df_t$origin_x_km <- df_t$origin_x/1000
df_t$origin_y_km <- df_t$origin_y/1000
mesh <- make_mesh(df_t, c('origin_x_km', 'origin_y_km'), cutoff = 5)
#plot(mesh)
land <- 'C:/Users/jcristia/Documents/GIS/MSc_Projects/MPA_connectivity/spatial/Coastline/landmask_FINAL.shp'
land_sf <- st_read(land)
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



#########################################
# Cross validation and log likelihood

start_time <- Sys.time()
m1_cv <- sdmTMB_cv(
  distance ~ 1 + pld + log(total_exposure), data = df_t,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID",
  k_folds = 5,
  #control = sdmTMBcontrol(nlminb_loop=3L, newton_loops = 3L)
)
end_time <- Sys.time()
m1_cv$fold_loglik
for (g in 1:length(m1_cv$max_gradients)) {
  if(max(m1_cv$models[[g]]$gradients)>0.01){
    print(g)# need to reference gradients within individual model
    m1_cv$models[[g]] <- sdmTMB::run_extra_optimization(m1_cv$models[[g]], nlminb_loops = 1L, newton_loops = 1L)
  }
}
# NOT SURE WHAT TO DO HERE
# How do we get new fold log likelihood values once we have run extra optimizations?

saveRDS(m1_cv, glue("m1_cv_{label}.rds"))
m1_cv <- readRDS(glue("m1_cv_{label}.rds"))

summary(m1_cv)
m1_cv$fold_loglik
m1_cv$sum_loglik

# Calculate coefficient of variance of log likelihood values
# (standard deviation / mean)
cv <- sd(m1_cv$fold_loglik) / mean(m1_cv$fold_loglik) * 100


#########################################
# Predict for holdout data

df_h$origin_x_km <- df_h$origin_x/1000
df_h$origin_y_km <- df_h$origin_y/1000
pred_ho <- predict(m1, newdata=df_h)

ggplot(pred_ho, aes(x=distance, y=exp(est)))+
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  scale_y_log10()+scale_x_log10()

rsq <- function (x, y) cor(x, y) ^ 2
r2_ho <- rsq(log(pred_ho$distance), pred_ho$est)
# compare to r2_m1?
# make sure the coefficients are the same?


#########################################
# Output plots

predictions1 %>%
  ggplot()+
  geom_sf(aes(color = exp(est)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Overall predicted")+
  facet_wrap(~month)

predictions1 %>%
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
predictions1 %>%
  ggplot()+
  geom_sf(aes(color = exp(omega_s)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial component")+
  facet_wrap(~month)

predictions1 %>%
  ggplot()+
  geom_sf(aes(color = exp(epsilon_st)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatiotemporal component")+
  facet_wrap(~month)


##########################################
# Exposure
# Compare to simple GLM and get a sense of how much influence it has.

m7 <- glm(distance ~ 1 + log(total_exposure), data=df_t, family= Gamma(link='log'))
summary(m7)
df_t$pred7 <- predict(m7)
ggplot(df_t, aes(distance, exp(pred7))) +
  geom_point() +
  scale_y_log10()+scale_x_log10()
ggplot(predictions1, aes(distance, exp(est))) +
  geom_point() +
  scale_y_log10()+scale_x_log10()
# so basically, exposure has a minimal effect. There is a bit of a trend for
# particles that traveled more than 1000m though.
# This make sense given the resolution of the hydro model and the resolution
# of exposure and particle distance spatial analysis.

# From Patrick: the coefficients aren't that different though, which is a good
# thing because it means there isn't some weird combination of exposure with
# space(?)


#########################################
# compile table to compare between model runs with increasing sample size

m1_df <- tidy(m1)
m1_df_rand <- tidy(m1, effects=c('ran_pars'))

# compile df
label_col = label
len_t = nrow(df_t)
len_h = nrow(df_h)
rt = runtime[[1]]  # might have to enter this in manually
intercept = m1_df$estimate[m1_df$term=='(Intercept)']
cof_pld = m1_df$estimate[m1_df$term=='pld']
se_pld = m1_df$std.error[m1_df$term=='pld']
cof_exp = m1_df$estimate[m1_df$term=='log(total_exposure)']
se_exp = m1_df$std.error[m1_df$term=='log(total_exposure)']
mrp = m1_df_rand$estimate[m1_df_rand$term=='range'][1]
disp = m1_df_rand$estimate[m1_df_rand$term=='phi']
spat_sd = m1_df_rand$estimate[m1_df_rand$term=='sigma_O']
spatemp_sd = m1_df_rand$estimate[m1_df_rand$term=='sigma_E']
# still to add: r2s of training data, everything for hold out data, log likelihood values, variation in log likelihood values

df_m1 <- data.frame(label, len_t, len_h, rt, intercept, cof_pld, se_pld, cof_exp, se_exp, mrp, disp, spat_sd, spatemp_sd)

tbl_res <- '~/GIS/MSc_Projects/MPA_connectivity/scripts/particles_model/tbl_res.csv'
if (!file.exists(tbl_res)) {
  write.csv(df_m1, tbl_res, row.names=FALSE)
} else {
  dfr <- read.csv(tbl_res)
  dfr <- rbind(dfr, df_m1)
  write.csv(dfr, tbl_res, row.names=FALSE)
}








##########################################
# Try different model setups
##########################################


##########################################
# Model - test with tweedie dist
# This doesn't run, even with a small sample.

df_samp <- sample_n(df_t, 1000)
mesh <- make_mesh(df_samp, c('origin_x_km', 'origin_y_km'), cutoff = 5)
land_sf <- st_read(land)
bmesh <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=FALSE, range_fraction=0.2)
m2 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_samp,
  spde = bmesh, family = tweedie(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)


##########################################
# Model - compare m1 to if I just coded 0 values as 1 instead of as 1-352

df_samp <- df_t
df_samp$distance[df_samp$distance < 353] <- 1

mesh <- make_mesh(df_samp, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bmesh <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=FALSE, range_fraction=0.2)
m3 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_samp,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m3$gradients)>0.01){
  m3 <- sdmTMB::run_extra_optimization(m3, nlminb_loops = 1L, newton_loops = 1L)
}
saveRDS(m3, "m3.rds")
m3 <- readRDS("m3.rds")

summary(m3)
predictions3 <- predict(m3)

# Calculate r-squared
rsq <- function (x, y) cor(x, y) ^ 2
r2_m3 <- rsq(log(predictions3$distance), predictions3$est) #r2 better in m1

ggplot(predictions3, aes(x=distance, y=exp(est)))+
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  scale_y_log10()+scale_x_log10()

predictions3$residuals <- residuals(m3)
qqnorm(predictions3$residuals, ylim=c(-5,5));abline(a = 0, b = 1)
plot(density(predictions3$residuals)) # this is a good one for showing that we need to do something with those 1s
ggplot(predictions3, aes(x = est, y = residuals)) + geom_point()


##########################################
# Model - try a hurdle model
# (1) does it leave the source location (yes/no binomial)?
# (2) Then, if it does, how far does it go (gamma distribution, exclude 0 distance)?

# create settle_home field
df_t <- df_t %>% mutate(settle_home = case_when(distance<353 ~ 0, distance>352 ~ 1))

# binomial model
m4 <- sdmTMB(
  settle_home ~ 1 + log(total_exposure), data = df_t,  # it doesn't make sense to include pld. We can do this with hurdle models right (?) https://seananderson.ca/2014/05/18/gamma-hurdle/
  spde = bmesh, family = binomial(link='logit'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
saveRDS(m4, "m4.rds")
m4 <- readRDS("m4.rds")
summary(m4)
predictions4 <- predict(m4)
predictions4$est_res <- predictions4$est / (predictions4$est + 1.0) # predict is from Sean's package. You can only do things in link space. Convert back to response space.
predictions4$est_res2 <- plogis(predictions4$est)  #USE THIS ONE. The one above is NOT the same. See plot below:
ggplot(predictions4, aes(x=est_res, y=est_res2))+
  geom_point()


# notice how the jitter makes it clear that there is a bit of trend
ggplot(predictions4, aes(x=settle_home, y=est_res2))+
  geom_jitter(width=0.1, height=0) +
  geom_smooth(method = "glm",
              method.args = list(family = binomial(link = "logit")),
              se = FALSE, alpha = 0.7, colour = "black")

ggplot(predictions4, aes(x=log(total_exposure), y=est_res2))+
  geom_jitter(width=0.1, height=0) +
  geom_smooth(method = "glm",
              method.args = list(family = binomial(link = "logit")),
              se = FALSE, alpha = 0.7, colour = "black")

ggplot(predictions4, aes(log(total_exposure), est_res2)) +
  geom_line(alpha = 0.8) +
  geom_point(aes(y = settle_home), alpha = 0.3,
             position = position_jitter(height = 0.05)) +
  ylab("Probability of leaving home")



# Overall predictions.
predictions4 %>%
  ggplot()+
  geom_point(aes(x=origin_x, y=origin_y,color = plogis(est)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Overall")+
  facet_wrap(~month)

# Just the fixed effects. We DO see that higher exposure means you are more likely to leave home.
predictions4 %>%
  ggplot()+
  geom_point(aes(x=origin_x, y=origin_y,color = plogis(est_non_rf)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Fixed")+
  facet_wrap(~month)

# And for just the spatial component. We see that there are some areas that are more likely to leave than predicted by fixed effects.
predictions4 %>%
  ggplot()+
  geom_point(aes(x=origin_x, y=origin_y,color = plogis(omega_s)))+
  scale_color_viridis_c()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial component")+
  facet_wrap(~month)








## gamma model

# set distances less than 353 to zero in new field and filter them out
# I will need this field later.
df_t$distance_adj <- ifelse(df_t$distance <353, 0, df_t$distance)
df_samp <- filter(df_t, distance_adj > 0)

# need a new mesh since we have a different sample size
mesh5 <- make_mesh(df_samp, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bmesh5 <- add_barrier_mesh(spde_obj = mesh5, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=FALSE, range_fraction=0.2)
m5 <- sdmTMB(
  distance_adj ~ 1 + pld + log(total_exposure), data = df_samp,
  spde = bmesh5, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m5$gradients)>0.01){
  m5 <- sdmTMB::run_extra_optimization(m5, nlminb_loops = 1L, newton_loops = 1L)
} # took 4 times
saveRDS(m5, "m5.rds")
m5 <- readRDS("m5.rds")
summary(m5)
predictions5 <- predict(m5, newdata=df_t) # need to use df_t because otherwise we have different size dataframes when combining models.

# Multiply prob by distance predictions
df_t$pred_hurd <- predictions4$est_res2 * exp(predictions5$est)

# Look at just predicted vs. observed for just hte gamma model
df_t$gamma <- predictions5$est
ggplot(df_t %>% filter(distance_adj>0), aes(distance_adj, exp(gamma))) +
  geom_point()+
  scale_y_log10()+scale_x_log10()

# Plot against observed
ggplot(df_t, aes(distance_adj, pred_hurd)) +
  geom_point()+
  scale_y_log10()+scale_x_log10()
# we get a WARNING message here. I need to figure out where those negative values came from.
# Remove them for now
df_samp2 <- filter(df_t, pred_hurd>=0)
# Turn any zeros to ones so I can take the log
df_samp2$distance_adj <- ifelse(df_samp2$distance_adj == 0, 1, df_samp2$distance_adj)
# Get R2
rsq <- function (x, y) cor(x, y) ^ 2
r2_m5 <- rsq(log(df_samp2$distance_adj), log(df_samp2$pred_hurd))
# r2 is a lot lower compared to m1. In the plot you can see we get more extreme predicted values.
# AIC
AIC(m4) + AIC(m5)
# Compare to m1
# It has a lower AIC. No idea if I did this right though and they may not be comparable. Maybe only comparable to other hurdle models?





##########################################
# Try a model with only distances greater than 500m (~ the res of the SSC model and the particle distance analysis)

df_samp3 <- filter(df_t, distance>=500)

mesh8 <- make_mesh(df_samp3, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bmesh8 <- add_barrier_mesh(spde_obj = mesh8, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=FALSE, range_fraction=0.2)

m8 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_samp3,
  spde = bmesh8, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m8$gradients)>0.01){
  m8 <- sdmTMB::run_extra_optimization(m8, nlminb_loops = 1L, newton_loops = 1L)
}
saveRDS(m8, "m8.rds")
m8 <- readRDS("m8.rds")
summary(m8)
predictions8 <- predict(m8)

ggplot(predictions8, aes(x=distance, y=exp(est)))+
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  scale_y_log10()+scale_x_log10()

# Calculate r-squared
rsq <- function (x, y) cor(x, y) ^ 2
r2_m8 <- rsq(log(predictions8$distance), predictions8$est)

# residuals
predictions8$residuals <- residuals(m8)
qqnorm(predictions8$residuals, ylim=c(-5,5));abline(a = 0, b = 1)
plot(density(predictions8$residuals))
ggplot(predictions8, aes(x = est, y = residuals)) + geom_point()





