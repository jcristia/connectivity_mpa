
# Calculate proportion of variance explained by the fixed and random effects

# This script was written with sdmTMB version 0.1.4.9003, which differs from the
# version I did the original analysis with for the other script files. At this
# point I'm not going to go back and change code. I need this version for the
# updated sdmTMB residual and R2 scripts, which I don't think was fully
# developed when I worked on this in 2021/2022.

# The primary code is in a scratch folder in the sdmTMB git repo:
# https://github.com/pbs-assess/sdmTMB/blob/main/scratch/r2.R
# It's not finished, but it seems to work.

# This follows the approach documented here:
# https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/j.2041-210x.2012.00261.x
# https://rdrr.io/cran/MuMIn/man/r.squaredGLMM.html
# https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/
# https://www.biorxiv.org/content/10.1101/2022.04.19.488709v1.full.pdf (see how they report here on page 12)
# https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/


library(tidyverse)
library(sdmTMB)
library(arrow)
library(patchwork)
library(sf)
library(glue)

source('residuals_JC.R')


########################################
########################################
# Inputs

# feather files proportion label
label <- '03'
# coastline shapefile
land <- 'landmask_FINAL.shp'



########################################
########################################
# Data

file_t = glue('sample_{label}_training.feather')
file_h = glue('sample_{label}_holdout.feather')
df_t <- read_feather(file_t) # training data
df_h <- read_feather(file_h) # holdout data
df_t$origin_x_km <- df_t$origin_x/1000 # km required for mesh
df_t$origin_y_km <- df_t$origin_y/1000
df_h$origin_x_km <- df_h$origin_x/1000
df_h$origin_y_km <- df_h$origin_y/1000

# datasets without 0 distances for gamma model
df_tg <- filter(df_t, (distance > 2500 ))
df_hog <- filter(df_h, (distance > 2500 ))
df_tg$mpa_part_id_orig <- factor(df_tg$mpa_part_id_orig)
df_hog$mpa_part_id_orig <- factor(df_hog$mpa_part_id_orig)


########################################
########################################
# Model used in paper

# Create mesh and barrier mesh
land_sf <- st_read(land, quiet=TRUE)
mesh <- make_mesh(df_tg, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bm <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                       proj_scaling = 1000, plot=FALSE, range_fraction=0.2)

# # sdmTMB model
# print('Running distance model')
# start_time <- Sys.time()
# m2 <- sdmTMB(
#   distance ~ 1 + log(pld) + log(total_exposure),
#   data = df_tg,
#   mesh = bm,
#   family = Gamma(link='log'),
#   time = 'month',
#   spatial = 'on',
#   spatiotemporal = "iid"
# )
# i <- 1
# while(i<30){ # max out at 30 optimizations
#   if(max(m2$gradients)>0.01){
#     print("Running extra optimization")
#     m2 <- sdmTMB::run_extra_optimization(m2, nlminb_loops = 1L, newton_loops = 1L)
#     print(max(m2$gradients))
#   }
#   i <- i+1
# }
# end_time <- Sys.time()
# runtime_g <- (end_time - start_time)[[1]]
# print(paste('Distance model runtime:', runtime_g))
# saveRDS(m2, glue("m_{label}_gamma.rds"))
# m2 <- readRDS(glue("m_{label}_gamma.rds"))
#
# # model summary
# m2
# tidy(m2, conf.int = TRUE)
# m2$sd_report
# tidy(m2, 'ran_pars', conf.int=TRUE)
# tidy(m2, 'ran_vals', conf.int=TRUE) # random intercepts by mpa ID
#
#
# # R2 - predicted vs. observed
# # Predict. In this case we want to calculate R2 with our hold out data. We are
# # most interested in how well our model can fit new data.
# predictions <- predict(m2, newdata=df_hog)
# #predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
# #predictions <- st_set_crs(predictions, value = "epsg:3005")
#
# # R-squared
# rsq <- function (x, y) cor(x, y) ^ 2
# r2_g <- rsq(log(predictions$distance), predictions$est)
# # R2 with holdout data
# predictions_ho <- predict(m2, newdata=df_hog)
# predictions_ho <- st_as_sf(predictions_ho, coords = c("origin_x", "origin_y"))
# predictions_ho <- st_set_crs(predictions_ho, value = "epsg:3005")
# r2_ho <- rsq(log(predictions_ho$distance), predictions_ho$est)




########################################
########################################
# Residual method type

# Before testing different models, first decide which residuals method to use.
# See in the notes here:
# https://github.com/pbs-assess/sdmTMB/blob/main/R/residuals.R
# Either use 'mle-leplace' or 'mle-mcmc'
# Mle-mcmc is theoretically preferred, but it is slower.
# Run with 100 iterations and see if it is close to mle-leplace.
# It takes a LONG time.

# NOTE!:
# For the sdmTMB residuals function, there doesn't seem to be a way to pass it
# the newdata predictions. In the residuals script there is even a spot where
# it gets set to null.
# This means you can only calculate marginal r2 for the training data.

# I've now changed the residuals script to allow for new data.


fit <- m2
fixef <- function(x) {
  b <- tidy(fit)
  stats::setNames(b$estimate, b$term)
}

fe <- fixef(fit)
X <- fit$tmb_data$X_ij
VarF <- var(as.vector(fixef(fit) %*% t(X[[1]]))) # variance from fixed-effects   # IMPORTANT TO INDEX [[1]]
b <- tidy(fit, "ran_par")
sigma_O <- b$estimate[b$term == "sigma_O"] # spatial standard deviation
sigma_E <- b$estimate[b$term == "sigma_E"] # spatiotemporal standard deviation
#sigma_G <- b$estimate[b$term == 'sigma_G'] # Random intercepts
VarO <- sigma_O^2
VarE <- sigma_E^2
#VarG <- sigma_G^2

# MLE-LAPLACE
res <- residuals(fit, type='mle-laplace')
length(res[is.infinite(res)]) # I get a few infinite values I need to remove. Not sure why this is happening but it is only a few out of hundreds of thousands.
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf)) # residual variance
r2_fix <- VarF/(VarF + VarO + VarE + VarR) # fixed effects
#r2_rin <- VarG/(VarF + VarG + VarO + VarE + VarR) # random intercepts
r2_spa <- VarO/(VarF + VarO + VarE + VarR) # spatial random effect
r2_spt <- VarE/(VarF + VarO + VarE + VarR) # spatiotemporal
df_lap <- data.frame(res='mle_laplace', varR=VarR, r2_fix=r2_fix, r2_spa=r2_spa, r2_spt=r2_spt)

# MLE-MCMC
res <- residuals(fit, type='mle-mcmc', mcmc_iter = 101, mcmc_warmup = 100)
saveRDS(res, glue("m_{label}_residuals_mlemcmc.rds"))
res <- readRDS(glue("m_{label}_residuals_mlemcmc.rds"))
length(res[is.infinite(res)])
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf)) # residual variance
r2_fix <- VarF/(VarF + VarO + VarE + VarR) # fixed effects
#r2_rin <- VarG/(VarF + VarG + VarO + VarE + VarR) # random intercepts
r2_spa <- VarO/(VarF + VarO + VarE + VarR) # spatial random effect
r2_spt <- VarE/(VarF + VarO + VarE + VarR) # spatiotemporal
df_mcmc <- data.frame(res='mle_mcmc', varR=VarR, r2_fix=r2_fix, r2_spa=r2_spa, r2_spt=r2_spt)

# This is apparently just observed minus predicted
res <- residuals(fit, type='response')
length(res[is.infinite(res)])
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf)) # residual variance
r2_fix <- VarF/(VarF + VarO + VarE + VarR) # fixed effects
#r2_rin <- VarG/(VarF + VarG + VarO + VarE + VarR) # random intercepts
r2_spa <- VarO/(VarF + VarO + VarE + VarR) # spatial random effect
r2_spt <- VarE/(VarF + VarO + VarE + VarR) # spatiotemporal
df_response <- data.frame(res='mle_response', varR=VarR, r2_fix=r2_fix, r2_spa=r2_spa, r2_spt=r2_spt)

# Interpretation:
# The laplace and mcmc methods nearly match, so I will just use the laplace method.







#############################################
#############################################
# Proportion of variance for different model set ups

# Model variations:
# Full model
# No spatiotemporal
# No spatial
# No MPA ID random effect
# Just PLD
# Just Exposure

# For each variation, calculate:
# PredObs R2
# Fixed effects R2 (i.e. marginal R2)
# Random intercept R2
# Spat R2
# Spatiotemporal R2
# (Sean did not have sigma_G random intercepts in his code, so I'm not 100% sure
# this is correct, however, the numbers appear to make sense)

# UPDATE 2023.03.27 - no longer have a random intercept in the model (mpa origin id)

# After the fact, I did all of this with and without the training data
# VarF and VarR data is compiled differently depending on which dataset is used.
training_data <- FALSE


fixef <- function(x) {
  b <- tidy(fit)
  stats::setNames(b$estimate, b$term)
}
rsq <- function (x, y) cor(x, y) ^ 2


########## FUll model
m_full <- readRDS(glue("m_{label}_gamma.rds"))
fit <- m_full

if(training_data){
  predictions <- predict(fit)
  } else{
predictions <- predict(fit, newdata=df_hog)}

r2_predobs <- rsq(log(predictions$distance), predictions$est)

# Test what var residuals look like when doing just observed minus predicted.
# Just response residuals
test_varR <- var(log(predictions$distance) - predictions$est)
# This differs a bit from the mle-laplace method. Guessing it is something with
# the structure of residuals across random "levels"(?)


fe <- fixef(fit)

if(training_data){
  X <- fit$tmb_data$X_ij
  VarF <- var(as.vector(fixef(fit) %*% t(X[[1]])))
} else {
  X <- predictions %>%
    mutate(pld_log = log(pld),
           exposure_log = log(total_exposure),
           intercept = 1) %>%
    select(intercept, pld_log, exposure_log)
  VarF <- var(as.vector(fixef(fit) %*% t(X)))}

b <- tidy(fit, "ran_par")
sigma_O <- b$estimate[b$term == "sigma_O"]
sigma_E <- b$estimate[b$term == "sigma_E"]
#sigma_G <- b$estimate[b$term == 'sigma_G']
VarO <- sigma_O^2
VarE <- sigma_E^2
#VarG <- sigma_G^2

if(training_data){
  res <- residuals_JC(fit, type='mle-laplace')
} else {
  res <- residuals_JC(fit, nd = df_hog, type='mle-laplace')
}
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf))

r2_fix <- VarF/(VarF + VarO + VarE + VarR)
#r2_rin <- VarG/(VarF + VarG + VarO + VarE + VarR)
r2_spa <- VarO/(VarF + VarO + VarE + VarR)
r2_spt <- VarE/(VarF + VarO + VarE + VarR)
df_mfull <- data.frame(model='full', r2_spt=r2_spt, r2_spa=r2_spa, r2_fix=r2_fix, r2_predobs=r2_predobs)



########## No spatiotemporal
#df_tg$mpa_part_id_orig <- factor(df_tg$mpa_part_id_orig)
m_nospt <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure),
  data = df_tg,
  mesh = bm,
  family = Gamma(link='log'),
  time = NULL,
  spatial = 'on',
  spatiotemporal = "off"
)
i <- 1
while(i<30){
  if(max(m_nospt$gradients)>0.01){
    print("Running extra optimization")
    m_nospt <- sdmTMB::run_extra_optimization(m_nospt, nlminb_loops = 1L, newton_loops = 1L)
    print(max(m_nospt$gradients))
  }
  i <- i+1
}

fit <- m_nospt
if(training_data){
  predictions <- predict(fit)
} else{
  predictions <- predict(fit, newdata=df_hog)}
r2_predobs <- rsq(log(predictions$distance), predictions$est)

if(training_data){
  X <- fit$tmb_data$X_ij
  VarF <- var(as.vector(fixef(fit) %*% t(X[[1]])))
} else {
  X <- predictions %>%
    mutate(pld_log = log(pld),
           exposure_log = log(total_exposure),
           intercept = 1) %>%
    select(intercept, pld_log, exposure_log)
  VarF <- var(as.vector(fixef(fit) %*% t(X)))}

b <- tidy(fit, "ran_par")
sigma_O <- b$estimate[b$term == "sigma_O"]
VarO <- sigma_O^2

if(training_data){
  res <- residuals_JC(fit, type='mle-laplace')
} else {
  res <- residuals_JC(fit, nd = df_hog, type='mle-laplace')
}
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf))
r2_fix <- VarF/(VarF + VarO + VarR)
r2_spa <- VarO/(VarF + VarO + VarR)
df_mnospt <- data.frame(model='no_spatiotemp', r2_predobs=r2_predobs, r2_fix=r2_fix, r2_spa=r2_spa)



# ########## No MPA ID random intercept
# # I'm removing this before removing the spatial random effect because I know
# # that this random intercept is largely spatial (MPA location), so I want to
# # know how much of the variation gets sucked up by the spatial random field.
# m_noran <- sdmTMB(
#   distance ~ 1 + log(pld) + log(total_exposure),
#   data = df_tg,
#   mesh = bm,
#   family = Gamma(link='log'),
#   time = NULL,
#   spatial = 'on',
#   spatiotemporal = "off"
# )
# i <- 1
# while(i<30){
#   if(max(m_noran$gradients)>0.01){
#     print("Running extra optimization")
#     m_noran <- sdmTMB::run_extra_optimization(m_noran, nlminb_loops = 1L, newton_loops = 1L)
#     print(max(m_noran$gradients))
#   }
#   i <- i+1
# }
#
# fit <- m_noran
# predictions <- predict(fit, newdata=df_hog)
# #predictions <- predict(fit)
# r2_predobs <- rsq(log(predictions$distance), predictions$est)
#
# X <- fit$tmb_data$X_ij
# VarF <- var(as.vector(fixef(fit) %*% t(X[[1]])))
# b <- tidy(fit, "ran_par")
# sigma_O <- b$estimate[b$term == "sigma_O"]
# VarO <- sigma_O^2
#
# res <- residuals(fit, type='mle-laplace')
# res_noinf <- res[!is.infinite(res)]
# VarR <- var(as.vector(res_noinf))
# r2_fix <- VarF/(VarF + VarO + VarR)
# r2_spa <- VarO/(VarF + VarO + VarR)
# df_mnorandint <- data.frame(model='no_randint', r2_predobs=r2_predobs, r2_fix=r2_fix, r2_spa=r2_spa)



########## No spatial
m_nospa <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure),
  data = df_tg,
  mesh = bm,
  family = Gamma(link='log'),
  time = NULL,
  spatial = 'off',
  spatiotemporal = "off"
)
i <- 1
while(i<30){
  if(max(m_nospa$gradients)>0.01){
    print("Running extra optimization")
    m_nospa <- sdmTMB::run_extra_optimization(m_nospa, nlminb_loops = 1L, newton_loops = 1L)
    print(max(m_nospa$gradients))
  }
  i <- i+1
}

fit <- m_nospa
if(training_data){
  predictions <- predict(fit)
} else{
  predictions <- predict(fit, newdata=df_hog)}
r2_predobs <- rsq(log(predictions$distance), predictions$est)

if(training_data){
  X <- fit$tmb_data$X_ij
  VarF <- var(as.vector(fixef(fit) %*% t(X[[1]])))
} else {
  X <- predictions %>%
    mutate(pld_log = log(pld),
           exposure_log = log(total_exposure),
           intercept = 1) %>%
    select(intercept, pld_log, exposure_log)
  VarF <- var(as.vector(fixef(fit) %*% t(X)))}

if(training_data){
  res <- residuals_JC(fit, type='mle-laplace')
} else {
  res <- residuals_JC(fit, nd = df_hog, type='mle-laplace')
}
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf))
r2_fix <- VarF/(VarF + VarR)
df_mnospa <- data.frame(model='no_spatial', r2_predobs=r2_predobs, r2_fix=r2_fix)



########## Just exposure
m_expo <- sdmTMB(
  distance ~ 1 + log(total_exposure),
  data = df_tg,
  mesh = bm,
  family = Gamma(link='log'),
  time = NULL,
  spatial = 'off',
  spatiotemporal = "off"
)
i <- 1
while(i<30){
  if(max(m_expo$gradients)>0.01){
    print("Running extra optimization")
    m_expo <- sdmTMB::run_extra_optimization(m_expo, nlminb_loops = 1L, newton_loops = 1L)
    print(max(m_expo$gradients))
  }
  i <- i+1
}

fit <- m_expo
if(training_data){
  predictions <- predict(fit)
} else{
  predictions <- predict(fit, newdata=df_hog)}
r2_predobs <- rsq(log(predictions$distance), predictions$est)

if(training_data){
  X <- fit$tmb_data$X_ij
  VarF <- var(as.vector(fixef(fit) %*% t(X[[1]])))
} else {
  X <- predictions %>%
    mutate(exposure_log = log(total_exposure),
           intercept = 1) %>%
    select(intercept, exposure_log)
  VarF <- var(as.vector(fixef(fit) %*% t(X)))}

if(training_data){
  res <- residuals_JC(fit, type='mle-laplace')
} else {
  res <- residuals_JC(fit, nd = df_hog, type='mle-laplace')
}
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf))
r2_fix <- VarF/(VarF + VarR)
df_expo <- data.frame(model='fix_exposure', r2_predobs=r2_predobs, r2_fix=r2_fix)



########## Just pld
m_pld <- sdmTMB(
  distance ~ 1 + log(pld),
  data = df_tg,
  mesh = bm,
  family = Gamma(link='log'),
  time = NULL,
  spatial = 'off',
  spatiotemporal = "off"
)
i <- 1
while(i<30){
  if(max(m_pld$gradients)>0.01){
    print("Running extra optimization")
    m_pld <- sdmTMB::run_extra_optimization(m_pld, nlminb_loops = 1L, newton_loops = 1L)
    print(max(m_pld$gradients))
  }
  i <- i+1
}

fit <- m_pld
if(training_data){
  predictions <- predict(fit)
} else{
  predictions <- predict(fit, newdata=df_hog)}
r2_predobs <- rsq(log(predictions$distance), predictions$est)

if(training_data){
  X <- fit$tmb_data$X_ij
  VarF <- var(as.vector(fixef(fit) %*% t(X[[1]])))
} else {
  X <- predictions %>%
    mutate(pld_log = log(pld),
           intercept = 1) %>%
    select(intercept, pld_log)
  VarF <- var(as.vector(fixef(fit) %*% t(X)))}

if(training_data){
  res <- residuals_JC(fit, type='mle-laplace')
} else {
  res <- residuals_JC(fit, nd = df_hog, type='mle-laplace')
}
res_noinf <- res[!is.infinite(res)]
VarR <- var(as.vector(res_noinf))
r2_fix <- VarF/(VarF + VarR)
df_pld <- data.frame(model='fix_pld', r2_predobs=r2_predobs, r2_fix=r2_fix)


################


df_all <- bind_rows(df_mfull, df_mnospt, df_mnospa, df_expo, df_pld)

# Then for PredObs, calculate the % change between each model and the full model
# to somewhat compare the predobs R2 with the marginal R2.

df_all <- df_all %>% mutate(r2_predobs_prop = (r2_predobs - r2_predobs[[1]]) / r2_predobs[[1]])

if(training_data){
write.csv(df_all, 'tbl_r2_trainingData.csv', row.names=FALSE)
}else{
  write.csv(df_all, 'tbl_r2_holdoutData.csv', row.names=FALSE)
}

###########################################

# Interpreting and general notes

# Patrick thinks an R2 above 0.2 is still good.
# Random fields explain 97+% of the variation.
# He thinks fixed effects explaining so little makes sense given the complexity of the coast.
# He also is impressed with how much the spatial random field explains.

# You'll see that once you remove spatial random fields, the variance explained
# by the fixed effects actually improves (0.08 to 4.4% with the training data).
# Apparenlty the spatial random fields tend to suck up more variation than they
# should (he has a reference for this is one of his papers, Thompson 2022).

# "There is, however, the risk that the spatial random fields can absorb
# variation in CPUE that is associated with climate variables, which would result
# in underestimating the projected impacts of climate change (Clayton et al., 1993)."

# AND THIS IS WHY we take this approach (i.e. looking at how the pseudo R2
# changes as we remove predictors, as opposed to just doing it for the full model)

# Regarding how I frame the overall results of the model: Patrick is actually
# quite positive about it. He likes Figure 8. He says that with the random
# effects included, it's actually a good model, because in a way the circulation
# model is included in those random effects. So, for this study area
# specifically, we can make decent predictions if we know the location. THIS is
# why he says that in stock assessment, Sean often just fits his models with
# space and time random effects. They don't care what might be causing
# spatial/temporal patterns, they just want to know if there are patterns. Then,
# they can say if we sampled at this location at this time, here is what we
# would expect to get.

# My model wouldn't necessarily work for predicting connectivity across the
# globe, but it DOES work for predicting connectivity on the coast of BC.



