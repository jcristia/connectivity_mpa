
# test script to show Jacob/Patrick

# NOTE: following my meeting with Patrick on Dec 17, 2021, I made a few changes in this
# script as I worked through it with him. I added some meeting notes as comments.
# It is a bit of a mess, but I am leaving this version as is and then cleaning things
# up in the next version of the script.
# See notes in Evenote Chapter 2/3 note for some more details.

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




######
# TO TRY: (based on notes from Jacob meeting)
# see notes in Evernote for more info

# nothing spatial:
df_tg$mpa_part_id_orig <- factor(df_tg$mpa_part_id_orig)
m3 <- sdmTMB(
  distance ~ 1 + log(pld) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)

# with exposure
# THIS IS THE BEST MODEL. I WILL MOVE FORWARD WITH THIS ONE
m3 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m3

# with mpa_area
m4 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + log(mpa_area),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m4

# with mpa_area and id
m5 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + log(mpa_area) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m5

# with mpa_area and id and space
m6 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + log(mpa_area) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
m2<-m6

# mpa_area (km), no space
df_tg$mpa_area_adj <- (df_tg$mpa_area / 1000000)
m7 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + mpa_area_adj + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m7

# mpa_area (km) squared
m8 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + (mpa_area_adj)^2 + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m8

# mpa_area (km) log
m9 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + log(mpa_area_adj) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m9


#######
# then also try with ocean groupings
df_tg$group_id <- factor(df_tg$group_id)
m10 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + (1|group_id),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m10

m11 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + (1|group_id) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = FALSE,
  fields = "IID"
)
m2<-m11

m12 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + (1|group_id) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
m2<-m12


# run m3 but with space, this one makes the most sense
df_tg$mpa_part_id_orig <- factor(df_tg$mpa_part_id_orig)
m3 <- sdmTMB(
  distance ~ 1 + log(pld) + log(total_exposure) + (1|mpa_part_id_orig),
  data = df_tg,
  spde = bm,
  family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
m2<-m3


# run real m2 again to compare the r2
# 0.4512 - so the mpa id does help a bit

# so in summary, up to this point...
# the ocean grouping didn't help
# mpa_area helped a very small amount
# space actually makes things slightly worse
# BUT, nothing changed enough to move away from what is obvious and easy to explain
# which is m3
# If the group_id or area significantly improved things, then I would find a way to use them,
# but they don't, so it isn't worth trying to include them and have to come up with an explanation




########
# From Patrick: try with splines (gam), which is now implemented.
# Sean has a brief example of this.
m13 <- sdmTMB(
  log(distance) ~ 1 + pld + s(total_exposure),
  data = df_tg,
  spde = bm,
  #family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
m2<-m13
#0.4514

m14 <- sdmTMB(
  log(distance) ~ 1 + s(pld, k=8) + s(total_exposure), # didn't work unless I specified K, but then it just tells me that k is automaticlly determined
  data = df_tg,
  spde = bm,
  #family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
m2<-m14
#0.4644

m15 <- sdmTMB(
  log(distance) ~ 1 + s(pld, k=8) + s(total_exposure) + s(mpa_part_id_orig, bs='re'), #https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/random.effects.html
  data = df_tg,
  spde = bm,
  #family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
m2<-m15
#0.475
# seems to be good.
# qqplot looks worse on low end, but slightly better on high end.
# Est vs. Obs looks same


########################################
########################################

# Printing out model summary:
m3
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
predictions <- predict(m2, m2$data %>% filter(pld==3))   # filtered for just 1 PLD
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


##########################################
##########################################
# Map of predicted minus observed
predictions$diff <- ((predictions$est) - log(predictions$distance))
predictions <- predictions %>% arrange(desc(diff)) # rearange so that shorter distances draw on top
predictions %>%
  ggplot()+
  geom_sf(aes(color = diff))+
  scale_color_gradient2()+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Predicted minus observed")+
  facet_wrap(~month)
# I'll need to think about this. It seems like drawing order really makes a difference
# It might be better to make this map in Arc for the finished product.

# For pred vs. obs map, don't do absolute values. It would be good to know which
# side of zero they are on. However, Patrick doesn't think this map is useful
# anyways. I would want to do percent change, so (pred-obs)/obs. Then I may want
# to filter out any outliers and make sure there aren't any weird spatial
# patterns. He thinks it looks ok. I would want to go with just 1 pld to clean
# things ups. However, I think I should just drop this plot.


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
