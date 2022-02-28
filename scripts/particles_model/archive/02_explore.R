# Building on the last script...
# Test different model setups and mesh resolutions.


library(tidyverse)
library(sdmTMB)
library(arrow)
library(patchwork)
library(sf)


########################################
# Read, clean up, sample

file = './sample_1.feather'
#df <- open_dataset(file, format='feather') # will probably need this method when using the full dataset
df <- read_feather(file)

df <- na.omit(df) # both distance and settle have na. I'm just working with the distance model now, so I can eliminate all of them.
df <- filter(df, total_exposure>0) #are there really sites with 0 exposure?
# there are ~4000 points out of 2.7 million that have 0 exposure. These are indeed
# incorrect, but it is not worth fixing. Most points are on the land border but
# others are not. There are also points on the border that have valid exposure
# values, so it is not clear what the issue is. There is likely some geometry
# precision problem in the fetch script.
df$distance[df$distance == 0] <- 1 # YES, we can have distances of zero. Ref points
# are sometimes the same, but for these, I'm going to make them all 1 though. Cell
# size was 500m, so the particles could have moved up to that distance. Make these
# one and then we can work with the gamma distribution.

# take a smaller sample?
#df_sub <- sample_n(df, 1000, replace = FALSE)
df_sub <- df



#########################################
# Create mesh and barrier mesh

df_sub$origin_x_km <- df_sub$origin_x/1000
df_sub$origin_y_km <- df_sub$origin_y/1000
mesh <- make_mesh(df_sub, c('origin_x_km', 'origin_y_km'), cutoff = 10)
plot(mesh)
# add barrier:
land <- 'C:/Users/jcristia/Documents/GIS/MSc_Projects/MPA_connectivity/spatial/Coastline/landmask_FINAL.shp'
land_sf <- st_read(land)
bmesh <- add_barrier_mesh(spde_obj = mesh, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=TRUE, range_fraction=0.2)
# range fraction? In his example it looks like we estimate the proportion of a
# triangle that will be land compared to water? Might need to play around with
# this. The significance will probably also change once we lower our cutoff.

# plot out triangle points
mesh_df_water <- bmesh$mesh_sf[bmesh$normal_triangles, ]
mesh_df_land <- bmesh$mesh_sf[bmesh$barrier_triangles, ]
ggplot(land_sf) +
  geom_sf() +
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  geom_sf(data = mesh_df_land, size = 1, colour = "green")





##########################################
# Models
##########################################




##########################################
#first run model with month as a fixed effect and only a spatial field

m1 <- sdmTMB(
  distance ~ 0 + pld + log(total_exposure) + factor(month), data = df_sub,
  spde = bmesh, family = Gamma(link='log'),
  #time = 'month',
  include_spatial = TRUE,
  spatial_only = TRUE, fields = "IID"
)
#So setting the intercept to zero in the model where we include month as a fixed effect makes sense when month is a factor data type. As a factor, it estimates a coefficient for each factor level (essentially making it a "grouping"??). Intercept of 0 and 1 are the same result I think, it is just more interpretable when you set it to 0. If it is 1, then the "Intercept" in the summary becomes the first month, and then the others are shown in relation to that. If you set it to 0 then we get them in relation to an intercept value of 0.

#this next part runs extra optimization if the model has not converged
if(max(m1$gradients)>0.01){
  m1 <- sdmTMB::run_extra_optimization(m1, nlminb_loops = 1L, newton_loops = 1L)
}

summary(m1)

predictions <- predict(m1)
predictions <- st_as_sf(predictions, coords = c("origin_x", "origin_y"))
predictions <- st_set_crs(predictions, value = "epsg:3005")

p1 <- predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(est)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Overall predicted")

p2 <- predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(est_non_rf)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Fixed effects component")

p3 <- predictions %>%
  ggplot()+
  geom_sf(aes(color = exp(omega_s)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial component")

p1+p2+p3



##########################################
#second option - with spatiotemporal random field

m2 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_sub,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This is how you do cross_validation
# you get a log likelihood (need to look up more about this)
# It will split data in 5 and run
m2_cv <- sdmTMB_cv(
  distance ~ 1 + pld + log(total_exposure), data = df_sub,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID",k_folds = 5
)
# will give me predictive power
# log lk - diff between pred and obv and take log (asking: what's the likelihood of getting that difference?)
# you want the model that has the highest summed loglik
# Sean sums his, but Patrick thinks averaging makes more sense.

# However, Patrick thinks it might be better to report an r2
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if(max(m2$gradients)>0.01){
  m2 <- sdmTMB::run_extra_optimization(m2, nlminb_loops = 1L, newton_loops = 1L)
}

summary(m2)

predictions2 <- predict(m2)
predictions2 <- st_as_sf(predictions2, coords = c("origin_x", "origin_y"))
predictions2 <- st_set_crs(predictions2, value = "epsg:3005")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# plot predicted vs. observed
ggplot(predictions2 %>% filter(distance>1), aes(x=distance, y=exp(est)))+
  geom_point() +
  scale_y_log10()+scale_x_log10()
# can then calc R2 on this

# We see that it has a hard time make predictions for ones that went 0 distance.
# So.... we could either:

# Do a tweedie model
# or
# do a hurdle model (ask, does it settle in source location? Then, if it doesn't, how far does it go?)
# So I would be excluding those with distances of 1 in the second model.
# can multiply the 2 together to get a prediction

# And actually, maybe try it by setting distances of 0 between 1 and 250
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

predictions2 %>%
  ggplot()+
  geom_sf(aes(color = exp(est)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Overall predicted")+
  facet_wrap(~month)

predictions2 %>%
  ggplot()+
  geom_sf(aes(color = exp(est_non_rf)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Fixed effects component")+
  facet_wrap(~month)

# https://pbs-assess.github.io/sdmTMB/articles/index-standardization.html
# omega_s
# Spatial random effects only
# "We can look at the spatial random effects that represent consistent
# deviations in space through time that are not accounted for by our fixed
# effects. In other words, these deviations represent consistent biotic and
# abiotic factors that are affecting biomass density but are not accounted for
# in the model."
predictions2 %>%
  ggplot()+
  geom_sf(aes(color = exp(omega_s)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial component")+
  facet_wrap(~month)

# epsilon_st
# Spatiotemporal random effects only
# "spatiotemporal random effects that represent deviation from the fixed effect
# predictions and the spatial random effect deviations. These represent biotic
# and abiotic factors that are changing through time and are not accounted for
# in the model.
predictions2 %>%
  ggplot()+
  geom_sf(aes(color = exp(epsilon_st)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatiotemporal component")+
  facet_wrap(~month)

#compare models####
AIC(m1, m2) #m2 has a lower AIC so is probably a better choice as both models are consistent with your hypothesis




##########################################
# add year as fixed effect

m3 <- sdmTMB(
  distance ~ 0 + pld + log(total_exposure) + factor(year), data = df_sub,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
if(max(m3$gradients)>0.01){
  m3 <- sdmTMB::run_extra_optimization(m3, nlminb_loops = 1L, newton_loops = 1L)
}
summary(m3)
AIC(m2, m3) # m2 still looks better



##########################################
# m2 with AR1 family

#Does AR1 mean that we think previous time steps can be predictive of future time steps? I don't think that is true in our case and we can probably stick with IID.
m4 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_sub,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "AR1"
)
if(max(m4$gradients)>0.01){
  m4 <- sdmTMB::run_extra_optimization(m4, nlminb_loops = 1L, newton_loops = 1L)
} # had to run this 5 times before it would converge
summary(m4) # The rho is 0.99, which the documentation says you should just switch to Random Walk in that case. (but doesn't a high value mean that we see autocorrelation?)
AIC(m2, m4)
# it seems to be better, so I guess I need to decide if it makes sense to do AR1



##########################################
# m2 with RW family
m5 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_sub,
  spde = bmesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "RW"
)
if(max(m5$gradients)>0.01){
  m5 <- sdmTMB::run_extra_optimization(m5, nlminb_loops = 1L, newton_loops = 1L)
} # had to run this 3 times before it would converge
summary(m5)
AIC(m2, m4, m5) # this one seems even better




##########################################
# see how a higher resolution bmesh changes things
##########################################

# my origin reference points were spaced 500 meters apart
# My mesh above is 10km. Compare models with 5, 1, 0.5 km.


mesh5 <- make_mesh(df_sub, c('origin_x_km', 'origin_y_km'), cutoff = 5)
bmesh5 <- add_barrier_mesh(spde_obj = mesh5, barrier_sf = land_sf,
                          proj_scaling = 1000, plot=TRUE, range_fraction=0.2)
# mesh1 <- make_mesh(df_sub, c('origin_x_km', 'origin_y_km'), cutoff = 1)
# bmesh1 <- add_barrier_mesh(spde_obj = mesh1, barrier_sf = land_sf,
#                            proj_scaling = 1000, plot=TRUE, range_fraction=0.2)

m6 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_sub,
  spde = bmesh5, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
# start_time <- Sys.time()
# m7 <- sdmTMB(
#   distance ~ 1 + pld + log(total_exposure), data = df_sub,
#   spde = bmesh1, family = Gamma(link='log'),
#   time = 'month',
#   include_spatial = TRUE,
#   fields = "IID"
# )
# if(max(m7$gradients)>0.01){
#   m7 <- sdmTMB::run_extra_optimization(m7, nlminb_loops = 1L, newton_loops = 1L)
# }
# end_time <- Sys.time()
# THIS crashes every time

summary(m6)
AIC(m2, m6)


# Try with 2.5
mesh25 <- make_mesh(df_sub, c('origin_x_km', 'origin_y_km'), cutoff = 2.5)
bmesh25 <- add_barrier_mesh(spde_obj = mesh25, barrier_sf = land_sf,
                           proj_scaling = 1000, plot=TRUE, range_fraction=0.2)
start_time <- Sys.time()
m7 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure), data = df_sub,
  spde = bmesh25, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)
end_time <- Sys.time()
end_time - start_time
if(max(m7$gradients)>0.01){
  m7 <- sdmTMB::run_extra_optimization(m7, nlminb_loops = 1L, newton_loops = 1L)
}
# take 7 minutes
# I'm thinking that 5km is a good cutoff to move forward with for more testing,
# but this might change after I talk to Patrick.




##########################################
# now test with progressively larger datasets
##########################################

# use the m6 model with a 5km cutoff














