# Patrick sent me this back after our first meeting.
# Most of the code is his.


library(tidyverse)
library(sdmTMB)
library(arrow)
library(patchwork)
library(sf)

file = './sample_1.feather'
df <- open_dataset(file, format='feather')
df <- read_feather(file)
df <- na.omit(df)

df <- filter(df, distance>0)
df <- filter(df, total_exposure>0) #are there really sites with 0 exposure?
df_sub <- df#sample_n(df, 10000, replace = FALSE)


df_sub$origin_x_km <- df_sub$origin_x/1000
df_sub$origin_y_km <- df_sub$origin_y/1000
mesh <- make_mesh(df_sub, c('origin_x_km', 'origin_y_km'), cutoff = 15)
# TODO: look into the mesh more: https://pbs-assess.github.io/sdmTMB/reference/make_mesh.html
plot(mesh)
#bmesh <- add_barrier_mesh(spde_obj = mesh, barrier_sf = coast_line, proj_scaling = 1000)


#first run model with month as a fixed effect and only a spatial field
m1 <- sdmTMB(
  distance ~ 0 + pld + log(total_exposure) + factor(month), data = df_sub,
  spde = mesh, family = Gamma(link='log'),
  #time = 'month',
  include_spatial = TRUE,
  spatial_only = TRUE, fields = "IID"
)

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


#second option - with spatiotemporal random field
m2 <- sdmTMB(
  distance ~ 1 + pld + log(total_exposure+1), data = df_sub,
  spde = mesh, family = Gamma(link='log'),
  time = 'month',
  include_spatial = TRUE,
  fields = "IID"
)

#this next part runs extra optimization if the model has not converged
if(max(m2$gradients)>0.01){
  m2 <- sdmTMB::run_extra_optimization(m2, nlminb_loops = 1L, newton_loops = 1L)
}

summary(m2)

predictions2 <- predict(m2)
predictions2 <- st_as_sf(predictions2, coords = c("origin_x", "origin_y"))
predictions2 <- st_set_crs(predictions2, value = "epsg:3005")

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

predictions2 %>%
  ggplot()+
  geom_sf(aes(color = exp(omega_s)))+
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  theme(legend.position = c(0,0), legend.justification = c(-0.05,-0.05))+
  ggtitle("Spatial component")+
  facet_wrap(~month)

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
