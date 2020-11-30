library(tidyverse)
library(glmmTMB)
library(rsample)      # data splitting 
library(mgcv)
library(AICcmodavg)
library(sdmTMB)
library(sf)
library(lubridate)
library(raster)
library(INLA)
library(tsibble)
library(zoo)
library(DHARMa)

source(here::here("R/sfc_as_cols.R"))



cb_poly <- st_read(here::here("data/cb_mesh_poly.kml")) %>% 
  sf::st_zm() %>% 
  st_transform(., crs = 4326) %>% 
  as_Spatial() 

rescale <- function(x) {
  x <- x - mean(x, na.rm = T)
  x <- x / (1 * sd(x))
  # x <- x/sum(x)
  return(x)
}

ggplot(cb_wq_modis) +
  geom_density(aes(x = daily_TSS, fill = sal_regime))

modis <- stack(here::here("data/modis_rrs_645_stacked.tif"))[[8]]
# plot(stack(here::here("data/modis_rrs_645_stacked.tif")))


cb_rast <- rasterize(cb_poly, modis) %>% 
  mask(modis, .) %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame() %>% 
  dplyr::rename(vals = 1,
                Longitude = 2,
                Latitude = 3)

# load wq-modis data
load(here::here("data/cb_wq_modis.rdata"))

cb_wq_modis %<>% 
  filter(!is.na(vals), !is.na(daily_TSS)) %>% 
  st_as_sf(.,coords = c("Longitude","Latitude"), crs = 3857) %>% 
  st_transform(crs = 4326) %>% 
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL) %>% 
  mutate(year = year(datetime),
         yrmn = as.yearmon(as.Date(datetime)),
         yrmn_id = as.numeric(yrmn),
         yrwk = as.numeric(yearweek(datetime))) 
  

# tss_split <- initial_split(cb_wq_modis, prop = .8)
# tss_train <- training(tss_split)
# tss_test  <- testing(tss_split)
# print(tss_test %>% slice(1:10))

tss_train <- cb_wq_modis %>% dplyr::filter(!year %in% c(2016:2020))
tss_test <- cb_wq_modis %>% dplyr::filter(year %in% c(2016:2020))



# comparison with GAM
mod1 <- gam(daily_TSS ~ s(vals) + 
              te(Longitude, Latitude, k = 20), 
            data = tss_train, method = "ML")
plot(mod1)
summary(mod1)
gam.check(mod1)
k.check(mod1)
sims <- simulateResiduals(mod1)
plot(sims)

# best model 
gam_pred <- predict(mod1, newdata = tss_test, type = "response")
gam_mod <- mod1

# make mesh
cb_mesh_seg <- inla.sp2segment(cb_poly)
cb_mesh <- inla.mesh.2d(boundary=cb_mesh_seg, cutoff=0.01,
                         max.edge=c(0.05, 0.4))
plot(cb_mesh, main = "")
cb_spde <- make_mesh(tss_train, c("Longitude","Latitude"), mesh = cb_mesh)

coords <- as.matrix(distinct(cb_wq_modis[,c("Longitude","Latitude")]))
points(coords, pch = 19, col = 2)

grf_mod <- sdmTMB(daily_TSS ~ vals,
                  data = tss_train,
                  spde = cb_spde,
                  spatial_only = TRUE)
qqnorm(residuals(grf_mod), ylim = c(-3,8));abline(a = 0, b = 1)

# test model
smod_test <- predict(grf_mod, tss_test)
# 
grf_pred <- smod_test$est

# RMSE of spatial model 
(spat_rmse <- sqrt(mean((tss_test$daily_TSS- smod_test$est)^2)))

plot_map_raster <- function(dat, column = "est") {
  ggplot(dat, aes_string("Longitude", "Latitude", fill = column)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_viridis_c()
}




# baseline regression model
reg_mod <- lm(daily_TSS ~ vals, data = tss_train)
reg_pred <- predict(reg_mod, newdata = tss_test)
reg_rmse <- sqrt(mean((tss_test$daily_TSS - reg_pred)^2))

spatial_pred <- predict(grf_mod, cb_rast) %>% 
  mutate(gam_pred = predict(mod1, newdata = cb_rast, type = "response"),
         lm_pred = predict(reg_mod, newdata = cb_rast))
plot_map_raster(spatial_pred, "est")
plot_map_raster(spatial_pred, "gam_pred")

# baseline GLM
glm_mod <- glm(daily_TSS ~ vals,family = Gamma("log"), data = tss_train)
glm_pred <- predict(glm_mod, newdata = tss_test, type = "response")
glm_rmse <- sqrt(mean((tss_test$daily_TSS - glm_pred)^2))


save(
     spatial_pred,

     gam_mod,
     gam_pred,

     grf_mod,
     grf_pred,

     reg_mod,
     reg_pred,

     glm_mod,
     glm_pred,

     tss_test,
     file = here::here("data/best_models.rdata"))
