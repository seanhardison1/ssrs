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
library(DHARMa)

source(here::here("R/sfc_as_cols.R"))



cb_poly <- st_read(here::here("data/inshore_cb_mesh_poly.kml")) %>% 
  sf::st_zm() %>% 
  st_transform(., crs = 4326) %>% 
  as_Spatial() 

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
load(here::here("data/ffs_cb_wq_modis.rdata"))

ffs_cb_wq_modis %<>% 
  filter(!is.na(daily_salinity), !is.na(vals), !is.na(daily_TSS)) %>% 
  st_as_sf(.,coords = c("Longitude","Latitude"), crs = 3857) %>% 
  st_transform(crs = 4326) %>% 
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL) %>% 
  mutate(year = year(datetime))

tss_split <- initial_split(ffs_cb_wq_modis, prop = .8)
tss_train <- training(tss_split)
tss_test  <- testing(tss_split)


# comparison with GAM
mod1 <- gam(daily_TSS ~ s(vals) + te(Longitude, Latitude, k = 20),
            family = Gamma("log"), data = tss_train, method = "ML")
summary(mod1)
gam.check(mod1)
k.check(mod1)
sims <- simulateResiduals(mod1)
plot(sims)


mod2 <- gam(daily_TSS ~ s(vals) + s(Longitude, k= 40),family = Gamma("log"),
            data = tss_train, method = "ML")
summary(mod2)
gam.check(mod2)
k.check(mod2)
sims <- simulateResiduals(mod2)
plot(sims)

mod3 <- gam(daily_TSS ~ te(vals, Longitude, k= 20),family = Gamma("log"),
            data = tss_train, method = "ML")
summary(mod3)
gam.check(mod3)
k.check(mod3)
sims <- simulateResiduals(mod3)
plot(sims)

test1 <- predict(mod1,newdata=tss_test, type = "response")
test2 <- predict(mod2,newdata=tss_test, type = "response")
test3 <- predict(mod3,newdata=tss_test, type = "response")

oos_pred <- 
  AIC(mod1, mod2, mod3) %>% 
  mutate(RMSE = c(sqrt(mean((tss_test$daily_TSS - test1)^2)),
                  sqrt(mean((tss_test$daily_TSS - test2)^2)),
                  sqrt(mean((tss_test$daily_TSS - test3)^2))))

g1_rmse <- sqrt(mean((tss_test$daily_TSS - test1)^2))
g2_rmse <- sqrt(mean((tss_test$daily_TSS - test2)^2))

# best model (lowest degrees of freedom, good out of sample accuracy)
gam_mod <- gam(daily_TSS ~ s(vals) + te(Longitude, Latitude, k = 20),
               family = Gamma("log"), data = tss_train, method = "ML")

gam_pred <- predict(gam_mod, newdata = tss_test, type = "response")
# make mesh
cb_mesh_seg <- inla.sp2segment(cb_poly)
cb_mesh <- inla.mesh.2d(boundary=cb_mesh_seg, cutoff=0.01,
                        max.edge=c(0.1, 0.3))
cb_spde <- make_mesh(tss_train, c("Longitude","Latitude"), mesh = cb_mesh)

# coords <- as.matrix(distinct(ffs_cb_wq_modis[,c("Longitude","Latitude")]))
# plot(cb_mesh)
# points(coords, pch = 19, col = 2)

ffs_grf_mod <- sdmTMB(daily_TSS ~ s(vals),
                  family = Gamma("log"),
                  data = tss_train,
                  spde = cb_spde, 
                  spatial_only = TRUE)

# test model
smod_test <- predict(ffs_grf_mod, tss_test) %>% 
  mutate(est = exp(est))

grf_pred <- smod_test$est

# RMSE of spatial model - better than GAMs
spat_rmse <- sqrt(mean((tss_test$daily_TSS - smod_test$est)^2))

# d <- cb_wq_modis
# d$residuals1 <- residuals(grf_mod)
# qqnorm(d$residuals1, ylim = c(-3,3));abline(a = 0, b = 1)

plot_map_raster <- function(dat, column = "est") {
  ggplot(dat, aes_string("Longitude", "Latitude", fill = column)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_viridis_c()
}

p1 <- predict(ffs_grf_mod, cb_rast) %>% 
  mutate(est = exp(est))
plot_map_raster(p1, "est")

# baseline regression model
reg_mod <- lm(daily_TSS ~ vals, data = tss_train)
reg_pred <- predict(reg_mod, newdata = tss_test)
reg_rmse <- sqrt(mean((tss_test$daily_TSS - reg_pred)^2))

# baseline GLM
glm_mod <- glm(daily_TSS ~ vals,family = Gamma("log"), data = tss_train)
glm_pred <- predict(glm_mod, newdata = tss_test, type = "response")
glm_rmse <- sqrt(mean((tss_test$daily_TSS - glm_pred)^2))

ffs_mod_rmse <- tibble(GLM = glm_rmse,
                   LR = reg_rmse,
                   GAM1 = g1_rmse,
                   GAM2 = g2_rmse,
                   `Spat-GLMM` = spat_rmse)

# save(ffs_mod_rmse,
#      ffs_gam_mod, 
#      ffs_grf_mod,
#      file = here::here("data/ffs_best_models.rdata"))
