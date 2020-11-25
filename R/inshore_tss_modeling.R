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
load(here::here("data/inshore_cb_wq_modis.rdata"))

inshore_cb_wq_modis %<>% 
  filter(!is.na(daily_salinity), !is.na(vals), !is.na(daily_TSS)) %>% 
  st_as_sf(.,coords = c("Longitude","Latitude"), crs = 3857) %>% 
  st_transform(crs = 4326) %>% 
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL) %>% 
  mutate(year = year(datetime))

tss_split <- initial_split(inshore_cb_wq_modis, prop = .7)
tss_train <- training(tss_split)
tss_test  <- testing(tss_split)


# comparison with GAM
mod1 <- gam(daily_TSS ~ s(vals) + s(Longitude, k = 25),
            family = Gamma("log"), data = tss_train, method = "ML")

summary(mod1)
gam.check(mod1)
k.check(mod1)
sims <- simulateResiduals(mod1, n = 5000)
plot(sims)


mod2 <- gam(daily_TSS ~ s(vals),family = Gamma("log"),
            data = tss_train, method = "ML")
summary(mod2)
gam.check(mod2)
k.check(mod2)
sims <- simulateResiduals(mod2)
plot(sims)



test1 <- predict(mod1,newdata=tss_test, type = "response")
test2 <- predict(mod2,newdata=tss_test, type = "response")

(oos_pred <- 
  AIC(mod1, mod2) %>% 
  mutate(RMSE = c(sqrt(mean((tss_test$daily_TSS - test1)^2)),
                  sqrt(mean((tss_test$daily_TSS - test2)^2)))))

g1_rmse <- sqrt(mean((tss_test$daily_TSS - test1)^2))
g2_rmse <- sqrt(mean((tss_test$daily_TSS - test2)^2))

# models are nearly identical. defaulting to simpler models
inshore_gam_mod <- gam(daily_TSS ~ s(vals),
               family = Gamma("log"), data = tss_train, method = "ML")

# sdmTMB

# make mesh
cb_mesh_seg <- inla.sp2segment(cb_poly)
cb_mesh <- inla.mesh.2d(boundary=cb_mesh_seg, cutoff=0.01,
                        max.edge=c(0.1, 0.3))
cb_spde <- make_mesh(tss_train, c("Longitude","Latitude"), mesh = cb_mesh)

# coords <- as.matrix(distinct(inshore_cb_wq_modis[,c("Longitude","Latitude")]))
# plot(cb_mesh)
# points(coords, pch = 19, col = 2)

inshore_grf_mod <- sdmTMB(daily_TSS ~ s(vals),
                  family = Gamma("log"),
                  data = tss_train,
                  spde = cb_spde, 
                  spatial_only = TRUE)

# test model
smod_test <- predict(grf_mod, tss_test) %>% 
  mutate(est = exp(est))

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

p1 <- predict(grf_mod, cb_rast) %>% 
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

inshore_mod_rmse <- tibble(GLM = glm_rmse,
                   LR = reg_rmse,
                   GAM1 = g1_rmse,
                   GAM2 = g2_rmse,
                   `Spat-GLMM` = spat_rmse)

save(inshore_mod_rmse,
     inshore_gam_mod, 
     inshore_grf_mod,
     file = here::here("data/inshore_best_models.rdata"))
