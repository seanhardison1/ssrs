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

rmse <- function(known, estimated){
  sqrt(mean((known - estimated)^2))
}

cb_poly <- st_read(here::here("data/cb_mesh_poly.kml")) %>% 
  sf::st_zm() %>% 
  st_transform(., crs = 4326) %>% 
  as_Spatial() 

cb_rast <- rasterize(cb_poly, modis) %>% 
  mask(modis, .) %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame() %>% 
  dplyr::rename(vals = 1,
                Longitude = 2,
                Latitude = 3)


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

cb_mesh_seg <- inla.sp2segment(cb_poly)
cb_mesh <- inla.mesh.2d(boundary=cb_mesh_seg, cutoff=0.01,
                        max.edge=c(0.075, 0.4))
# plot(cb_mesh, main = "")
cb_spde <- make_mesh(tss_train, c("Longitude","Latitude"), mesh = cb_mesh)

performance <- NULL
for (i in 1:10){
  tss_split <- initial_split(cb_wq_modis, prop = .8)
  tss_train <- training(tss_split)
  tss_test  <- testing(tss_split)
  known <- tss_test$daily_TSS
  
  # comparison with GAM
  message("GAM")
  mod1 <- gam(daily_TSS ~ s(vals) + 
                te(Longitude, Latitude, k = 20), 
              data = tss_train, method = "ML")
  gam_pred <- predict(mod1, newdata = tss_test, type = "response")

  message("spatial GLMM")
  grf_mod <- sdmTMB(daily_TSS ~ vals,
                     data = tss_train,
                     spde = cb_spde,
                     spatial_only = TRUE)
  smod_pred <- predict(grf_mod, tss_test) %>%
    pull(est)
  
  message("GLMs")
  reg_mod <- lm(daily_TSS ~ vals, data = tss_train)
  reg_pred <- predict(reg_mod, newdata = tss_test)
  
  glm_mod <- glm(daily_TSS ~ vals, family = Gamma("log"), data = tss_train)
  glm_pred <- predict(glm_mod, newdata = tss_test, type = "response")
  
  perf_i <- tibble(spat_mod_perf = rmse(smod_pred, known),
         gam_mod_perf = rmse(gam_pred, known),
         # gam_mod_perf2 = rmse(gam_pred2, known),
         lm_mod_perf = rmse(reg_pred, known),
         glm_mod_perf = rmse(glm_pred, known),
         k = i)
  
  assign("performance", rbind(performance, perf_i))
  print(perf_i)
}

save(performance, file = here::here("data/model_performance.rdata"))

smod_pred_t2 <- predict(grf_mod, tss_train) 

ggplot(smod_pred_t2) +
  geom_point(aes(y = daily_TSS, x = est))
