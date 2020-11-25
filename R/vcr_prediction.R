library(tidyverse)
library(glmmTMB)
library(rsample)      # data splitting 
library(mgcv)
library(AICcmodavg)
library(sdmTMB)
library(sf)
library(raster)
library(INLA)
load(here::here("data/vcr-wq.rdata"))
load(here::here("data/vcr_wq_modis.rdata"))
load(here::here("data/best_models.rdata"))
load(here::here("data/inshore_best_models.rdata"))

## Spatially comprehensive models--------------

# Using the GAM to predict for the VCR
vcr_gam_pred <- predict(gam_mod, newdata = vcr_wq_modis, type = "response")

vcr_gam_pred_df <- vcr_wq_modis %>% 
  mutate(gam_pred = vcr_gam_pred)

sqrt(mean((vcr_wq_modis$daily_TSS - vcr_gam_pred)^2))

ggplot(vcr_gam_pred_df) +
  geom_point(aes(x = gam_pred, y = daily_TSS))

# Using the spatial GLMM to predict for the VCR
vcr_spat_pred_df <- predict(grf_mod, newdata = vcr_wq_modis) %>% 
  mutate(est = exp(est))

ggplot(vcr_spat_pred_df) +
  geom_point(aes(x = est, y = daily_TSS))

## Inshore models----------

# Using the inshore GAM to predict for the VCR
vcr_gam_pred_inshore <- predict(inshore_gam_mod, newdata = vcr_wq_modis, type = "response")

vcr_gam_pred_inshore_df <- vcr_wq_modis %>% 
  mutate(gam_pred = vcr_gam_pred_inshore )

sqrt(mean((vcr_wq_modis$daily_TSS - vcr_gam_pred_inshore)^2))

ggplot(vcr_gam_pred_inshore_df) +
  geom_point(aes(x = gam_pred, y = daily_TSS))

# Using the spatial GLMM to predict for the VCR
vcr_spat_pred_inshore_df <- predict(inshore_grf_mod, newdata = vcr_wq_modis) %>% 
  mutate(est = exp(est))

ggplot(vcr_spat_pred_inshore_df) +
  geom_point(aes(x = est, y = daily_TSS))

