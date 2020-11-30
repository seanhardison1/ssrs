library(tidyverse)
library(sf)
library(raster)
library(mgcv)
library(anytime)
library(glmmTMB)
library(lubridate)
library(ggeffects)
source(here::here("R/sfc_as_cols.R"))

# read water quality data
cb_wq <- read.csv(here::here("data/WaterQualityWaterQualityCBSeg2003.csv"))

# load modis reflectance time series
load(here::here("data/modis_reflectance_ts.rdata"))

# pts <- cb_wq %>% 
#   group_by(Station, Latitude, Longitude) %>% 
#   dplyr::summarise(n = n()) %>% 
#   arrange(desc(n)) %>%
#   as.data.frame() %>% 
#   dplyr::slice(1:40) %>% 
#   st_as_sf(coords = c("Longitude","Latitude"), crs = 4326)
# write_sf(pts, here::here("data/station_locs2.shp"))

raw <- modis_reflec %>% 
  mutate(Longitude = Longitude*-1,
         datetime = as.Date(datetime, format = "%m/%d/%Y")) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(.,"+init=epsg:3857") %>% 
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL)

# sus_pix <- quantile(raw$vals, .977)

unique_refs <- modis_reflec %>% 
  dplyr::select(1,2) %>% 
  distinct() %>% 
  mutate(Longitude = Longitude*-1) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(.,"+init=epsg:3857")

cb_wq_modis <- cb_wq %>% 
  mutate(Longitude = round(Longitude, 4),
         Latitude = round(Latitude, 4),
         datetime = as.Date(SampleDate, format = "%m/%d/%Y")) %>% 
  filter(Parameter %in% c("TSS", "SALINITY"),
         Depth < 2) %>%
  # shrink data
  dplyr::select(Station, datetime, Longitude, Latitude, Parameter,
                MeasureValue, SampleDate, sal_prof = CBSeg2003, Depth, SampleReplicateType, Layer) %>% 
  
  # pivot variables to columns
  pivot_wider(., names_from = Parameter, values_from = MeasureValue, values_fn = mean) %>% 
  filter(!is.na(Longitude)) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(.,"+init=epsg:3857") %>% 
  st_snap(.,unique_refs, tolerance = 1000) %>% 
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL) %>% 
  left_join(.,raw, by = c("datetime","Longitude","Latitude")) %>% 
  filter(!is.na(vals)) %>% 
  group_by(sal_prof, Station, Longitude, Latitude, datetime, vals) %>% 
  dplyr::summarise(daily_TSS = mean(TSS, na.rm = T),
                   daily_salinity = mean(SALINITY, na.rm = T)) %>% 
  filter(!is.na(daily_TSS)) %>% 
         # vals < sus_pix) %>% 
  mutate(sal_regime = factor(ifelse(str_detect(sal_prof, "PH"),"PH","MH")),
         week = week(datetime))

save(cb_wq_modis, file = here::here("data/cb_wq_modis.rdata"))


