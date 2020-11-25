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
load(here::here("data/vcr_reflectance_ts.rdata"))
source(here::here("R/sfc_as_cols.R"))

raw <- vcr_reflec %>% 
  mutate(Longitude = Longitude*-1,
         datetime = as.Date(datetime, format = "%m/%d/%Y")) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(.,"+init=epsg:3857") %>% 
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL)

unique_refs <- vcr_reflec %>% 
  dplyr::select(1,2) %>% 
  distinct() %>% 
  mutate(Longitude = Longitude*-1) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(.,"+init=epsg:3857")

vcr_wq_modis <- wq_summarized %>% 
  st_set_crs(4326) %>% 
  sfc_as_cols(names = c("Longitude","Latitude")) %>% 
  mutate(Longitude = round(Longitude, 4),
         Latitude = round(Latitude, 4),
         datetime = as.Date(measureDate, format = "%m/%d/%Y")) %>% 
  dplyr::select(tss, station, datetime, Longitude, Latitude, depth) %>% 

  filter(!is.na(Longitude),
         station %in% c("SS", "QI", "SSI", "NM", "LC", "MI")) %>% 
  st_transform(.,"+init=epsg:3857") %>% 
  st_snap(.,unique_refs, tolerance = 1000) %>% 
  sfc_as_cols(names = c("Longitude","Latitude")) %>% 
  st_set_geometry(NULL) %>% 
  left_join(.,raw) %>% 
  filter(!is.na(vals)) %>%
  group_by(station, Longitude, Latitude, datetime, vals) %>%
  dplyr::summarise(daily_TSS = mean(tss, na.rm = T)) %>%
  filter(!is.na(daily_TSS)) %>%
  mutate(week = week(datetime)) %>% 
  st_as_sf(.,coords = c("Longitude","Latitude"), crs = 3857) %>% 
  st_transform(crs = 4326) %>% 
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL) 

save(vcr_wq_modis, file = here::here("data/vcr_wq_modis.rdata"))

# library(leaflet)

# station_geom <- 
#   wq_summarized %>% 
#   dplyr::select(station) %>% 
#   distinct() 
# 
# leaflet() %>% 
#   addTiles() %>% 
#   addMarkers(data = station_geom, popup = ~station)

tss_vis_vcr <- 
  wq %>% 
  filter(!is.na(tss),
         tss < 200, tss > 0)

tss_vis_cb <- 
  cb_wq %>% 
  mutate(Longitude = round(Longitude, 4),
         Latitude = round(Latitude, 4),
         datetime = as.Date(SampleDate, format = "%m/%d/%Y")) %>% 
  filter(Parameter %in% c("TSS", "SALINITY"),
         Depth <= 20) %>%
  # shrink data
  dplyr::select(Station, datetime, Longitude, Latitude, Parameter, 
                MeasureValue, SampleDate, sal_prof = CBSeg2003, Depth, SampleReplicateType, Layer) %>%
  
  # pivot variables to columns
  pivot_wider(., names_from = Parameter, values_from = MeasureValue, values_fn = mean) %>% 
  filter(!is.na(TSS),
         TSS < 200)

ggplot(tss_vis_vcr) +
  geom_density(aes(tss)) +
  geom_density(data = tss_vis_cb, aes(TSS))
