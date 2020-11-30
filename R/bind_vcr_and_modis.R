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

vcr_wq_modis <- NULL
for (i in c(-1,0,1)){
  
  raw2 <- raw %>% 
    mutate(datetime = datetime + i)
  
  vcr_wq_modis2 <- wq_summarized %>% 
    st_set_crs(4326) %>% 
    sfc_as_cols(names = c("Longitude","Latitude")) %>% 
    mutate(Longitude = round(Longitude, 4),
           Latitude = round(Latitude, 4),
           datetime = as.Date(measureDate, format = "%m/%d/%Y")) %>% 
    dplyr::select(tss, station, datetime, Longitude, Latitude, depth) %>% 
    
    filter(!is.na(Longitude),
           station %in% c("SS", "QI", "SSI", "NM", "LC", "MI")) %>% 
    # (depth < 2 | is.na(depth))) %>% 
    st_transform(.,"+init=epsg:3857") %>% 
    st_snap(.,unique_refs, tolerance = 1000) %>% 
    sfc_as_cols(names = c("Longitude","Latitude")) %>% 
    st_set_geometry(NULL) %>% 
    left_join(.,raw2) %>% 
    filter(!is.na(vals)) %>%
    group_by(station, Longitude, Latitude, datetime, depth, vals) %>%
    dplyr::summarise(daily_TSS = mean(tss, na.rm = T)) %>%
    filter(!is.na(daily_TSS)) %>%
    mutate(week = week(datetime)) %>% 
    st_as_sf(.,coords = c("Longitude","Latitude"), crs = 3857) %>% 
    st_transform(crs = 4326) %>% 
    sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
    st_set_geometry(NULL) 
  
  assign('vcr_wq_modis', rbind(vcr_wq_modis2, vcr_wq_modis))
}



save(vcr_wq_modis, file = here::here("data/vcr_wq_modis.rdata"))