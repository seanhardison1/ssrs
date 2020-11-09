library(tidyverse)
library(sf)
library(raster)
library(mgcv)
library(anytime)
library(glmmTMB)
library(GGally)
library(ggeffects)
source(here::here("R/sfc_as_cols.R"))
cb_wq <- read.csv(here::here("data/WaterQualityWaterQualityCBSeg2003.csv"))


raw <- read.csv(here::here("data/10yr_modis_Rrs_555.csv")) %>% 
  dplyr::select(-imageID,-.geo) %>%
  gather(loc, Rrs_555, -system.index,-timeMillis) %>% 
  separate(.,col = loc, into = c("Longitude","Latitude"), sep = "_", remove = F) %>% 
  mutate(Longitude = as.numeric(str_remove(Longitude, "X\\."))/-1000000,
         Latitude = as.numeric(Latitude)/1000000,
         datetime = as.Date(anytime(as.numeric(timeMillis)/1000), format = "%m/%d/%Y")) %>% 
  dplyr::select(-timeMillis) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(.,"+init=epsg:3857") 

ggplot(raw) +
  geom_line(aes(x = datetime, y = Rrs_555, color = loc)) +
  guides(color = F)

station_locs <- cb_wq %>% 
  
  # measurements above pycnocline and at surface
  filter(Layer %in% c("S", "AP")) %>% 
  
  # identify salinity regimes and create datetime column
  mutate(sal_prof = ifelse(str_detect(CBSeg2003,"PH"), "PH", "MH"),
         datetime = as.Date(SampleDate, format = "%m/%d/%Y")) %>%
  filter(!is.na(Longitude), between(datetime, min(raw$datetime), max(raw$datetime))) %>%
  filter(Parameter %in% c("TSS", "SALINITY"),
         Depth < 4) %>% 
  
  # shrink data
  dplyr::select(Station, datetime, Longitude, Latitude, Parameter, 
                MeasureValue, SampleDate, sal_prof = CBSeg2003, Depth, SampleReplicateType, Layer) %>%
  
  # pivot variables to columns
  pivot_wider(., names_from = Parameter, values_from = MeasureValue, values_fn = mean) %>% 
 
  # bind with MODIS time series
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(.,"+init=epsg:3857") %>%
  st_snap(.,raw, tolerance = 5000) %>%
  sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
  st_set_geometry(NULL) %>%
  left_join(.,raw %>%
              sfc_as_cols(.,names = c("Longitude","Latitude")) %>%
              st_set_geometry(NULL),
            by = c("datetime","Longitude","Latitude")) %>%
  
  # return daily summary by station
  group_by(sal_prof, Station, datetime, Rrs_555) %>% 
  dplyr::summarise(daily_TSS = mean(TSS, na.rm = T),
                   daily_salinity = mean(SALINITY, na.rm = T))



ggplot() +
  geom_point(data = station_locs, 
             aes(x = Rrs_555, y = daily_TSS)) +
  scale_y_continuous(trans = "log10")

lmdf <- station_locs %>% 
  filter(!is.na(Rrs_555)) %>% 
  mutate(log_TSS = log10(daily_TSS))

lmdf %>% 
  ungroup() %>%  
  dplyr::select(daily_salinity, log_TSS, Rrs_555) %>% 
  ggpairs()

ggplot(data = lmdf) +
  geom_point(aes(x = Rrs_555, y = log_TSS, color = sal_prof)) +
  geom_smooth(aes(x = Rrs_555, y = log_TSS), method = "lm") +
  facet_wrap(~sal_prof, scales = "free")

mod <- glmmTMB(log_TSS ~ Rrs_555 + sal_prof + (1|Station), data = lmdf)
summary(mod)
mod2 <- lm(log_TSS ~ Rrs_555, data = lmdf)
mod3 <- gam(log_TSS ~ s(Rrs_555), data = lmdf )
plot(mod3)
summary(mod3)

AICcmodavg::aictab(list(mod,mod2))

summary(mod)
dream::validate(list(mod))
r555_mod <- ggpredict(mod, type = "random", terms = c("Rrs_555","Station")) %>% 
  as.data.frame()

ggplot(r555_mod) +
  geom_line(aes(x = x, y = predicted, color = group))
