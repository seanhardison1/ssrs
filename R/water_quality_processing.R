library(tidyverse)
library(lubridate)
library(tsibble)
library(sf)

wq <- read_csv(here::here('data/WQ_Integrated.csv'), skip = 21)

wq_summarized <- wq %>% 
  dplyr::select(station,
                scttemp, 
                dotemp,
                chlor_a,
                secchi,
                do, 
                measureDate,
                Latitude, 
                Longitude,
                pom,
                tss,
                tssqual,
                phaeopigments_a,
                phaeopigments_aqual,
                depth) %>% 
  filter(!is.na(Latitude) , !is.na(Longitude), (tssqual != "Q" | is.na(tssqual))) %>% 
  mutate(Longitude = ifelse(Longitude > 0, Longitude * -1, Longitude),
         phaeopigments_a = ifelse(phaeopigments_aqual == "Q", NA, phaeopigments_a)) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude")) %>% 
  group_by(station, measureDate) %>% 
  dplyr::summarise_at(vars(scttemp:do, pom, tss, phaeopigments_a, depth), mean, na.rm = T) 

save(wq_summarized, file = here::here("data/vcr-wq.rdata"))
