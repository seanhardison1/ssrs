library(XML)
library(anytime)
library(magrittr)
library(tidyverse)

g1 <- xmlParse(file = here::here("data/ffs_modis_rrs_645_2000_2010.kml"))
g2 <- xmlParse(file = here::here("data/ffs_modis_rrs_645_2010_2020.kml"))
xml_data_1 <- xmlToList(g1) 
xml_data_2 <- xmlToList(g2) 

# parsing XML
ffs_modis_reflec <- NULL
for (j in c(1,2)){
  message(j)
  xml_data <- get(paste0("xml_data_",j))
  
  for (i in 1:length(xml_data$Document)){
    unlisted <- unlist(xml_data$Document[i]$Placemark)
    
    locs <- unlisted[names(unlisted) == "ExtendedData.Data..attrs.name"]
    locs <- locs[c(-1, -(length(locs)), (-(length(locs)) + 1))]
    img_df <- tibble(loc = locs) %>% 
      separate(.,col = loc, into = c("Longitude","Latitude"), sep = "_", remove = T) %>% 
      mutate(Longitude = as.numeric(str_remove(Longitude, "X\\."))/-10000,
             Latitude = as.numeric(Latitude)/10000)
    
    vals <- unlisted[names(unlisted) == "ExtendedData.Data.value"]
    vals <- as.numeric(vals[c(-1, -(length(vals)), (-(length(vals)) + 1))])
    
    dt <- anytime(as.numeric(unlisted[length(unlisted) - 2])/1000)
    
    img_df %<>% 
      mutate(vals = vals,
             datetime = dt)
    
    assign("ffs_modis_reflec", rbind(ffs_modis_reflec, img_df))
  }
}


save(ffs_modis_reflec, file = here::here("data/ffs_modis_reflectance_ts.rdata"))

vcr1 <- xmlParse(file = here::here("data/vcr_modis_rrs_645_2000_2020.kml"))
xml_data <- xmlToList(vcr1) 

# parsing XML
vcr_reflec <- NULL
for (i in 1:length(xml_data$Document)){
  unlisted <- unlist(xml_data$Document[i]$Placemark)
  
  locs <- unlisted[names(unlisted) == "ExtendedData.Data..attrs.name"]
  locs <- locs[c(-1, -(length(locs)), (-(length(locs)) + 1))]
  img_df <- tibble(loc = locs) %>% 
    separate(.,col = loc, into = c("Longitude","Latitude"), sep = "_", remove = T) %>% 
    mutate(Longitude = as.numeric(str_remove(Longitude, "X\\."))/-10000,
           Latitude = as.numeric(Latitude)/10000)
  
  vals <- unlisted[names(unlisted) == "ExtendedData.Data.value"]
  vals <- as.numeric(vals[c(-1, -(length(vals)), (-(length(vals)) + 1))])
  
  dt <- anytime(as.numeric(unlisted[length(unlisted) - 2])/1000)
  
  img_df %<>% 
    mutate(vals = vals,
           datetime = dt)
  
  assign("vcr_reflec", rbind(vcr_reflec, img_df))

}

save(vcr_reflec, file = here::here("data/vcr_reflectance_ts.rdata"))
