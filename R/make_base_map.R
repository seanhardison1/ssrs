library(sf)
library(ggmap)
library(ggplot2)
library(leaflet)

# Note: You won't be able to use ggmap without first getting an
# API key to query google
# See https://github.com/dkahle/ggmap

dabay <- c(-77,36.5, -75,  39.5)


cbmap_query <- get_map(dabay,source = "google",maptype = "roadmap") 

save(cbmap_query, file = here::here("data/cb_ggmap.rdata"))
