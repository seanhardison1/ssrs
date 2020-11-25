library(rgee)
library(tidyverse)
library(cptcity)
library(raster)
library(stars)
library(sf)

ee_Initialize()

cb_loc <- ee$Geometry$Point(c(-76.377094, 37.198504))

s2 <- ee$ImageCollection("COPERNICUS/S2_SR")

cb_poly <- st_read(here::here("data/lower_cb_poly.kml")) %>% 
  sf::st_zm() %>% 
  st_transform(., crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>% 
  st_geometry() %>% 
  sf_as_ee()

getQABits <- function(image, qa) {
  # Convert decimal (character) to decimal (little endian)
  qa <- sum(2^(which(rev(unlist(strsplit(as.character(qa), "")) == 1))-1))
  # Return a single band image of the extracted QA bits, giving the qa value.
  image$bitwiseAnd(qa)$lt(1)
}

s2_clean <- function(img) {
  # Select only band of interest, for instance, B2,B3,B4,B8
  img_band_selected <- img$select("B[2-4|8]")
  
  # quality band
  qa <- img$select("QA60")
  
  # Select pixels to mask
  quality_mask <- getQABits(qa, "110000000000")
  
  # Mask pixels with value zero.
  img_band_selected$updateMask(quality_mask)
}

s2_cb <- s2$
  filterBounds(cb_loc)$
  filter(ee$Filter$lte("CLOUDY_PIXEL_PERCENTAGE", 10))$
  filter(ee$Filter$date('2017-01-01', '2020-09-01'))$
  filter(ee$Filter$calendarRange(6, field = "month"))$
  map(s2_clean)

nimages <- s2_cb$size()$getInfo()
ic_date <- ee_get_date_ic(s2_cb)


Map$centerObject(cb_loc,zoom = 8)
s2_img_list <- list() 
for (index in seq_len(nimages)) {
  py_index <- index - 1
  s2_img <- ee$Image(s2_cb$toList(1, py_index)$get(0))
  s2_img_list[[index]] <- Map$addLayer(
    eeObject = s2_img,
    visParams = list(min = -0.1, max = 0.8, palette = cpt("grass_ndvi", 10)),
    name = ic_date$id[index]
  )
}
Reduce('+', s2_img_list)
