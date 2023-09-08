library(tidyverse)
library(sf)

locs <- st_read("projects/drake-farm/data/trait_locations.gpkg")

sampled_already <- c(59, 27, 29, 31, 33, 18, 6,5,15,16,42,40,57)

locs %>%
  filter(!X2012Flag %in% sampled_already) %>%
  st_write("projects/drake-farm/data/trait_locs_friday0728.gpkg")
