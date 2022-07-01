# initial processing

library(tidyverse)
library(sf)

# this is the elevation DEM grid
grid <- st_read("data/grid.gpkg")

# previous sampling locations
c_locations <- st_read("data/redrakesoilsampling/DrakeCarbon198loc.shp")

# selecting only grid cells that overlap with previously sampled soil C locations

sampled_grid <- grid %>%
  mutate(n = rowSums(st_intersects(grid, c_locations, sparse = F))) %>%
  filter(n > 0) %>%
  st_join(c_locations)

sampled_centroids <- sampled_grid %>%
  st_centroid()

st_write(sampled_grid, "data/sampled_grid.gpkg", delete_dsn = TRUE)
st_write(sampled_centroids, "data/sampled_centroids.gpkg", delete_dsn = TRUE)

st_coordinates(sampled_centroids) %>%
  as_tibble() %>%
  write_csv("data/centroid_utms.csv")

# then make a georeferenced PDF in QGIS

# the_chosen_plots <- sample(grid$dem2018, nrow(grid)/10)
# 
# sampled_grid <- grid %>%
#   filter(dem2018 %in% the_chosen_plots) 
# 
# sampled_centroids <-
#   sampled_grid %>%
#   st_centroid()
# 
# st_write(sampled_grid, "data/sampled_grid.gpkg", delete_dsn = TRUE)
# st_write(sampled_centroids, "data/sampled_centroids.gpkg", delete_dsn = TRUE)


# now we have the carbon sampling locations

