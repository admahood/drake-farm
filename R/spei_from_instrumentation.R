library(tidyverse)
library(SPEI)
library(sf)
library(terra)

# data input ===================================================================
dem <- raster::raster("data/dem2018.tif")
demt <- terra::rast("data/dem2018.tif")


# precip =======================================================================
precip_sensors <- readxl::read_xlsx("data/precip_temp_insitu/precip.xlsx", 
                                    sheet = "sensors") %>%
  dplyr::select(-time) %>%
  dplyr::filter(`@H` != "elevation", `@H` != "ID") %>%
  pivot_longer(names_to = "sensor", values_to = "value", -`@H`) %>%
  pivot_wider(names_from = `@H`) %>%
  st_as_sf(coords = c("x", "y"), crs = raster::crs(dem, asText=TRUE)) %>%
  mutate(dem = raster::extract(dem,.)) %>%
  filter(!is.na(dem)); precip_sensors


precip_data <- readxl::read_xlsx("data/precip_temp_insitu/precip.xlsx", 
                                 sheet = "rain") %>%
  pivot_longer(names_to = "sensor", values_to = "precip", -time) %>%
  filter(precip > -9999) %>%
  mutate(yearmonth = str_sub(time, 4,10)) %>%
  group_by(yearmonth, sensor) %>%
  summarise(precip = sum(precip)) %>%
  ungroup() %>%
  mutate(date = as.Date(paste0("01.", yearmonth), "%d.%m.%Y")) %>%
  filter(date > as.Date("2009-12-31"), date < as.Date("2017-01-01"),
         sensor %in% precip_sensors$sensor); precip_data

ppt_joined <- precip_sensors %>%
  right_join(precip_data)

# temp =========================================================================

temp_sensors <- readxl::read_xlsx("data/precip_temp_insitu/temp.xlsx", 
                                    sheet = "sensors") %>%
  dplyr::select(-date) %>%
  dplyr::filter(`@H` != "elevation", `@H` != "ID") %>%
  pivot_longer(names_to = "sensor", values_to = "value", -`@H`) %>%
  pivot_wider(names_from = `@H`) %>%
  st_as_sf(coords = c("x", "y"), crs = raster::crs(dem, asText=TRUE))

temp_sensors_drake <- temp_sensors %>%
  mutate(dem = raster::extract(dem,.)) %>%
  filter(!is.na(dem)); temp_sensors


tmax_data <- readxl::read_xlsx("data/precip_temp_insitu/temp.xlsx", 
                                 sheet = "tmax_tidy") %>%
  pivot_longer(names_to = "sensor", values_to = "tmax", -date) %>%
  filter(tmax > -9999);tmax_data

tmin_data <- readxl::read_xlsx("data/precip_temp_insitu/temp.xlsx", 
                               sheet = "tmin_tidy") %>%
  pivot_longer(names_to = "sensor", values_to = "tmin", -date) %>%
  filter(tmin > -9999);tmin_data

temp_data <- left_join(tmax_data, tmin_data) %>%
  mutate(tmean = (tmax + tmin)/2,
         yearmonth = str_sub(date, 4,10)) %>%
  group_by(sensor, yearmonth) %>%
  summarise(tmean = mean(tmean)) %>%
  ungroup()
tmp <- temp_sensors %>%
  right_join(temp_data) %>%
  na.omit()



# kriging temp
ym <- "03.2014" 

d <- filter(tmp, yearmonth == ym)

# data prep for spatial process model
x = st_coordinates(d)[,c(1,2)]

z = dplyr::select(d, dem) %>%
  st_set_geometry(NULL)
y = d %>%
  filter(yearmonth == ym) %>%
  pull(tmean)

print(paste(ym))
mod <- fields::spatialProcess(x=x, y=y, Z=z)


krik <- fields::Krig(x=x, Y=y, z=z,df = 6)
surface(krik)
print(summary(mod))

elv_df<-dem%>%
  as.data.frame(xy=TRUE)

lat <- elv_df$y
lon <- elv_df$x
elev <- elv_df[,3]
xps <-cbind(lon, lat)

## Predict using the spatial process model fitted
yp = predict(mod,
             x=xps,Z=elev)

# convert to a raster
spmod_rast_sm <- raster::rasterFromXYZ(data.frame(lon=lon,lat=lat,  
                                                  prediction = yp),
                                       crs = raster::crs(elevation)) %>%
  terra::rast()
names(spmod_rast_sm) <- paste0(month_grp, "_", sample_year)
result_sm[[counter]] <- spmod_rast_sm
