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

temp_data %>%
  mutate(date = as.Date(paste0(yearmonth,"01"), "%m.%Y%d")) %>%
  group_by(date) %>%
  summarise(tmean = mean(tmean)) %>%
  ungroup() %>%
  filter(date>as.Date("2012-05-01"),
         date<as.Date("2013-06-01")) %>%
  pull(tmean) %>%
  mean()

temp_data %>%
  mutate(date = as.Date(paste0(yearmonth,"01"), "%m.%Y%d")) %>%
  group_by(date) %>%
  summarise(tmean = mean(tmean)) %>%
  ungroup() %>%
  filter(date>as.Date("2013-05-01"),
         date<as.Date("2014-06-01")) %>%
  pull(tmean) %>%
  mean()



ptemp <- left_join(tmax_data, tmin_data) %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  # mutate(tmean = (tmax + tmin)/2,
  #        yearmonth = str_sub(date, 4,10)) %>%
  # group_by(yearmonth) %>%
  # summarise(tmean = mean(tmean)) %>%
  # ungroup() %>%
  # mutate(date = as.Date(paste0(yearmonth, "01"),"%m.%Y%d")) %>%
  filter(date > as.Date("2012-04-01"),
         date < as.Date("2015-01-01")) %>%
  mutate(month = lubridate::month(date),
         year = lubridate::year(date)) %>%
  pivot_longer(cols = c(tmax, tmin)) %>%
ggplot() +
  geom_point(aes(x=date,y=value, color = name), alpha=0.1) +
  geom_hline(yintercept = 30, lty=3) +
  geom_hline(yintercept = -10, lty=3) +
  geom_smooth(aes(x=date,y=value, color = name), se=F) +
  geom_labelvline(xintercept = as.Date("2013-04-29"), col="grey40", label = "CRP", hjust=.95, lty=2) +
  geom_labelvline(xintercept = as.Date("2014-05-01"), col="grey40", label = "CRP", hjust=0.05, lty=2)+
  theme_classic() +
  theme(legend.position =c(0,0),
        legend.justification = c(0,0),
        legend.background = element_rect(fill=NA),
        legend.title = element_blank())

ggsave(plot=ptemp, filename = "figs/tempsummary.png", width = 7, height=3, bg="white")


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
