# compare sp models to microclima mods

source("R/drake_data_prep.R")

mc_files <- Sys.glob("data/microclima*.tif")

soil_temp_summary_jf <- sentek %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "probe_depth",
               values_to = "soil_temp_c") %>%
  mutate(probe = str_sub(probe_depth, 1,2),
         depth = str_sub(probe_depth, 4,7)) %>%
  filter(depth == "30cm",
         `Date/Time` > as.Date("2012-01-01"),
         `Date/Time` < as.Date("2014-12-01")) %>%
  mutate(year = lubridate::year(`Date/Time`),
         month = lubridate::month(`Date/Time`),
         day = lubridate::day(`Date/Time`),
         hour = lubridate::hour(`Date/Time`)) %>%
  filter(month == 1 | month == 2,
         ) %>%
  group_by(year, probe) %>%
  summarise(max_t = max(soil_temp_c),
            min_t = min(soil_temp_c),
            mean_t = mean(soil_temp_c))%>% 
  pivot_wider(names_from = "year", id_cols = "probe",
              values_from =c("max_t", "min_t", "mean_t"))

d <- sentek_locations %>%
  dplyr::rename(probe = Probe) %>%
  left_join(soil_temp_summary_jf) %>%
  mutate(terra::extract(x=terra::rast(mc_files[1]), 
                 y=terra::vect(.))) %>%
  mutate(terra::extract(x=soil_temp_rasts_13,
                        y=terra::vect(.))) %>%
  mutate(terra::extract(x=soil_temp_rasts_14,
                        y=terra::vect(.)))


ggplot(d, aes(x = mean_t_2013, y=jf_13_tmean)) +
  geom_point()

ggplot(d, aes(x = pre_seed_jf_temp_c, y=jf_13_tmean)) +
  geom_point() 

ggplot(d, aes(x = mean_t_2013, y=pre_seed_jf_temp_c)) +
  geom_point()


temp_compare <- ggpubr::ggarrange(
soil_temp_rasts_13$pre_seed_jf_temp_c %>%
  as.data.frame(xy=TRUE) %>%
  ggplot(aes(x=x,y=y,fill=pre_seed_jf_temp_c)) +
  geom_raster()+
  scale_fill_viridis(name = "Soil Temp (C)")+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Soil Temperatures","Spatial Process Model")+
  theme_classic()+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill= "white", color = "black"))
,
terra::rast(mc_files[1])[[3]]%>%
  as.data.frame(xy=TRUE) %>%
  ggplot(aes(x=x,y=y,fill=jf_13_tmean)) +
  geom_raster()+
  scale_fill_viridis(name = "Air Temp (C)")+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Near-Surface (10cm) Air Temperatures","Microclimate Model")+
  theme_classic()+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill= "white", color = "black"))
,
nrow=2)

ggsave(temp_compare, filename = "figs/temp_compare.png", height = 10, width = 8)
