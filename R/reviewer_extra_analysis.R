# soil moisture analysis by treatment for reviewer #1

# A potential confounding factor with the study design is that the first 
# planting was sown into winter wheat that was harvested and the second planting
# was sown into winter wheat that was not harvested due to low grain yield. The 
# second planting had more soil moisture, which was attributed to climate. But 
# could the residual biomass of unharvested winter wheat have contributed to 
# higher soil moisture in addition to climate? To this end, the results should 
# be interpreted in the context of different antecedent environmental conditions 
# or climate and soil moisture conditions throughout the manuscript (as stated
# in the first sentence of the discussion).

library(tidyverse)
library(sf)

# get probe locations and mgmt units
mgmt_strips <- st_read("data/strips/") %>%
  mutate(mgmt = c("seeded_2013", "seeded_2014"))
sentek_locations <- st_read(
  "data/past_data/Sentek probe locations/sentek_probes.shp")
ggplot(mgmt_strips[2,]) +
  geom_sf(aes(fill = mgmt)) +
  geom_sf(data=sentek_locations)

key_mgmt_probe <- st_intersection(sentek_locations, mgmt_strips[2,]) %>%
  st_set_geometry(NULL) %>%
  dplyr::select(Probe, mgmt) %>%
  right_join(st_set_geometry(sentek_locations, NULL)) %>%
  replace_na(list(mgmt = "seeded_2013")) %>%
  dplyr::select(probe = Probe, mgmt)

# sentek input
e <- readxl::read_xlsx(
  "data/past_data/E.wc.hr.2002-18_Nezat.xlsx",col_types=c("date",rep("numeric",16))) %>%
  # filter(timestamp > as.Date("2012-01-01"),
  #        timestamp < as.Date("2014-01-01")) %>%
  dplyr::select(-numericdate) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,1,2),
         depth = str_sub(variable,4,6) %>% str_c("cm")) %>%
  dplyr::select(-variable) %>%
  filter(value > -9999)

n <- readxl::read_xlsx(
  "data/past_data/ESDrakeN_2002_2017_daily.xlsx") %>%
  # filter(`Date Time` > as.Date("2012-01-01"),
  #        `Date Time` < as.Date("2014-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999)%>%
  dplyr::select(-variable, -measure) %>%
  dplyr::rename(timestamp = 1)

s <- readxl::read_xlsx(
  "data/past_data/ESDrakeS_2002_2017_daily.xlsx") %>%
  # filter(`Date Time` > as.Date("2012-01-01"),
  #        `Date Time` < as.Date("2014-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999) %>%
  dplyr::select(-variable, -measure) %>%
  dplyr::rename(timestamp = 1)

lines_df <- data.frame(timestamp = c(as_datetime("2014-05-01"),
                                     as_datetime("2013-05-01")),
                       mgmt = c("seeded_2014",
                                "seeded_2013"),
                       value = NA)
wheat <- data.frame(timestamp = c(as_datetime("2012-07-15"),
                                  as_datetime("2013-07-15")),
                       mgmt1 = c("wheat_harvest"),
                       mgmt = c("seeded_2013", "seeded_2014"),
                       value = NA)                                    
bind_rows(s,n,e) %>%
  left_join(key_mgmt_probe) %>%
  filter(depth == "60cm") %>%
  mutate(year = lubridate::year(timestamp),
         month = lubridate::month(timestamp)) %>%
  filter(year>2011, year <2019) %>%
  ggplot(aes(x=timestamp, y=value, color = mgmt)) +
  geom_smooth() +
  facet_grid(mgmt~year, scales = "free_x") +
  geom_vline(data = lines_df, aes(xintercept = timestamp,
                                  color = mgmt)) +
  geom_vline(data = wheat, aes(xintercept = timestamp,
                                  color = mgmt1))

# just looking at planting year
bind_rows(s,n,e) %>%
  left_join(key_mgmt_probe) %>%
  filter(depth == "30cm") %>%
  mutate(year = lubridate::year(timestamp),
         month = lubridate::month(timestamp)) %>%
  filter(year== 2013 & mgmt == "seeded_2013"|
         year == 2014 & mgmt == "seeded_2014") %>%
  ggplot(aes(x=timestamp, y=value, color = mgmt)) +
  geom_smooth() +
  facet_wrap(~year, scales = "free_x") +
  geom_vline(data = lines_df, 
             aes(xintercept = timestamp))

