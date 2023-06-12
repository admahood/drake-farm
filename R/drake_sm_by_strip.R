# looking at soil moisture among strips
library(tidyverse)
library(sf)
# sensor locations ================

sl <- st_read("data/past_data/Sentek probe locations/sentek_probes.shp")
plot(sl[,1], add=TRUE)

ggplot(sl) +
  geom_sf_label(aes(label=Probe))

# soil_moisture data ==========================================================
a_es_drake_e <- readxl::read_xlsx(
  "data/past_data/E.wc.hr.2002-18_Nezat.xlsx",col_types=c("date",rep("numeric",16))) %>%
  dplyr::select(-numericdate) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,1,2),
         depth = str_sub(variable,4,6) %>% str_c("cm")) %>%
  dplyr::select(-variable) %>%
  filter(value > -9999); a_es_drake_e

# ggplot(a_es_drake_e, aes(x=value)) +
#   geom_histogram()

a_es_drake_n <- readxl::read_xlsx(
  "data/past_data/ESDrakeN_2002_2017_daily.xlsx") %>%
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
  dplyr::rename(timestamp = 1); a_es_drake_n

a_es_drake_s <- readxl::read_xlsx(
  "data/past_data/ESDrakeS_2002_2017_daily.xlsx") %>%
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
  dplyr::rename(timestamp = 1); a_es_drake_s

sm_full <- bind_rows(a_es_drake_n, a_es_drake_s) %>%
  mutate(strip = ifelse(probe %in% c("C1", "C2", "C3", "C4", 
                                     "A1" , "A2", "A3", "A4"), "herb", "shrub")) %>%
  filter(depth %in% c("30cm", "60cm", "90cm"))

ggplot(sm_full, aes(x=timestamp, y=value, color = strip)) +
  geom_line(aes(group = probe)) +
  geom_smooth(se=F) +
  facet_wrap(~depth, nrow = 3) +
  geom_vline(xintercept = as.Date("2013-05-15") %>% lubridate::as_datetime())
