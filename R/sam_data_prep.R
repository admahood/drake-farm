# prep for SAM model

source("R/drake_data_prep.R")

raw <- read_csv("data/drake_veg_data_2022 - cover_2022(1).csv")

# new groups - further cinched down
sp_list <- read_csv("data/drake_veg_data_2022 - species_list.csv") %>%
  dplyr::select(species_code = code,group, introduced, perennial, woody, graminoid,
                rhizomatous, pp = photosynthetic_pathway)

surface_cover <- raw %>%
  filter(str_sub(species_code,1,1)=="_") %>%
  group_by(plot) %>%
  mutate(n_subplots = length(unique(subplot)),
         n_subplots = if_else(n_subplots > 6 , 8, 4)) %>%
  ungroup()

plant_cover <- raw %>%
  filter(str_sub(species_code,1,1)!="_") %>%
  mutate(strip_type = str_sub(plot, 7,10)) %>%
  group_by(plot) %>%
  mutate(n_subplots = length(unique(subplot)),
         n_subplots = if_else(n_subplots > 6 , 8, 4)) %>%
  ungroup()

# Hmsc-specific data wrangling =================================================

groups <- read_csv("data/drake_veg_data_2022 - species_list.csv") %>%
  dplyr::select(species_code = code, group)

# prepping raw cover data
dv<- plant_cover %>%
  filter(subplot < 9) %>%
  mutate(plot_sub = str_c(plot, "_", subplot)) %>%
  dplyr::select(species_code, cover_pct, plot_sub) %>%
  dplyr::left_join(groups) %>%
  group_by(plot_sub, group) %>%
  summarise(cover_pct = sum(cover_pct)) %>%
  ungroup() %>%
  pivot_wider(names_from = group, values_from = cover_pct, values_fill = 0) %>%
  arrange(plot_sub) %>%
  tibble::column_to_rownames("plot_sub")

# converting to occurrence matrix
Y <- dv %>%
  as.matrix
Y[Y>0] <-1