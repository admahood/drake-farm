# drake beta diversity
# check native div vs non-native div
# predicting seeding success based on pre-treamtnet climate
# what variables will be the best to use for pre-treatment climate
# landscape factors that will drive cRP treatment success -
# splitting the strips
# Hmsc analysis
source("R/drake_data_prep.R")

library(tidyverse)
library(Hmsc)
require(snow)
library(ggpubr)
library(ggcorrplot)
library(ggthemes)
library(ggtext)

# import diversity data ========================================================
# veg data at https://docs.google.com/spreadsheets/d/1sYD0lucZ0X81ebDllucSBeCB2rOWE2v6lzRCvGlLLkI/edit?usp=sharing

raw <- read_csv("data/drake_veg_data_2022 - cover_2022(1).csv")
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

# prepping raw cover data
dv<- plant_cover %>%
  dplyr::select(species_code, cover_pct, plot, n_subplots) %>%
  group_by(plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct)/(n_subplots)) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(names_from = species_code, 
              values_from = cover_pct, 
              values_fill = 0) %>%
  arrange(plot) %>%
  tibble::column_to_rownames("plot") 

vegan::diversity(dv,"invsimpson") %>%
  as_tibble(rownames = "plot") %>%
  dplyr::rename(shannon = value) %>%
  mutate(plot = str_sub(plot,3,5) %>% as.numeric()) %>%
  left_join(xdata) %>%
  ggplot(aes(x=strip_type, y=shannon)) +
  geom_boxplot()
    

pa <- dv
pa[pa>0] <- 1

adespatial::beta.div(Y = dv, nperm = 9999) -> bd
adespatial::beta.div(Y = pa, nperm = 9999) -> bdpa

lcbd_df <- tibble(LCBD = bd$LCBD,
                p_LCBD = bd$p.adj,
                plot = rownames(dv)) %>%
  mutate(strip_type = str_sub(plot,7,10))

lcbd_df <- tibble(LCBD = bdpa$LCBD,
                  p_LCBD = bdpa$p.adj,
                  plot = rownames(dv)) %>%
  mutate(strip_type = str_sub(plot,7,10))

summary(lcbd_df)
ggplot(lcbd_df, aes(x=strip_type, y=LCBD)) +
  geom_boxplot() +
  geom_jitter()

scbd_df <- tibble(SCBD = bdpa$SCBD, species_code = colnames(dv)) %>%
  left_join(sp_list)

ggplot(scbd_df, aes(x=SCBD, y=species_code)) +
  geom_bar(stat="identity")

ggplot(scbd_df, aes(x=SCBD, y=introduced)) +
  geom_boxplot() +
  geom_jitter()


