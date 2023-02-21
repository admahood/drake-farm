#hmsc 4 slurm

# Hmsc analysis
source("R/drake_data_prep.R")

library(tidyverse)
library(Hmsc)
require(snow)
library(ggpubr)
library(ggcorrplot)
library(ggthemes)
library(ggtext)
# rareness vs dominance?
# talk about persistence after establishment
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

# calculating prevalence
prevalence<- colSums(Y) %>%
  as_tibble(rownames = "Species") %>%
  dplyr::rename(prevalence = value) %>%
  arrange(desc(prevalence))

# getting site data from plant cover df


XData<-left_join(
  plant_cover %>%
    mutate(plot_sub = str_c(plot, "_", subplot)) %>%
    filter(subplot <9) %>%
    group_by(plot_sub) %>%
    summarise(strip_type = first(strip_type),
              plot = first(plot)) %>%
    ungroup() %>%
    mutate(strip_number = str_sub(plot, 12,13))
  ,
  surface_cover %>%
    mutate(species_code = str_remove_all(species_code,"_")) %>%
    mutate(plot_sub = str_c(plot, "_", subplot)) %>%
    filter(subplot <9) %>%
    dplyr::select(species_code, cover_pct, plot_sub) %>%
    pivot_wider(names_from = species_code, values_from = cover_pct, values_fill = 0)
) %>%
  arrange(plot_sub) %>%
  column_to_rownames("plot_sub") %>%
  mutate(plot_raw = str_sub(plot, 3,5) %>% as.numeric) %>%
  left_join(xdata %>% 
              dplyr::rename(plot_raw = plot),
            by = c("plot_raw", "strip_type")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(fa = abs(180 - abs(aspect - 225)),
         strip_number = as.factor(strip_number.x))

# Formula(s) ===================================================================
XFormula_pre <- ~ 
  soil_moisture_pre_jf_30cm + 
  soil_moisture_pre_ma_30cm + 
  soil_moisture_pre_son_30cm + 
  soil_temp_pre_jf + 
  soil_temp_pre_ma + 
  soil_temp_pre_son + 
  air_temp_jf_pre + 
  air_temp_mam_pre + 
  air_temp_son_pre + 
  twi +
  # soil_texture +
  bare + 
  strip_type +
  total_n_top_15cm_2012



# maybe do one with only spei+topo

# creating a trait table for height 

heights <- plant_cover %>%
  left_join(groups) %>%
  group_by(group) %>%
  summarise(height = median(height_cm, na.rm=TRUE)) %>%
  ungroup() %>%
  na.omit() %>%
  left_join(sp_list) %>%
  mutate_if(is.character, as.factor) %>%
  dplyr::select(-species_code) %>%
  group_by(group) %>%
  summarise_all(first) %>%
  ungroup()

traits <- data.frame(group = colnames(Y)) %>%
  left_join(heights) %>%
  tibble::column_to_rownames("group") %>%
  na.omit()


t_formula <- ~ height + 
  introduced +
  perennial +
  woody +
  graminoid +
  rhizomatous +
  pp

studyDesign <- data.frame(plot = as.factor(XData$plot),
                          strip_number = as.factor(XData$strip_number))
rLpl <- HmscRandomLevel(units = studyDesign$plot)
rLsn <- HmscRandomLevel(units = studyDesign$strip_number)

# The models ===================================================================

mod = Hmsc(Y = Y, 
           XData = XData, 
           XFormula = XFormula_pre,
           distr="probit",
           studyDesign = studyDesign,
           TrData = traits,
           TrFormula = t_formula,
           ranLevels = list("plot" = rLpl, "strip_number" = rLsn))

day <- format(Sys.time(), "%b_%d")

nChains = 4
run_type = "test"
run_type = "oblas"
if (run_type == "test"){
  #with this option, the vignette evaluates in ca. 1 minute in adam's laptop
  thin = 1
  samples = 100
  transient = ceiling(thin*samples*.5)
  hmsc_file <- "data/hmsc/hmsc_probit_subplot_test.Rda"
  
}
if (run_type == "oblas"){
  nChains = 4
  thin = 500
  samples = 1000
  transient = ceiling(thin*samples*.5)
  hmsc_file <- paste0("data/hmsc/hmsc_probit_subplot_rr_",day,".Rda")
}


t0 <- Sys.time();print(t0)
if(!dir.exists("data/hmsc"))dir.create("data/hmsc")
if(!file.exists(hmsc_file)){
  
  ttest0 <- Sys.time()
  mtest <- sampleMcmc(mod, samples = 30, nChains = nChains, 
                      nParallel = nChains, verbose=F)
  ttest1 <- Sys.time()
  t_per_iter<-(ttest1 - ttest0)/(30)
  print(paste(t_per_iter, "seconds per iteration"))
  estimated_time <- t_per_iter * (transient*3)
  print(paste((as.numeric(estimated_time)/60)/60, "hours estimated"))
  
  
  m = sampleMcmc(mod, 
                 thin = thin,
                 samples = samples,
                 transient = transient,
                 useSocket = FALSE,
                 nChains = nChains,
                 nParallel = nChains)
  hmsc_file <- str_replace(hmsc_file, run_type, 
                           paste0(run_type, "_",
                                  format(round(Sys.time()- t0)) %>% 
                                    str_replace_all(" ", "_")))
  save(m, file=hmsc_file)
}else{load(hmsc_file)}