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
  soil_texture +
  bare + 
  strip_type +
  total_n_top_15cm_2012 

XFormula_spei <- ~ 
  march_spei06 +
  # march_spei12 +
  march_spei24 +
  june_spei06 +
  # june_spei12 +
  june_spei24 +
  fa +
  twi +
  soil_texture +
  bare + 
  slope +
  strip_type+
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

mod_spei = Hmsc(Y = Y, 
           XData = XData, 
           XFormula = XFormula_spei,
           distr="probit",
           studyDesign = studyDesign,
           TrData = traits,
           TrFormula = t_formula,
           ranLevels = list("plot" = rLpl, "strip_number" = rLsn))

day <- format(Sys.time(), "%b_%d")

nChains = 4
run_type = "rolls_royce"
run_type = "mid"
if (run_type == "test"){
  #with this option, the vignette evaluates in ca. 1 minute in adam's laptop
  thin = 1
  samples = 100
  transient = ceiling(thin*samples*.5)
  hmsc_file <- "data/hmsc/hmsc_probit_subplot_test.Rda"

}
if (run_type == "mid"){
  thin = 10
  samples = 1000
  transient = ceiling(thin*samples*.5)
  hmsc_file <- paste0("data/hmsc/hmsc_probit_subplot_mid_spei",day,".Rda")
}
if (run_type == "rolls_royce"){
  
  thin = 500
  samples = 1000
  transient = ceiling(thin*samples*.5)
  hmsc_file <- paste0("data/hmsc/hmsc_probit_subplot_rr_spei",day,".Rda")
}



t0 <- Sys.time()
dir.create("data/hmsc")
if(!file.exists(hmsc_file)){
  m = sampleMcmc(mod, 
                 thin = thin,
                 samples = samples,
                 transient = transient,
                 # useSocket = FALSE,
                 nChains = nChains,
                 nParallel = nChains)
  print(Sys.time()-t0)
  save(m, file=hmsc_file)
}else{load(hmsc_file)}

t0 <- Sys.time()
dir.create("data/hmsc")
if(!file.exists(hmsc_file)){
  mspei = sampleMcmc(mod_spei, 
                 thin = thin,
                 samples = samples,
                 transient = transient,
                 # useSocket = FALSE,
                 nChains = nChains,
                 nParallel = nChains)
  print(Sys.time()-t0)
  save(mspei, file=hmsc_file)
}else{load(hmsc_file)}

# plotting =======================================================
source("R/hmsc_plotting_functions.R")
load("data/hmsc/hmsc_probit_subplot_rrJan_27.Rda")
library(RColorBrewer)
ggplot_convergence(m,omega = T, gamma=T) %>%
  ggsave(plot=.,filename = "figs/convergence.png", width=5, height=5)
ggplot_convergence(mspei,omega = T, gamma=T)

ggplot_fit(m, which="named")%>%
  ggsave(plot=.,filename = "figs/r2.png", width=5, height=5)
  
ggplot_fit(mspei, which="named")

vp_cols <- c(brewer.pal(12, "Set3"),
             brewer.pal(9, "Set1"),
             "black", "white","burlywood4")

ggplot_vp(m, cols = vp_cols) %>%
  ggsave(filename = "figs/vp.png", plot=., width=10, height=10, bg="white")
ggplot_vp(mspei, cols = vp_cols) %>%
  ggsave(filename = "figs/vp_spei.png", plot=., width=10, height=10, bg="white")
ggplot_beta(m, grouping_var = "introduced")%>%
  ggsave(filename = "figs/beta.png", plot=., width=10, height=10, bg="white")
ggplot_beta(mspei, grouping_var = "introduced")
ggplot_vp(mspei, cols = vp_cols)
ggplot_gamma(m)
ggplot_gamma(mspei)

ggplot_omega(m, hc_method = "centroid", dots=F)%>%
  ggsave(filename = "figs/omega.png", plot=., width=15, height=15, bg="white")

# big multipanel ==================

ggarrange(p_beta, vp +
            theme(plot.background = element_rect(color="black")) ,
          diag_all, nrow=3, ncol=1, labels = c("(a)","(b)","(c)")) %>%
  ggarrange(pcor1, widths = c(1,2.5), labels = c("", "(d)")) %>%
  ggsave(., filename = "figs/multipanel_jsdm.png", height=10, width=15.45, bg="white")

ggarrange(p_beta, vp +
            theme(plot.background = element_rect(color="black")) ,
          diag_all, 
          widths = c(1.5,1.5,1),
          nrow=1, ncol=3, labels = c("(a)","(b)","(c)")) %>%
  ggarrange(pcor1, heights = c(1,2), labels = c("", "(d)"),nrow=2) %>%
  ggsave(., filename = "figs/multipanel_jsdm_vertical.png", height=16.45, width=11, bg="white")



#smaller mulipanel

ggarrange(p_beta, 
          diag_all, 
          nrow=2, ncol=1, labels = c("(a)","(b)")) %>%
  ggarrange(vp + theme(plot.background = element_rect(color="black")),
            nrow=1, ncol=2, widths = c(1,1.5),labels = c("", "(c)")) %>%
  ggsave(., filename = "figs/multipanel_jsdm_nocor.png", height=7, width=12, bg="white")

# gradients ===================

peci_gradient = constructGradient(m, focalVariable = "peci_sub_cv")


predY_peci = predict(m, XData=peci_gradient$XDataNew, 
                     studyDesign=peci_gradient$studyDesignNew, 
                     ranLevels=peci_gradient$rLNew, expected=TRUE)

plotGradient(m, peci_gradient, pred=predY_peci, measure="S")

n_runs <- nChains*samples

pred_df_grazing <- do.call("rbind", predY_peci) %>%
  as_tibble() %>%
  mutate(peci_cover = rep(peci_gradient$XDataNew$peci_sub_cv,
                          n_runs),
         run = rep(1:n_runs,each=20)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -c(peci_cover,run)) %>%
  left_join(prevalence) %>%
  filter(prevalence > 1) %>%
  left_join(spp_list, by = c("Species" = "species")) %>%
  arrange(habit, desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) 

grazing_native <- pred_df_grazing  %>%
  ggplot(aes(x=peci_cover, y=cover, color = habit)) +
  geom_line(alpha = 0.03, aes(group=run), key_glyph="rect")+
  facet_wrap(~Species_f, ncol=7)+
  xlab("PECI Cover") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = c(.9,0.05),
        legend.direction = "horizontal",
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.title = element_blank())

ggsave(grazing_native, filename = "figs/gradient_peci_cover_subplot.png", width = 20, height = 15)



# ggmcmc
mc <- convertToCodaObject(m)
sbeta<-ggs(mc$Beta)
sgamma <- ggs(mc$Gamma)
sv <- ggs(mc$V)


ggs_caterpillar(sv)

ggs_caterpillar(sgamma)
ggs_crosscorrelation(sgamma)%>%
  ggsave(filename = "figs/ccg.png", plot=., 
         width=20, height=20, bg="white")

ggs_crosscorrelation(sv %>%
                       filter(str_sub(Parameter, 1, 22) == 
                                "V[soil_moisture_pre_jf"))%>%
  ggsave(filename = "figs/ccv.png", plot=., 
         width=10, height=10, bg="white")

sbeta_x<-sbeta%>%
  separate(col = "Parameter", into = c("var","xx1", "sp","xx2"),sep = " ",
           remove = F, convert=T,extra = "drop", fill = "warn") %>%
  dplyr::select(-xx1,-xx2)

sbeta_sm <- filter(sbeta, str_sub(Parameter, 1, 9) == "B[soil_mo" )

ggs_caterpillar(sbeta_x %>%
                  filter(var == "B[soil_moisture_pre_jf_30cm")) %>%
  ggsave(filename = "figs/caterpillar_smjf.png", plot=., 
         width=7, height=5, bg="white")

ggs_crosscorrelation(sbeta_x %>%
                  filter(sp == "kosc")) %>%
  ggsave(filename = "figs/cc_kosc.png", plot=., 
         width=7, height=5, bg="white")

ggs_caterpillar(sbeta_x %>%
                  filter(sp == "pasm")) %>%
  ggsave(filename = "figs/caterpillar_pasm.png", plot=., 
         width=7, height=5, bg="white")


prh<-ggs_Rhat(sbeta)

prh %>%
  ggsave(filename = "figs/rhat.png", plot=., 
         width=7, height=35, bg="white")
ggsave(filename = "figs/rhat_gamma.png", plot=ggs_Rhat(sgamma), 
       width=7, height=15, bg="white")
ggsave(filename = "figs/rhat_v.png", plot=ggs_Rhat(sv), 
       width=7, height=25, bg="white")
ggs_crosscorrelation(sbeta_sm) %>%
  ggsave(filename = "figs/cc.png", plot=., 
         width=20, height=40, bg="white")

