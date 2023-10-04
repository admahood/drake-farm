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
  soil_moisture_pre_djf_30cm +
  soil_moisture_pre_ma_30cm +
  soil_moisture_pre_son_30cm + 
  soil_moisture_post_mj_30cm +
  soil_moisture_post_jas_30cm +
  soil_moisture_post_ond_30cm + 
  soil_temp_post_jas + 
  soil_temp_post_ond + 
  soil_temp_post_mj + 
  soil_temp_pre_djf + 
  soil_temp_pre_ma + 
  soil_temp_pre_son + 
  twi +
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

df_crp <- data.frame(group = c("pasm", "navi", "bocu", "mesa", "scsc","atca", "bogr", "pavi"),
                    seeding_intensity = c(2.0,0.8, 0.7, 0.4, 0.4, 0.3, 0.2, 0.2))

traits <- data.frame(group = colnames(Y)) %>%
  left_join(heights) %>%
  left_join(df_crp) %>%
  replace_na(list(seeding_intensity = 0)) %>%
  tibble::column_to_rownames("group") %>%
  na.omit() %>%
  mutate(seeded = ifelse(seeding_intensity >0, "yes", "no"))

t_formula <- ~ height + 
  introduced +
  perennial +
  graminoid +
  rhizomatous +
  pp +
  seeded

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

day <- format(Sys.time(), "%b_%d_%Y")

nChains = 4
run_type = "rr"
run_type = "mid"
run_type = "oblas"
run_type = "test"
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
  hmsc_file <- paste0("data/hmsc/hmsc_probit_subplot_mid_",day,".Rda")
}
if (run_type == "rr"){
  
  thin = 2000
  samples = 1000
  transient = ceiling(thin*samples*.5)
  hmsc_file <- paste0("data/hmsc/hmsc_probit_subplot_rr_",day,".Rda")
}
if (run_type == "rrp"){
  nChains = 4
  thin = 200
  samples = ceiling(4000/nChains) 
  transient = ceiling(thin*samples*.5)
  total_iterations <- ((thin*samples) + transient)*nChains
  iterations_per_chain <- total_iterations/nChains
  hmsc_file <- paste0("data/hmsc/hmsc_probit_subplot_rrp_",day,".Rda")
}
if (run_type == "oblas"){
  nChains = 1
  thin = 1000
  samples = 1000
  transient = ceiling(thin*samples*.5)
  hmsc_file <- paste0("data/hmsc/hmsc_probit_subplot_oblas_",thin, "_",day,"a.Rda")
}


t0 <- Sys.time();print(t0)
if(!dir.exists("data/hmsc")) dir.create("data/hmsc")
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
  print(Sys.time()-t0)
  save(m, file=hmsc_file)
}else{load(hmsc_file)}

# plotting =======================================================
source("R/hmsc_plotting_functions.R")
# load("data/hmsc/hmsc_probit_subplot_rrp_Feb_14.Rda")
library(RColorBrewer)

ggplot_ess(m, omega=T)
ggplot_convergence(m,omega = T, gamma=T) %>%
  ggsave(plot=.,filename = "figs/convergence.png", width=5, height=5)

ggplot_fit(m, which="named")%>%
  ggsave(plot=.,filename = "figs/r2.png", width=5, height=5)
  
vp_cols <- c(brewer.pal(12, "Set3"),
             brewer.pal(9, "Set1"),
             "black", "white","burlywood4")

ggplot_vp(m, cols = vp_cols) %>%
  ggsave(filename = "figs/vp.png", plot=., width=10, height=10, bg="white")

ggplot_beta(m, grouping_var = "introduced") %>%
  ggsave(filename = "figs/beta.png", plot=., width=10, height=10, bg="white")

ggplot_gamma(m)

lut_gensp <- c('BROM', 'BASC', 'SATR', 'PASM', 'BOCU', 'MESA', 'BRAS', 
               'NAVI', 'LACT', 'FORB', 'CIAR', 'COAR', 'PAVI', 'SCSC', 
               'ATCA', 'BOGR', 'BRIN', 'SECE')
names(lut_gensp) <- colnames(m$Y)

ggplot_omega(m,lut_gensp=lut_gensp,hc.method = "centroid")
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

