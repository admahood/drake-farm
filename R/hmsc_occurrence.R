# Hmsc analysis
source("R/drake_data_prep.R")

library(tidyverse)
library(Hmsc)
require(snow)
library(ggpubr)
library(ggthemes)
library(ggtext)

# import diversity data ========================================================
# veg data at https://docs.google.com/spreadsheets/d/1sYD0lucZ0X81ebDllucSBeCB2rOWE2v6lzRCvGlLLkI/edit?usp=sharing

raw <- read_csv("data/drake_veg_data_2022 - cover_2022(1).csv")
sp_list <- read_csv("data/drake_veg_data_2022 - species_list.csv") %>%
  dplyr::select(species_code = code, introduced, perennial, woody, annual, graminoid,
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


# creating a trait table for height =============================================

heights <- plant_cover %>%
  group_by(species_code) %>%
  summarise(height = mean(height_cm, na.rm=TRUE)) %>%
  ungroup() %>%
  na.omit() %>%
  left_join(sp_list) %>%
  mutate_if(is.character, as.factor)

# Hmsc-specific data wrangling =================================================

# prepping raw cover data
dv<-plant_cover %>%
  filter(subplot < 9) %>%
  mutate(plot_sub = str_c(plot, "_", subplot)) %>%
  dplyr::select(species_code, cover_pct, plot_sub) %>%
  pivot_wider(names_from = species_code, values_from = cover_pct, values_fill = 0) %>%
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
  mutate(fa = abs(180 - abs(aspect - 225)))
  
XFormula <- ~ bare + 
  slope + #strip_type +
  fa +
  post_seed_jja_temp_c +
  post_seed_son_temp_c +
  pre_seed_jf_temp_c +
  pre_seed_mam_temp_c +
  pre_seed_son_temp_c+
  jf_pre_air_temp_c_tmean +
  mam_pre_air_temp_c_tmean +
  son_pre_air_temp_c_tmean +
  twi +
  soil_unit_name +
  total_n_top_15cm_2012 

traits <- data.frame(species_code = colnames(Y)) %>%
  left_join(heights) %>%
  tibble::column_to_rownames("species_code") %>%
  na.omit()


t_formula <- ~ height + 
  introduced +
  perennial +
  woody +
  annual +
  graminoid +
  rhizomatous +
  pp

studyDesign <- data.frame(plot = as.factor(XData$plot))
rL <- HmscRandomLevel(units = levels(studyDesign$plot))

# The models ===================================================================

mod = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr="probit",
           studyDesign = studyDesign,
           TrData = traits,
           TrFormula = t_formula,
           ranLevels = list("plot" = rL))



nChains = 4
test.run = TRUE
if (test.run){
  #with this option, the vignette evaluates in ca. 1 minute in adam's laptop
  thin = 1
  samples = 100
  transient = ceiling(thin*samples*.5)
  hmsc_file <- "data/hmsc/hmsc_probit_subplot_test.Rda"
}else{

  thin = 100
  samples = 1000
  transient = ceiling(thin*samples*.5)
  hmsc_file <- "data/hmsc/hmsc_probit_subplot.Rda"
}

t0 <- Sys.time()
dir.create("data/hmsc")
if(!file.exists(hmsc_file)){
  m = sampleMcmc(mod, thin = thin,
                 samples = samples,
                 transient = transient,
                 adaptNf = rep(ceiling(0.4*samples*thin),1),
                 nChains = nChains,
                 nParallel = nChains)
  print(Sys.time()-t0)
  save(m, file=hmsc_file)
}else{load(hmsc_file)}

# getting the posteriors =======================================================

mpost <- convertToCodaObject(m)
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
VP <- computeVariancePartitioning(m)

# model convergence, diagnostics ===============================================

# psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf%>%
#   as_tibble() %>% dplyr::rename(psrf_v = `Point est.`)
# 
# ess.v <- effectiveSize(mpost$V)%>%
#   as_tibble() %>% dplyr::rename(ess_v = value)

ess.beta <- effectiveSize(mpost$Beta) %>%
  as_tibble() %>% 
  dplyr::mutate(variable = "Effective Sample Size")
psrf.beta <- gelman.diag(mpost$Beta, multivariate=FALSE)$psrf%>%
  as_tibble() %>% 
  dplyr::rename(value = `Point est.`) %>%
  dplyr::mutate(variable = "Gelman Diagnostic")




diag_all<-bind_rows(ess.beta, psrf.beta) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  theme_classic() +
  facet_wrap(~variable, scales="free") +
  ggtitle("Convergence Diagnostics")+
  theme(plot.background = element_rect( color ="black"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 1, face = "bold"))



ggsave(diag_all,filename = paste0("figs/gelman_ess.png"), 
       width = 5.5, height=3.5, bg="white")
MF$TjurR2%>% mean(na.rm=T)
# explanatory power

ggarrange(
  ggplot(as.data.frame(MF),aes(x=(RMSE))) + geom_histogram(),
  ggplot(as.data.frame(MF),aes(x=(TjurR2))) + geom_histogram(),
  ggplot(as.data.frame(MF),aes(x=(AUC))) + geom_histogram())

# plot the variance partitioning ===============================================

mf_df <- data.frame(Species = colnames(m$Y),
                    R2 = MF$TjurR2,
                    AUC = MF$AUC,
                    RMSE = MF$RMSE) %>%
  left_join(prevalence)
mean(mf_df%>% filter(prevalence>7) %>% pull(R2), na.rm=T)
ggplot(mf_df, aes(x=prevalence, y=R2)) +
  geom_point()

sbquants <- summary(mpost$Beta)$quantiles %>%
  as_tibble(rownames = "variable") %>% 
  mutate(sign = `2.5%` * `97.5%`) %>%
  filter(sign>0) %>%
  separate(variable,
           into = c("variable", "species"),
           sep = ",") %>%
  mutate(variable = str_sub(variable, 3,nchar(variable)-5),
         species = str_sub(species, 2,nchar(species)-6) %>% trimws) %>%
  filter(variable!= "(Intercept)") %>%
  dplyr::select(variable,species,`2.5%`,`50%`,`97.5%`) %>%
  arrange(variable)


vp_df <- VP$vals%>%
  as_tibble(rownames = "variable") %>%
  pivot_longer(cols=names(.)[2:ncol(.)], 
               names_to = "Species", 
               values_to = "value") %>%
  left_join(prevalence) %>%
  na.omit()

vp_summary <- vp_df %>%
  group_by(variable) %>%
  summarise(value_pct = mean(value) * 100) %>%
  ungroup() 



# 
# vp_order <- vp_df %>% filter(variable == "Random: sample") %>%
#   filter(origin=="I") %>%
#   left_join(prevalence) %>%
#   arrange(prevalence, origin) %>%
#   mutate(Species_f = factor(Species, levels = .$Species)) %>%
#   dplyr::select(Species, Species_f, origin) %>%
#   rbind(vp_order_n)# %>%
#   #left_join(mf_df)


vp_order <- vp_df %>%
  filter(variable == "bare") %>%
  arrange(value) %>%
  mutate(Species_f = factor(Species, levels = .$Species)) %>%
  dplyr::select(Species, Species_f) 

vp <- left_join(vp_df, vp_order) %>% 
  ggplot(aes(x=value,y=Species_f, fill = variable)) +
  geom_bar(stat="identity", color="black")+
  theme_classic() +
  ylab("Species") +
  xlab("Proportion of Variance Explained") +
  # scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "right",
        legend.text = element_markdown(),
        legend.title = element_blank(),
        # legend.justification = c(1,0),
        legend.background = element_rect(color="black")) +
  ggtitle("Variance Partitioning")


ggsave(vp, filename=paste0("figs/variance_partitioning.png"),
       height = 11.5, width = 12)

# Environmental filters ==================================================

postBeta <- getPostEstimate(m, parName = "Beta")

covNamesNumbers <- c(TRUE, TRUE)
covNames = character(m$nc)
for (i in 1:m$nc) {
  sep = ""
  if (covNamesNumbers[1]) {
    covNames[i] = paste(covNames[i], m$covNames[i], sep = sep)
    sep = " "
  }
  if (covNamesNumbers[2]) {
    covNames[i] = paste(covNames[i], sprintf("(C%d)", i), sep = sep)
  }
}
covNames

means <- postBeta$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c(covNames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Mean")

supported <- postBeta$support %>% 
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = covNames) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], 
               names_to = "Species", 
               values_to = "Support") %>%
  filter(Support >0.89|Support<0.11,
         env_var != "intercept") %>%
  left_join(means, by = c("env_var", "Species"))%>%
  mutate(sign = ifelse(Mean>0, "+", "-"))%>%
  left_join(vp_order)#

p_beta<-supported %>%
  mutate(env_var = replace(env_var, env_var == "fa", "aspect")) %>%
  ggplot(aes(x=env_var,y=reorder(Species_f,Species), fill = Mean, color = sign)) +
  geom_tile(lwd=.5) +
  theme_pubclean()+
  scale_fill_steps2() +
  scale_color_manual(values = c(("red"), ("blue"))) +
  guides(color = "none")+
  theme(axis.text.x = element_text(angle=45, vjust=1,hjust = 1),
        axis.title = element_blank(),
        legend.position = "left",
        plot.background = element_rect(color="black"),
        plot.title = element_text(hjust = 1, face = "bold")) +
  ggtitle("Environmental Filters")

ggsave(p_beta,filename= paste0("figs/betas_binomial_subplot.png"), 
       bg="white", width=10, height=8)




plotBeta(m, post = postBeta, param = "Support",
         supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))

# traits =======================================================================
trNames = character(m$nt)
trNamesNumbers = c(T,T)
for (i in 1:m$nt) {
  sep = ""
  if (trNamesNumbers[1]) {
    trNames[i] = paste(trNames[i], m$trNames[i], sep = sep)
    sep = " "
  }
  if (trNamesNumbers[2]) {
    trNames[i] = paste(trNames[i], sprintf("(T%d)", i), 
                       sep = sep)
  }
}

postGamma = getPostEstimate(m, parName="Gamma")
plotGamma(m, post=postGamma, param="Support", supportLevel = 0.89, 
          covNamesNumbers = c(T,F), trNamesNumbers = c(T,F), colorLevels = 3)


means_gamma <- postGamma$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c(covNames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Trait", values_to = "Mean")

lut_gamma <- trNames
names(lut_gamma) <- unique(means_gamma$Trait)

supported_gamma <- postGamma$support %>% 
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = covNames) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], 
               names_to = "Trait", 
               values_to = "Support") %>%
  filter(Support >0.89|Support<0.11,
         env_var != "intercept") %>%
  left_join(means_gamma, by = c("env_var", "Trait"))%>%
  mutate(sign = ifelse(Mean>0, "+", "-"),
         Trait = lut_gamma[Trait])

p_gamma<-supported_gamma %>%
  mutate(env_var = replace(env_var, env_var == "fa", "aspect")) %>%
  ggplot(aes(x=env_var,y=(Trait), fill = Mean, color = sign)) +
  geom_tile(lwd=.5) +
  theme_pubclean()+
  scale_fill_steps2() +
  scale_color_manual(values = c(("red"), ("blue"))) +
  guides(color = "none")+
  theme(axis.text.x = element_text(angle=45, vjust=1,hjust = 1),
        axis.title = element_blank(),
        legend.position = "left",
        plot.background = element_rect(color="black"),
        plot.title = element_text(hjust = 1, face = "bold")) +
  ggtitle("Environmental Filters")

# species associations =========================================================
OmegaCor = computeAssociations(m)
supportLevel = 0.89


hmdf_mean <- OmegaCor[[1]]$mean %>%
  as.matrix
hmdf_support <- OmegaCor[[1]]$support %>%
  as.matrix

# avg association strengths
OmegaCor[[1]]$mean %>%
  abs() %>%
  rowSums() %>%
  as_tibble(rownames = "Species") %>%
  arrange(desc(value)) %>%
  left_join(prevalence) %>%
  filter(prevalence>10) %>%
  print(n=20) %>%
  write_csv("figs/residual_correlation_abs.csv")


# switch colors around
pmat <-OmegaCor[[1]]$support
pmat[pmat < (1- supportLevel)] <- 0.99

pcor1<- ggcorrplot::ggcorrplot(hmdf_mean,type = "lower",
                               hc.order = TRUE,hc.method = "single",
                               colors = c("red", "white", "blue"),
                               p.mat = pmat,
                               pch = 20,
                               sig.level = supportLevel,
                               title = "Species Associations") +
  theme(plot.background = element_rect(color="black"),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(color="black"),
        plot.title = element_text(hjust = 1, face = "bold"))

ggsave(pcor1,
       filename=paste0("figs/species_associations.png"),
       bg="white", width=10, height = 10)


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
