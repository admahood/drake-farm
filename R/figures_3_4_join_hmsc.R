# joining hmsc mods together
# devtools::install_github("hmsc-r/HMSC")
library(Hmsc)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggtext)
library(RColorBrewer)
source("R/hmsc_plotting_functions.R")
mods <- list.files("data/hmsc", full.names=T, pattern = "1000")
ml<- list()
for(i in 1:length(mods)){
  load(mods[i]);ml[[i]] <- m
}
mj <- do.call(c, ml)

# Hm$XData %>% names() -> nn
# nn<-nn[c(8:20)]

ggplot_convergence(mj, gamma=T,omega=T) %>%
  ggsave(plot=.,filename = "figs/convergence.png", width=5, height=5)
ggplot_fit(mj, which = "named")%>%
  ggsave(plot=.,filename = "figs/r2.png", width=5, height=5)

vp_cols <- c(brewer.pal(12, "Set3"),
             brewer.pal(9, "Set1"),
             "black", "white","burlywood4")
ggplot_vp(mj, cols = vp_cols) %>%
  ggsave(filename = "figs/vp.png", plot=., width=10, height=10, bg="white")

ggplot_gamma2(mj) %>%
  ggsave(filename = "figs/gamma2.png", 
         plot=., width=10, height=5, bg="white")
ggplot_gamma(mj) %>%
  ggsave(filename = "figs/gamma.png", 
         plot=., width=5, height=5, bg="white")
ivars <- c("air_temp_son_pre", "bare", "soil_moisture_pre_jf_30cm",
           "soil_moisture_pre_ma_30cm", "soil_moisture_pre_son_30cm",
           "soil_temp_pre_jf", "soil_temp_pre_ma", "soil_temp_pre_son",
           "total_n_top_15cm_2012", "twi")

lut_ivars <-c("Ta son", "Bare", "Ms jf", "Ms ma", "Ms son", "Ts jf",
                      "Ts ma", "Ts son", "Ns", "TWI")
lut_ivars <-c("Fall AT", "Bare Ground", "Winter SM", "Spring SM", "Fall SM", "Winter ST",
              "Spring ST", "Fall ST", "Soil N", "TWI")

names(lut_ivars) <-  c("air_temp_son_pre", "bare", "soil_moisture_pre_jf_30cm",
               "soil_moisture_pre_ma_30cm", "soil_moisture_pre_son_30cm",
               "soil_temp_pre_jf", "soil_temp_pre_ma", "soil_temp_pre_son",
               "total_n_top_15cm_2012", "twi")
lut_gensp <- c('BROM', 'BASC', 'SATR', 'PASM', 'BOCU', 'MESA', 'BRAS', 
               'NAVI', 'LACT', 'FORB', 'CIAR', 'COAR', 'PAVI', 'SCSC', 
               'ATCA', 'BOGR', 'BRIN', 'SECE')
names(lut_gensp) <- colnames(mj$Y)

source("R/hmsc_plotting_functions.R")
ggplot_beta2_drake(mj, 
                   lut_gensp = lut_gensp,
                   included_variables = ivars,
                   lut_ivars = lut_ivars)%>%
  ggsave(plot = ., filename = "figs/beta2.png", width = 12, height=5, bg = "white")

ggplot_omega(mj,hc.method = "centroid", hc.order = TRUE, 
             lut_gensp=lut_gensp) %>%
  ggsave(filename = "figs/omega.png", plot=., width=7, height=7, bg="white")


