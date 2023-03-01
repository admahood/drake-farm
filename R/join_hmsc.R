# joining hmsc mods together
# devtools::install_github("hmsc-r/HMSC")
library(Hmsc)
library(tidyverse)
library(ggpubr)
library(ggtext)
source("R/hmsc_plotting_functions.R")
mods <- list.files("data/hmsc", full.names=T, pattern = "750")

load(mods[1]); mod1<-m; rm(m)
load(mods[2]); mod2<-m

mj<- c(mod1, mod2)

Hm$XData %>% names() -> nn
nn<-nn[c(8:20)]

ggplot_convergence(mj, gamma=T,omega=T) %>%
  ggsave(plot=.,filename = "figs/convergence.png", width=5, height=5)
library(RColorBrewer)
vp_cols <- c(brewer.pal(12, "Set3"),
             brewer.pal(9, "Set1"),
             "black", "white","burlywood4")
ggplot_vp(mj, cols = vp_cols) %>%
  ggsave(filename = "figs/vp.png", plot=., width=10, height=10, bg="white")

ggplot_gamma(mj)
ggplot_beta2_drake(mj)%>%
  ggsave(plot = ., filename = "figs/beta2.png", width = 18, height=7)
ggplot_beta(mod1,grouping_var = "introduced") 
ggplot_omega(mj) %>%
  ggsave(filename = "figs/omega.png", plot=., width=15, height=15, bg="white")

load(mods[3])
ggplot_ess(m, omega=T)

