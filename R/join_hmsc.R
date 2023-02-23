# joining hmsc mods together
# devtools::install_github("hmsc-r/HMSC")
library(Hmsc)
library(tidyverse)
library(ggpubr)
library(ggtext)
source("R/hmsc_plotting_functions.R")
mods <- list.files("data/hmsc", full.names=T, pattern = "oblas")

load(mods[1]); mod1<-m; rm(m)
load(mods[2]); mod2<-m

mj<- c(mod1, mod2)

ggplot_convergence(mj, omega=T)
ggplot_gamma(mj)
ggplot_beta(mj,grouping_var = "introduced")
ggplot_omega(mj)
