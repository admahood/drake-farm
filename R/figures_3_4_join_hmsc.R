# joining hmsc mods together
# devtools::install_github("hmsc-r/HMSC")
library(Hmsc)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggtext)
library(RColorBrewer)
source("R/hmsc_plotting_functions.R")
load(list.files("data/hmsc", full.names=T, pattern = "rr_14"))

mj<-m

# variable definitions ==================
m$XFormula %>% paste() %>% str_replace_all(" \\+", ", ")
ivars_pre <- c("soil_moisture_pre_djf_30cm",  "soil_moisture_pre_ma_30cm",
           "soil_moisture_pre_son_30cm",  "twi", "total_n_top_15cm_2012",
           "soil_temp_pre_djf",  "soil_temp_pre_ma",  "soil_temp_pre_son")
ivars_post <- c("soil_temp_post_jas",  "soil_temp_post_ond",  "soil_temp_post_mj",
                "soil_moisture_post_mj_30cm",  
                "soil_moisture_post_jas_30cm",  "soil_moisture_post_ond_30cm")
ivars_other <- c("twi",  "bare",  "strip_type",  "total_n_top_15cm_2012")

lut_ivars <-c("3. Winter SM (DJF)", "5. Spring SM (MA)", "1. Fall SM (SON)", 
              "1. Spring SM (MJ)", "3. Summer SM (JAS)","5. Fall SM (OND)",
              "4. Summer ST (JAS)", "6. Fall ST (OND)", "2. Spring ST (MJ)", 
              "4. Winter ST (DJF)","6. Spring ST (MA)", "2. Fall ST (SON)",
              "TWI", "Bare Ground", "CRP Year", "Soil N")

names(lut_ivars) <-   c("soil_moisture_pre_djf_30cm",  "soil_moisture_pre_ma_30cm",
                        "soil_moisture_pre_son_30cm",  "soil_moisture_post_mj_30cm",  
                        "soil_moisture_post_jas_30cm",  "soil_moisture_post_ond_30cm",  
                        "soil_temp_post_jas",  "soil_temp_post_ond",  "soil_temp_post_mj",  
                        "soil_temp_pre_djf",  "soil_temp_pre_ma",  "soil_temp_pre_son",  
                        "twi",  "bare",  "strip_type",  "total_n_top_15cm_2012")

lut_tvars <-c("Winter SM (Pre)", "3. Spring SM (Pre)", "Fall SM (Pre)", 
              "Spring SM (Post)", "4. Summer SM (Post)","Fall SM (Post)",
              "5. Summer ST (Post)", "6. Fall ST (Post)", "Spring ST (Post)", 
              "1. Winter ST (Pre)","2. Spring ST (Pre)", "Fall ST (Pre)",
              "7. TWI", "Bare Ground", "CRP Year", "Soil N")

names(lut_tvars) <-   c("soil_moisture_pre_djf_30cm",  "soil_moisture_pre_ma_30cm",
                        "soil_moisture_pre_son_30cm",  "soil_moisture_post_mj_30cm",  
                        "soil_moisture_post_jas_30cm",  "soil_moisture_post_ond_30cm",  
                        "soil_temp_post_jas",  "soil_temp_post_ond",  "soil_temp_post_mj",  
                        "soil_temp_pre_djf",  "soil_temp_pre_ma",  "soil_temp_pre_son",  
                        "twi",  "bare",  "strip_type",  "total_n_top_15cm_2012")

lut_gensp <- c('BROM', 'BASC', 'SATR', 'PASM', 'BOCU', 'MESA', 'BRAS', 
               'NAVI', 'LACT', 'FORB', 'CIAR', 'COAR', 'PAVI', 'SCSC', 
               'ATCA', 'BOGR', 'BRIN', 'SECE')
names(lut_gensp) <- colnames(mj$Y)
lut_gensp_f <- c('Annual Bromes', 'Bassia scoparia', 'Salsola tragus', 
                 'Pascopyrum smithii', 'Bouteloua curtipendula', 'Medicago sativa',
                 'Introduced Brassicaceae', 
               'Nassella viridulis', 'Introduced Chicoriidae', 'Native Forbs',
               'Cirsium arvense', 'Convolvulus arvensis', 'Panicum virgatum', 
               'Schizachrium scoparium', 
               'Atriplex canescens', 'Bouteloua gracilis', 'Bromus inermis', 
               'Secale cereale')
names(lut_gensp_f) <- colnames(mj$Y)
# the plotting ==============================

ggplot_convergence(mj, gamma=T,omega=T) %>%
  ggsave(plot=.,filename = "figs/convergence.png", width=5, height=5)
pf <- ggplot_fit(mj, which = "named", sp_names = lut_gensp_f)
ggsave(plot=pf,filename = "figs/r2.png", width=5, height=5)
ggsave(plot=pf,filename = "figs/r2.pdf", width=5, height=5, dpi=600)

vp_cols <- c(brewer.pal(12, "Set3"),
             brewer.pal(9, "Set1"),
             "black", "white","burlywood4")
ggplot_vp(mj, cols = vp_cols) %>%
  ggsave(filename = "figs/vp.png", plot=., width=10, height=10, bg="white")

ggplot_gamma2(mj) %>%
  ggsave(filename = "figs/gamma2.png", 
         plot=., width=10, height=5, bg="white")

ggplot_gamma(mj, lut_vars = lut_tvars, title = "Trait Associations") %>%
  ggsave(filename = "figs/gamma.png", 
         plot=., width=5, height=5, bg="white")
ggplot_gamma(mj, lut_vars = lut_tvars, title = "Trait Associations") %>%
  ggsave(filename = "figs/gamma.pdf", 
         plot=., width=5, height=5, bg="white", dpi=600)


pb1<- ggplot_beta2_drake(mj, 
                   lut_gensp = lut_gensp,
                   included_variables = ivars_pre,
                   lut_ivars = lut_ivars) +
  ggtitle("Conditions Before Seeding")

ggsave(plot = pb1, filename = "figs/beta2_before.png", width = 12, height=6, bg = "white")
ggsave(plot = pb1, filename = "figs/beta2_before.pdf", width = 12, height=6, bg = "white", dpi=600)

pb2 <- ggplot_beta2_drake(mj, 
                   lut_gensp = lut_gensp,
                   included_variables = ivars_post,
                   lut_ivars = lut_ivars) +
  ggtitle("Conditions After Seeding")
ggsave(plot = pb2, filename = "figs/beta2_after.png", width = 10, height=6, bg = "white")
ggsave(plot = pb2, filename = "figs/beta2_after.pdf",dpi=600, width = 10, height=6, bg = "white")

# pb3 <- ggplot_beta2_drake(mj, 
#                    lut_gensp = lut_gensp_f,
#                    included_variables = ivars_other,
#                    lut_ivars = lut_ivars) +
#   ggtitle("Other Variables")
# ggsave(plot = pb3, filename = "figs/beta2_other.png", width = 9, height=6, bg = "white")
# 
# 
# ggarrange(pb1, pb2, pb3, ncol = 1, nrow = 3, labels = "AUTO") %>%
#   ggsave(plot = ., filename = "figs/beta2_multipanel.png", width = 9, height=18, bg = "white")

p_omega <- ggplot_omega(mj,hc.method = "centroid", hc.order = TRUE, 
             lut_gensp=lut_gensp_f) +
  geom_rect(aes(xmin = 0.5, xmax = 2.5, ymin = 0.5, ymax = 2.5),
            color = "black", fill = "transparent", lwd =1) +
  geom_rect(aes(xmin = 9.5, xmax = 17.5, ymin = 9.5, ymax = 17.5),
            color = "black", fill = "transparent", lwd =1) +
  geom_text(label = "Group 2", x=7, y=14, size=10) +
  geom_text(label = "Group 1", x=3, y=4, size=10);p_omega

ggsave(filename = "figs/omega.png", plot= p_omega, width=9, height=9, bg="white")
ggsave(filename = "figs/omega.pdf", plot= p_omega, width=9, height=9, bg="white", dpi=600)
