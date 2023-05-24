# diversity for the different strips

source("R/drake_data_prep.R")
library(vegan)
library(ggtext)
library(brms)

sp_list <- read_csv("data/drake_veg_data_2022 - species_list.csv") %>%
  dplyr::select(code, introduced)


introduced <- sp_list %>% filter(introduced == 'yes') %>% pull(code)
native <- sp_list %>% filter(introduced == 'no') %>% pull(code)

pc_all <- plant_cover %>%
  dplyr::select(plot, species_code, cover_pct, n_subplots) %>%
  group_by(plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(names_from = "species_code", 
              values_from = "cover_pct",values_fill = 0) %>%
  tibble::column_to_rownames("plot")

# bromus pasm interaction =====================================

pc_all_sub <- raw_veg %>%
  dplyr::select(plot, species_code, cover_pct, subplot) %>%
  mutate(plot = str_c(plot,"__", subplot)) %>%
  pivot_wider(names_from = "species_code", values_from = "cover_pct", id_cols = "plot", values_fill = 0) %>%
  tidyr::separate(plot, c("plot", "subplot"), sep = "__") %>%
  filter(subplot != "999")
glimpse(pc_all_sub)
# models examining pasm brte interaction -- maybe bootstrap, lmer with plot, or do a bayesian thing to get it real nice
# library(lme4)
# m1 <- pc_all_sub %>%
#   mutate(brte_pa = ifelse(brte>0,1,0),
#          ia_brass = sial+deso,
#          n_forbs = stpa + assp+chbe+sotr+upf_12+upfr) %>%
#   dplyr::rename(pasm = elre, scsc = upg4, pavi = sphe, bare = `_bare`) %>%
#   glmer(brte_pa ~ pasm + bocu + bogr + scsc  + ia_brass +
#           kosc + satr + mesa  + pavi + n_forbs+ bare+(1|plot), 
#         data = ., family = "binomial",nAGQ = 3) 
# peff<- ggeffects::ggpredict(m1) %>% 
#   plot(facet=TRUE) + 
#   scale_color_manual(values = c(rep("black",12))) +
#   ylab("P(Bromus tectorum)") +
#   xlab("Percent Cover");peff
# broom.mixed::tidy(m1)
# ggsave(plot=peff,filename = "figs/abundance_vs_brte_pa.png", width =5, height=5, bg="white")

m1b <- pc_all_sub %>%
  mutate(brte_pa = ifelse(brte>0,1,0),
         ia_brass = sial+deso,
         n_forbs = stpa + assp+chbe+sotr+upf_12+upfr) %>%
  dplyr::rename(pasm = elre, scsc = upg4, pavi = sphe, bare = `_bare`) %>%
  brm(brte_pa ~ pasm + bocu + bogr + scsc  + ia_brass +
          kosc + satr + mesa  + pavi + n_forbs+ bare+(1|plot), 
        data = ., family = "bernoulli") 
m2b <- pc_all_sub %>%
  dplyr::rename(pasm = elre, scsc = upg4, pavi = sphe, bare = `_bare`) %>%
  mutate( ia_brass = sial+deso,
          n_forbs = stpa + assp+chbe+sotr+upf_12+upfr) %>%
  mutate_at(c(5:9,11:ncol(.)), function(x)ifelse(x>0,1,0)) %>%
  brm(mvbind(pasm , bocu , bogr , scsc  , ia_brass ,
             kosc , satr , mesa  , pavi , n_forbs) ~ brte + bare+(1|plot), 
      data = ., family = "bernoulli") 

summary(m1b)
p <- conditional_effects(m1b, effects = "pasm", spaghetti = TRUE, ndraws=300) %>%
  plot(line_args=list(colour="black", lwd = 2), theme = theme_classic())
p_pasm <- p$pasm +  
  ylab("P(Cheatgrass)") +
  xlab("Western Wheat Cover (%)") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown())

p2 <- conditional_effects(m2b, effects = "brte", spaghetti = TRUE, ndraws=300) %>%
  plot(line_args=list(colour="black", lwd = 2), theme = theme_classic())
p_brte <- p2$pasm.pasm_brte +  
  ylab("P(Western Wheat)") +
  xlab("Cheatgrass Cover (%)") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown())


ggarrange(p_pasm, p_brte) %>%
  ggsave(plot = ., filename = "figs/pasm_brte.pdf", height =3.5, width=8, bg="white")


ce1 <- conditional_effects(m1b)
p11 <- data.frame(PASM = ce1$pasm$estimate__,  
           BOCU = ce1$bocu$estimate__, 
           BOGR = ce1$bogr$estimate__ , 
           SCSC = ce1$scsc$estimate__ ,
           BRAS = ce1$ia_brass$estimate__,
           BASC = ce1$kosc$estimate__, 
           SATR = ce1$satr$estimate__, 
           MESA  = ce1$mesa$estimate__, 
           PAVI = ce1$pavi$estimate__, 
           FORB = ce1$n_forbs$estimate__,
           Cover =ce1$pasm$effect1__) %>%
  pivot_longer(cols = names(.)[1:ncol(.)-1]) %>%
  ggplot(aes(x=Cover, y=value, color = name)) +
  geom_line()+
  scale_color_discrete(name = "Species Code")+
  ylab("P(Cheatgrass)")+
  xlab("Spp. Cover") +
  theme_classic()

ce2 <- conditional_effects(m2b)
p22 <- data.frame(PASM = ce2$pasm.pasm_brte$estimate__,  
           BOCU = ce2$bocu.bocu_brte$estimate__, 
           BOGR = ce2$bogr.bogr_brte$estimate__, 
           SCSC = ce2$scsc.scsc_brte$estimate__,
           BRAS = ce2$iabrass.iabrass_brte$estimate__,
           BASC = ce2$kosc.kosc_brte$estimate__, 
           SATR = ce2$satr.satr_brte$estimate__, 
           MESA  = ce2$mesa.mesa_brte$estimate__, 
           PAVI = ce2$pavi.pavi_brte$estimate__, 
           FORB = ce2$nforbs.nforbs_brte$estimate__,
           Cover =ce2$pasm.pasm_brte$effect1__) %>%
  pivot_longer(cols = names(.)[1:ncol(.)-1]) %>%
  ggplot(aes(x=Cover, y=value, color = name)) +
  geom_line() +
  xlab("Cheatgrass Cover")+
  ylab("P(Spp. Occurrence)")+
  scale_color_discrete(name = "Species Code") +
  theme_classic() +
  theme(legend.position = "none")

ggarrange(p11, p22, p_pasm, p_brte, nrow=2, ncol=2, common.legend = TRUE) %>%
  ggsave(plot=., filename = "figs/pasm_brte_4pan.png", height=7, width=7, bg="white")

# ggsave(plot = p_pasm, filename = "figs/pasm_brte.png", height=3.5, width=3.5, bg="white")
# conditional_effects(m1b, effects = "bocu", spaghetti = TRUE, ndraws=300)
# conditional_effects(m1b, effects = "ia_brass", spaghetti = TRUE, ndraws=300)
# conditional_effects(m1b, effects = "bare", spaghetti = TRUE, ndraws=300)
# conditional_effects(m1b, effects = "ia_brass", spaghetti = TRUE, ndraws=300)
# conditional_effects(m1b, effects = "ia_brass", spaghetti = TRUE, ndraws=300)

eff<- ggeffects::ggpredict(m1b) 

peff <- eff %>% 
  plot(facet=TRUE) + 
  scale_color_manual(values = c(rep("black",12))) +
  scale_fill_manual(values = c(rep("grey30",12)))+
  ylab("P(*Bromus . tectorum* .  )") +
  xlab("Percent Cover") +
  theme(axis.title.y = element_markdown())
# ggsave(plot=peff,filename = "figs/abundance_vs_brte_pa.png", width =5, height=5, bg="white")

ggarrange(p_pasm, peff, nrow=1, ncol=2, labels = c("b","c")) %>%
ggsave(plot=., filename = "figs/pasm_2pan.png", width =8.5, height=3.5, bg="white")



summary(m2b)

broom.mixed::tidy(m2b) %>%
  filter(effect == "fixed") %>%
  dplyr::select(-component, -group, -effect) %>%
  filter(term != "(Intercept)", term != "bare")

p <- conditional_effects(m2b, effects = "brte",
                         spaghetti = T, ndraws=300) %>%
  plot(line_args=list(colour="black", lwd = 2), theme = theme_classic(), ask=F)
ggarrange(plotlist = p, nrow = 2, ncol = 5, labels = c("d", rep("",9))) %>%
  ggsave(plot = ., filename = "figs/brte_vs_all.png", width = 8, height = 3.5, bg="white")
p_pasm <- p$pasm +  
  ylab("P(*Bromus . tectorum* . )") +
  xlab("*Pascopyrum . smithii* . Cover (%)") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown())

ggeffects::ggpredict(m2b) %>% plot()
summary(m1); summary(m2)
summary(m1b)
broom.mixed::tidy(m1b) %>% 
  filter(term != "(Intercept)", term != "bare")%>%
  filter(effect == "fixed") %>%
  dplyr::select(-component, -group, -effect)
broom.mixed::tidy(m2)

pc_native <- plant_cover %>%
  filter(species_code %in% native) %>%
  dplyr::select(plot, species_code, cover_pct, n_subplots) %>%
  group_by(plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
  unique() %>%
  pivot_wider(names_from = "species_code", 
              values_from = "cover_pct",values_fill = 0) %>%
  tibble::column_to_rownames("plot")

pc_notnative <- plant_cover %>%
  filter(species_code %in% introduced) %>%
  dplyr::select(plot, species_code, cover_pct, n_subplots) %>%
  group_by(plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
  unique() %>%
  pivot_wider(names_from = "species_code", 
              values_from = "cover_pct",values_fill = 0) %>%
  tibble::column_to_rownames("plot")

vegan::diversity(pc_all) %>%
  as_tibble(rownames = "plot") %>%
  tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
  ggplot(aes(x=type, y=value)) +
  geom_boxplot() +
  ggtitle("Total Shannon Diversity")

ggarrange(
vegan::diversity(pc_native) %>%
  as_tibble(rownames = "plot") %>%
  tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
  mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
  ggplot(aes(x=CRP_year, y=value, fill=CRP_year)) +
  geom_violin(draw_quantiles = 0.5, trim=F)+
  theme_classic() +
  guides(fill="none") +
  labs(x = "CRP Application Year", y = "Native Shannon Diversity")
,
vegan::specnumber(pc_native) %>%
  as_tibble(rownames = "plot") %>%
  tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
  mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
  ggplot(aes(x=CRP_year, y=value, fill=CRP_year)) +
  geom_violin(draw_quantiles = 0.5, trim=F)+
  theme_classic() +
  guides(fill="none") +
  labs(x = "CRP Application Year", y = "Native Species Richness")
, 
nrow=1, ncol=2) %>%
  ggsave("figs/diversity_by_strip.png", plot=., width=8.5, height=3.5, bg="white")

vegan::diversity(pc_notnative) %>%
  as_tibble(rownames = "plot") %>%
  tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
  ggplot(aes(x=type, y=value)) +
  geom_boxplot()+
  ggtitle("Introduced Shannon Diversity")


# which species is the most abundant by plot

dom_spp <- plant_cover %>%
  dplyr::select(cover_pct, species_code, plot, n_subplots) %>%
  group_by(plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
  unique()%>%
  ungroup() %>%
  group_by(plot) %>%
  filter(cover_pct == max(cover_pct)) %>%
  dplyr::select(plot, dominant_sp = species_code)

vegan::specnumber(pc_native) %>%
  as_tibble(rownames = "plot") %>%
  tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_",remove = F) %>%
  mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
  left_join(dom_spp) %>%
  mutate(dom_sp = ifelse(dominant_sp %in% c("satr", "elre", "brte", "kosc", "atca", "bocu"),
                         dominant_sp, "other")) %>%
  ggplot(aes(x=dom_sp, y=value)) +
  geom_violin(draw_quantiles = 0.5, trim=F)+
  theme_classic() +
  guides(fill="none") +
  labs(x = "Dominant_sp", y = "Native Shannon Diversity")

vegan::diversity(pc_native) %>%
  as_tibble(rownames = "plot") %>%
  tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_",remove = F) %>%
  mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
  left_join(dom_spp) %>%
  mutate(dom_sp = ifelse(dominant_sp %in% c("satr", "elre", "brte", "kosc", "atca", "bocu"),
                         dominant_sp, "other")) %>%
  ggplot(aes(x=dom_sp, y=value)) +
  geom_violin(draw_quantiles = 0.5, trim=F)+
  theme_classic() +
  guides(fill="none") +
  labs(x = "Dominant_sp", y = "Native Shannon Diversity")
table(dom_spp$dominant_sp)
