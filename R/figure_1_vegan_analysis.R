# setup ========================================================================
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggpubr)
library(Hmsc)

# data ingest ==================================================================
raw <- read_csv("data/drake_veg_data_2022 - cover_2022(1).csv")


lut_spp <- c('COAR', 'SATR', 'BOCU', 'SPHE', 'SCSC', 'BRTE', 'BASC',
                    'NAVI', 'SIAL', 'MESA', 'UPF_12', 'PASM', 'ATCA', 'POLA',
                    'SECE', 'DESO', 'BOGR', 'BRIN', 'ERNA', 'CHBE', 'UPFR', 
                    'HEVI', 'HEAN', 'BRSE', 'STPA', 'MEOF', 'TAOF', 'CIAR', 
                    'ASSP', 'ELEL', 'TRDU', 'LASE', 'TRTE', 'SOTR')
names(lut_spp) <- c('coar', 'satr', 'bocu', 'sphe', 'upg4', 'brte', 'kosc', 'navi',
             'sial', 'mesa', 'upf_12', 'elre', 'atca', 'pola', 'sece', 'deso',
             'bogr', 'brin', 'erna', 'chbe', 'upfr', 'hevi', 'hean', 'brse',
             'stpa', 'meof', 'taof', 'ciar','assp', 'elel', 'trdu', 'lase',
             'trte', 'sotr')
surface_cover <- raw %>%
  filter(str_sub(species_code,1,1)=="_")

plant_cover <- raw %>%
  filter(str_sub(species_code,1,1)!="_") %>%
  mutate(strip_type = str_sub(plot, 7,10)) %>%
  group_by(plot) %>%
  mutate(n_subplots = length(unique(subplot)),
         n_subplots = if_else(n_subplots > 6 , 8, 4)) %>%
  mutate(species_code = lut_spp[species_code])



# nmds abundance-based =========================================================

# creating the community matrix
dv<-plant_cover %>%
  group_by(plot, species_code) %>%
  summarise(cover = (sum(cover_pct)/n_subplots)) %>%
  ungroup() %>%
  unique %>%
  pivot_wider(names_from = species_code, values_from = cover, values_fill = 0) %>%
  tibble::column_to_rownames("plot")

set.seed(1234)
nmds<- dv %>%
  wisconsin() %>%
  vegan::metaMDS(trymax=100)

stressplot(nmds)

site_scores <- as.data.frame(vegan::scores(nmds)$sites) %>%
  as_tibble(rownames = "plot") %>%
  mutate(strip_type = str_sub(plot, 7,10),
         CRP_year = ifelse(strip_type=="shru", "2013", "2014"),
         strip_number = str_sub(plot, 12,13))

ef <- envfit(nmds, dv, na.rm = T, permutations = 9999)
sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.001)

p_ab <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_point(size=2, aes(color = CRP_year)) +
  geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  stat_ellipse(aes(color = CRP_year), key_glyph="rect") +
  theme_classic() +
  theme(panel.background = element_rect(fill="transparent", color = "black"),
        legend.position = "none") +
  ggtitle("Abundance-Based")

# presence/absence-based =======================================================

nmds_pa<- dv %>%
  vegan::decostand(method = "pa") %>%
  vegan::metaMDS(trymax=100)
stressplot(nmds_pa)

site_scores_pa <- as.data.frame(vegan::scores(nmds_pa)$sites) %>%
  as_tibble(rownames = "plot") %>%
  mutate(strip_type = str_sub(plot, 7,10),
         CRP_year = ifelse(strip_type=="shru", "2013", "2014"),
         strip_number = str_sub(plot, 12,13))

ef_pa <- envfit(nmds_pa, dv %>%
                  vegan::decostand(method = "pa"),
                na.rm = T, permutations = 9999)
sp_pa <-as.data.frame(ef_pa$vectors$arrows*sqrt(ef_pa$vectors$r))
species_pa <- as.data.frame(cbind(sp_pa, p=ef_pa$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.001)

p_oc <- ggplot(site_scores_pa, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_point(size=2, aes(color = CRP_year)) +
  geom_segment(data = species_pa,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_text_repel(data = species_pa,size=4, aes(label = species), color = "grey40") +
  stat_ellipse(aes(color = CRP_year), key_glyph="rect") +
  theme_classic() +
  scale_color_discrete(name = "CRP\nYear")+
  theme(panel.background = element_rect(fill="transparent", color = "black"),
        legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_rect(fill=NA)) +
  ggtitle("Occurrence-Based")


ggarrange(p_ab, p_oc, labels = "auto", widths = c(1,1.09)) %>%
  ggsave(filename = "figs/initial_veg_look.png", width = 12, height=6, bg="white")

# native diversity

sp_list <- read_csv("data/drake_veg_data_2022 - species_list.csv") %>%
  dplyr::select(code, introduced)


introduced <- sp_list %>% filter(introduced == 'yes') %>% pull(code)
native <- sp_list %>% filter(introduced == 'no') %>% pull(code)

# pc_all <- plant_cover %>%
#   dplyr::select(plot, species_code, cover_pct, n_subplots) %>%
#   group_by(plot, species_code) %>%
#   summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
#   unique() %>%
#   pivot_wider(names_from = "species_code", 
#               values_from = "cover_pct",values_fill = 0) %>%
#   tibble::column_to_rownames("plot")

pc_native <- raw %>%
  filter(str_sub(species_code,1,1)!="_") %>%
  mutate(strip_type = str_sub(plot, 7,10)) %>%
  group_by(plot) %>%
  mutate(n_subplots = length(unique(subplot)),
         n_subplots = if_else(n_subplots > 6 , 8, 4)) %>%
  filter(species_code %in% native) %>%
  dplyr::select(plot, species_code, cover_pct, n_subplots) %>%
  group_by(plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
  unique() %>%
  pivot_wider(names_from = "species_code", 
              values_from = "cover_pct",values_fill = 0) %>%
  tibble::column_to_rownames("plot")

# pc_notnative <- plant_cover %>%
#   filter(species_code %in% introduced) %>%
#   dplyr::select(plot, species_code, cover_pct, n_subplots) %>%
#   group_by(plot, species_code) %>%
#   summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
#   unique() %>%
#   pivot_wider(names_from = "species_code", 
#               values_from = "cover_pct",values_fill = 0) %>%
#   tibble::column_to_rownames("plot")
# 
# vegan::diversity(pc_all) %>%
#   as_tibble(rownames = "plot") %>%
#   tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
#   ggplot(aes(x=type, y=value)) +
#   geom_boxplot() +
#   ggtitle("Total Shannon Diversity")

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
  nrow=2, ncol=1) %>%
  ggsave("figs/diversity_by_strip.png", plot=., width=8.5, height=3.5, bg="white")

ggarrange(p_ab, p_oc, labels = "auto", widths = c(1,1.09)) %>%
  ggarrange(., ggarrange(
    vegan::diversity(pc_native) %>%
      as_tibble(rownames = "plot") %>%
      tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
      mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
      ggplot(aes(x=CRP_year, y=value, fill=CRP_year)) +
      geom_violin(draw_quantiles = 0.5, trim=F)+
      theme_classic() +
      scale_y_continuous(breaks =c(0,1,2))+
      guides(fill="none") +
      labs(x = "CRP Application Year", y = "Native\nShannon Diversity")
    ,
    vegan::specnumber(pc_native) %>%
      as_tibble(rownames = "plot") %>%
      tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
      mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
      ggplot(aes(x=CRP_year, y=value, fill=CRP_year)) +
      geom_violin(draw_quantiles = 0.5, trim=F)+
      theme_classic() +
      guides(fill="none") +
      labs(x = "CRP Application Year", y = "Native\nSpecies Richness")
    , 
    nrow=2, ncol=1, labels = c("      c", "      d")), nrow=1, ncol=2, widths = c(3,1.25)) %>%
  ggsave("figs/diversity_nmds.png", plot=., width=12, height=4, bg="white")

# # extra analysis not in paper
# vegan::diversity(pc_notnative) %>%
#   as_tibble(rownames = "plot") %>%
#   tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_") %>%
#   ggplot(aes(x=type, y=value)) +
#   geom_boxplot()+
#   ggtitle("Introduced Shannon Diversity")
# 
# 
# # which species is the most abundant by plot
# 
# dom_spp <- plant_cover %>%
#   dplyr::select(cover_pct, species_code, plot, n_subplots) %>%
#   group_by(plot, species_code) %>%
#   summarise(cover_pct = sum(cover_pct)/n_subplots) %>%
#   unique()%>%
#   ungroup() %>%
#   group_by(plot) %>%
#   filter(cover_pct == max(cover_pct)) %>%
#   dplyr::select(plot, dominant_sp = species_code)
# 
# vegan::specnumber(pc_native) %>%
#   as_tibble(rownames = "plot") %>%
#   tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_",remove = F) %>%
#   mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
#   left_join(dom_spp) %>%
#   mutate(dom_sp = ifelse(dominant_sp %in% c("satr", "elre", "brte", "kosc", "atca", "bocu"),
#                          dominant_sp, "other")) %>%
#   ggplot(aes(x=dom_sp, y=value)) +
#   geom_violin(draw_quantiles = 0.5, trim=F)+
#   theme_classic() +
#   guides(fill="none") +
#   labs(x = "Dominant_sp", y = "Native Shannon Diversity")
# 
# vegan::diversity(pc_native) %>%
#   as_tibble(rownames = "plot") %>%
#   tidyr::separate(plot, c("d", "n", "type", "strip"),sep = "_",remove = F) %>%
#   mutate(CRP_year = ifelse(type == "herb", "2014", "2013")) %>%
#   left_join(dom_spp) %>%
#   mutate(dom_sp = ifelse(dominant_sp %in% c("satr", "elre", "brte", "kosc", "atca", "bocu"),
#                          dominant_sp, "other")) %>%
#   ggplot(aes(x=dom_sp, y=value)) +
#   geom_violin(draw_quantiles = 0.5, trim=F)+
#   theme_classic() +
#   guides(fill="none") +
#   labs(x = "Dominant_sp", y = "Native Shannon Diversity")
# table(dom_spp$dominant_sp)
# 
