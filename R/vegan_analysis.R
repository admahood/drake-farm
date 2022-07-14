# setup ========================================================================
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggpubr)
library(Hmsc)

# data ingest ==================================================================
raw <- read_csv("data/drake_veg_data_2022 - cover_2022(1).csv")


surface_cover <- raw %>%
  filter(str_sub(species_code,1,1)=="_")

plant_cover <- raw %>%
  filter(str_sub(species_code,1,1)!="_") %>%
  mutate(strip_type = str_sub(plot, 7,10)) %>%
  group_by(plot) %>%
  mutate(n_subplots = length(unique(subplot)),
         n_subplots = if_else(n_subplots > 6 , 8, 4))

# unique(plant_cover$species_code) %>% sort



# looking at heights ===========================================================
ggplot(plant_cover, aes(x=species_code, y=height_cm, fill=strip_type)) +
  geom_boxplot() +
  facet_wrap(~species_code, scales = "free")

# nmds =========================================================================

dv<-plant_cover %>%
  group_by(plot, species_code) %>%
  summarise(cover = (sum(cover_pct)/n_subplots)) %>%
  ungroup() %>%
  unique %>%
  pivot_wider(names_from = species_code, values_from = cover, values_fill = 0) %>%
  tibble::column_to_rownames("plot")

# abundance-based

nmds<- dv %>%
  wisconsin() %>%
  vegan::metaMDS(trymax=10000)


site_scores <- as.data.frame(vegan::scores(nmds)) %>%
  as_tibble(rownames = "plot") %>%
  mutate(strip_type = str_sub(plot, 7,10),
         strip_number = str_sub(plot, 12,13))

ef <- envfit(nmds, dv, na.rm = T, permutations = 9999)
sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.01)

p_ab <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  geom_text(size=4, aes(label = strip_number, color = strip_type)) +
  stat_ellipse(aes(color = strip_type)) +
  theme_classic() +
  theme(panel.background = element_rect(fill="transparent", color = "black")) +
  ggtitle("abundance-based")

# presence/absence-based


nmds_pa<- dv %>%
  vegan::decostand(method = "pa") %>%
  vegan::metaMDS(trymax=10000)

# sp_scores_pa <- as.data.frame(scores(nmds_pa)$species)%>%
#   as_tibble(rownames = "species") 

site_scores_pa <- as.data.frame(vegan::scores(nmds_pa)) %>%
  as_tibble(rownames = "plot") %>%
  mutate(strip_type = str_sub(plot, 7,10),
         strip_number = str_sub(plot, 12,13))

ef_pa <- envfit(nmds_pa, dv %>%
                  vegan::decostand(method = "pa"),
                na.rm = T, permutations = 9999)
sp_pa <-as.data.frame(ef_pa$vectors$arrows*sqrt(ef_pa$vectors$r))
species_pa <- as.data.frame(cbind(sp_pa, p=ef_pa$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.01)

p_oc <- ggplot(site_scores_pa, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_segment(data = species_pa,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_text_repel(data = species_pa,size=4, aes(label = species), color = "grey40") +
  geom_text(size=4, aes(label = strip_number, color = strip_type)) +
  stat_ellipse(aes(color = strip_type)) +
  theme_classic() +
  theme(panel.background = element_rect(fill="transparent", color = "black")) +
  ggtitle("occurrence-based")


ggarrange(p_ab, p_oc, common.legend = TRUE,legend = "bottom") %>%
  ggsave(filename = "figs/initial_veg_look.png", width = 12, height=6, bg="white")
