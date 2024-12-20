# drake strip comparison

source("R/drake_data_prep.R")

long_d <- xdata %>%
  pivot_longer(cols = dplyr::select_if(.,is.numeric) %>% names(),
               names_to = "variable",
               values_to = "value")

p<-long_d %>%
  ggplot(aes(x = strip_type, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free")
  

ggsave(p, filename = "figs/big_facet.png",
       width=20,height=20, bg="white")

all_cover <- plant_cover %>%
  group_by(strip_type, plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct, na.rm=T)/n_subplots) %>%
  ungroup() %>%
  unique()

occurrences <- plant_cover %>%
  group_by(plot, species_code) %>%
  summarise(prevalence = n()) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(id_cols = "plot",
              names_from = "species_code",
              values_from = "prevalence",
              values_fill = 0) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "species_code",
               values_to = "prevalence") %>%
  mutate(strip_type = str_sub(plot, 7,10)) %>%
  group_by(strip_type, species_code) %>%
  summarise(prevalence = sum(prevalence)) %>%
  ungroup() %>%
  group_by(species_code) %>%
  mutate(totalprev = sum(prevalence),
         strip_type = ifelse(strip_type == "herb", 2014, 2013) %>% as.factor) %>%
  ungroup() %>%
  filter(totalprev>0)
  
sp_list <- read_csv("data/drake_veg_data_2022 - species_list.csv") %>%
  dplyr::select(genus, species,species_code = code,group, introduced, perennial, woody, graminoid,
                rhizomatous, pp = photosynthetic_pathway) %>%
  mutate(fg = paste0(ifelse(introduced=="yes", "I", "N"),
                     ifelse(perennial == "yes", "P", "A"),
                     ifelse(graminoid == "yes", "G", "D"))) %>%
  dplyr::select(genus, species, species_code, fg) %>%
  mutate(gensp = paste(genus, species))

p1 <- occurrences %>%
  left_join(sp_list)%>%
  mutate(species_code = ifelse(species_code == "elre", "pasm", species_code),
         species_code = ifelse(species_code == "kosc", "basc", species_code),
         species_code = ifelse(species_code == "upg4", "scsc", species_code),
         species_code = ifelse(species_code == "sphe", "pavi", species_code),
         gensp = ifelse(gensp == "cf Aster d_081_herb_05", "Symphiotricum sp.", gensp),
         gensp = ifelse(gensp == "Rhizomatous_perennial_forb d_141_herb_09", "Unknown rhizomatous\nperennial forb.", gensp),
         gensp = ifelse(gensp == "Perennial_forb d_012_herb_01", "Mirabilis sp.", gensp),
         gensp = ifelse(gensp == "cf Aster d_081_herb_05", "Symphiotricum sp.", gensp),
         gensp = ifelse(gensp == "cf Aster d_081_herb_05", "Symphiotricum sp.", gensp),
         
         species_code = str_to_upper(species_code)) %>%
  filter(str_sub(fg, 1,1) == "N") %>%
  ggplot(aes(x = prevalence, y = gensp, fill=strip_type)) +
  geom_bar(stat = 'identity', position = "dodge") +
  # facet_wrap(~fg, ncol=1, scales = "free") +
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text.y = element_text(face = "italic")) +
  ggtitle("Native Species") +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_fill_manual(values = c("chocolate4", "turquoise3"), 
                    name = "CRP\nYear");p1

p2 <- occurrences %>%
  left_join(sp_list)%>%
  mutate(species_code = ifelse(species_code == "elre", "pasm", species_code),
         species_code = ifelse(species_code == "kosc", "basc", species_code),
         species_code = str_to_upper(species_code)) %>%
  filter(str_sub(fg, 1,1) == "I")%>%
  ggplot(aes(x = prevalence, y = gensp, fill=strip_type)) +
  geom_bar(stat = 'identity', position = "dodge") +
  # facet_wrap(~fg, ncol=1, scales = "free") +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  xlab("Prevalence by subplot") +
  ggtitle("Introduced Species") +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))+
  scale_fill_manual(values = c("chocolate4", "turquoise3"), 
                    name = "CRP\nYear");p2

ggsave(ggarrange(p1, p2, ncol=1, nrow=2, labels = "auto", common.legend = TRUE,
                 legend = "bottom"), filename = "figs/figure_s2_prevalence.png",
       width=10,height=10, bg="white")
