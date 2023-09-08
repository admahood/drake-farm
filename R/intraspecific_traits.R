# intraspecific traits

library(tidyverse)
library(janitor)

dh <- readxl::read_xlsx("data/drake_veg_data.xlsx",sheet = "height")
dl <- readxl::read_xlsx("data/drake_veg_data.xlsx",sheet = "LDMC") %>%
  filter(species != "brte")
ds <- readxl::read_xlsx("data/drake_veg_data.xlsx",sheet = "seed_mass") %>%
  mutate(plot = ifelse(plot == "-", NA, plot))

ggplot(dl) +
  geom_density(aes(x=LDMC, fill=species), alpha=0.5)

ggplot(ds) +
  geom_density(aes(x=log(seed_mass_mg), fill=species), alpha=0.5) 

ggplot(dh) +
  geom_density(aes(x=height, fill=species), alpha=0.5) +
  facet_wrap(~species, scales = "free")

dl %>%
  filter(species == "navi") %>%
  arrange(desc(LDMC))




# plot-level variation?

da <- read_csv("data/past_data/Drake_Carbon_2001_2012.csv") %>%
  janitor::clean_names() %>% 
  dplyr::select(plot = x2012smp_no, mgt = mgt_strip, cosasp, elev_2012, 
                sc = x12c_0_6in,  ssoc = x12soc_0_6, sn = x12n_0_6in,
                slope, potsolrad, curv, lnscasink, wetsink, lnscafill, wetfill)

dxx <- left_join(dl %>%
                   mutate(species = str_to_lower(str_sub(species, 1,2))), dh %>%
                   mutate(species = str_to_lower(str_sub(species, 1,2))),
                 by =c("plot", "species", "individual")) %>%
  left_join(da)

ggplot(dxx, aes(x=height, y = LDMC, color = species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~species)

library(lmerTest)
lmer(LDMC ~ height*species + (1|plot/species/individual), data = dxx) %>% summary


dl <- dl %>%
  left_join(da)

ggplot(dl) +
  geom_point(aes(x=LDMC, y=as.factor(plot), color = mgt)) +
  facet_wrap(~species, scales = "free")

ggplot(dl) +
  geom_boxplot(aes(x=species, y=LDMC, fill =mgt))

lmer(LDMC ~ mgt + sc + sn +ssoc + cosasp+  (1|plot),
     data = dl %>% filter(species == "pasm")) -> pamod; pamod %>% summary

lmer(LDMC ~ mgt +sc+ sn+ ssoc+(1|plot),
     data = dl %>% filter(species == "basc")) %>% summary

lmer(LDMC ~ mgt + sc + sn +ssoc + cosasp+ (1|plot),
     data = dl %>% filter(species == "bocu")) ->bamod;bamod %>% summary
lmer(LDMC ~ mgt + sc + sn +ssoc + cosasp+ (1|plot),
     data = dl %>% filter(species == "satr")) %>% summary
lmer(LDMC ~ sc + sn +ssoc + cosasp+ (1|plot),
     data = dl %>% filter(species == "navi")) %>% summary

ggeffects::ggpredict(bamod, terms = c("sc")) -> partialeff

plot(partialeff)

ggeffects::ggpredict(bamod) -> partialeff

plot(partialeff$ssoc)
plot(partialeff$sn)

