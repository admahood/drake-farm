# figure 1 panel a and b

# figuring out spei

library(tidyverse)
library(SPEI)
library(ggpubr)
library(geomtextpath)
library(vroom)
?SPEI::spei()

# maybe two separate models for strip types

# workaround, just downloaded the global modeled stuff... ======================
spei<- bind_rows(read_delim("data/spei_database/spei06.csv",delim = ";") %>% mutate(scale = "06_month"),
                 read_delim("data/spei_database/spei12.csv",delim = ";") %>% mutate(scale = "12_month"),
                 read_delim("data/spei_database/spei24.csv",delim = ";") %>% mutate(scale = "24_month")) %>%
  dplyr::rename(date = `days since 1900-1-1`)

p1<- spei %>%
  filter(date > as.Date("2010-01-01"),
         scale == "12_month") %>%
  mutate(scale = ifelse(scale == "12_month", "12 Month SPEI", "24 Month SPEI")) %>%
  ggplot(aes(x=date, y=spei)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(2,-2,1, -1), color = "grey", lty=3) +
  geom_labelvline(xintercept = as.Date("2013-04-29"), col="grey40", label = "CRP", hjust=.95, lty=2) +
  geom_labelvline(xintercept = as.Date("2014-05-01"), col="grey40", label = "CRP", hjust=0.05, lty=2) +
  geom_line() +
  xlab("Year") +
  ylab("12 Month SPEI")+
  theme_classic() 


precip <- read_csv('data/monthly_precip.csv') %>%
  pivot_longer(cols = names(.)[2:length(.)],names_to = "year") %>%
  mutate(year = as.numeric(year),
         date = as.Date(paste(year, month, "01", sep="-"), "%Y-%b-%d"),
         nmonth = str_sub(date, 6,7) %>% as.numeric()) %>%
  filter(year>2012, year<2015) %>%
  arrange(date) %>%
  group_by(year) %>%
  mutate(cumsum = cumsum(value)) %>%
  ungroup() 

p_raw<-ggplot(precip, aes(x=as.factor(nmonth), y=value, fill=as.factor(year))) +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic() +
  geom_labelvline(xintercept = 4.5, col="grey40", label = "CRP", hjust=.95, lty=2) +
  scale_color_discrete(name = "Year") +
  xlab("Month") +
  ylab("Precipitation (mm)") +
  scale_fill_manual(values = c("chocolate4", "turquoise3")) +
  theme(legend.position = c(1,1),
        legend.title = element_blank(),
        legend.justification = c(1,1)) 

ggarrange(p1, p_raw, ncol=1, nrow=2, heights = c(1,1)) %>%
  ggsave(plot=., filename = "figs/precip_spei.png", height=4, width=7, bg="white")