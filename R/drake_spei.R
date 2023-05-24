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
  
ggsave(plot= p1 ,filename = "figs/spei_1deg.png",  width=7, height=4.35, bg="white")


# adding in soil moisture ======================================================

es_drake_n <- readxl::read_xlsx(
  "data/past_data/ESDrakeN_2002_2017_daily.xlsx") %>%
  filter(`Date Time` > as.Date("2012-01-01"),
         `Date Time` < as.Date("2016-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999); es_drake_n

es_drake_s <- readxl::read_xlsx(
  "data/past_data/ESDrakeS_2002_2017_daily.xlsx") %>%
  filter(`Date Time` > as.Date("2012-01-01"),
         `Date Time` < as.Date("2016-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999); es_drake_s


soil_moisture_summaries <- bind_rows(es_drake_n, es_drake_s) %>%
  mutate(date = as.Date(`Date Time`),
         month = lubridate::month(date),
         year = lubridate::year(date),
         ym = as.Date(paste(year, month, "01", sep= "-")))
library(ggnewscale)

p2 <- ggplot(soil_moisture_summaries, aes(x=`Date Time`, y=value))+
  geom_line(alpha = 0.25, color="black", aes(group = probe))+
  geom_vline(xintercept = as.Date("2013-05-01") %>% lubridate::as_datetime(), col="red", lty=2) +
  geom_vline(xintercept = as.Date("2014-05-01") %>% lubridate::as_datetime(), col="red", lty=2) +
  geom_line(data = spei %>%
              filter(date > as.Date("2012-01-01"),date < as.Date("2016-01-01")) %>%
              mutate(probe = "NA", depth = "spei"), 
            aes(x=lubridate::as_datetime(date), y = spei, 
                color = scale),
            alpha=0.95, linewidth=1) +
  scale_color_brewer(palette = "Dark2")  +
  facet_wrap(~depth, nrow=6, scales = "free_y")+
  theme_classic()
ggsave(filename = "figs/spei_sentek_compare.png", plot = p2, width=7, height=12,bg="white")

# looking at precip summaries ==================================================

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

p_cum<-ggplot(precip, aes(x=nmonth, y=cumsum, color=as.factor(year))) +
  geom_line() +
  geom_vline(xintercept = 5) +
  theme_classic()

p_raw<-ggplot(precip, aes(x=as.factor(nmonth), y=value, fill=as.factor(year))) +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic() +
  geom_labelvline(xintercept = 4.5, col="grey40", label = "CRP", hjust=.95, lty=2) +
  scale_color_discrete(name = "Year") +
  xlab("Month") +
  ylab("Precipitation (mm)") +
  theme(legend.position = c(1,1),
        legend.title = element_blank(),
        legend.justification = c(1,1)) 

ggarrange(p_cum, p_raw,ncol = 1, nrow=2) %>%
  ggsave(plot=., filename = "figs/precip_summaries.png", height=4, width=8, bg="white")

ggarrange(p1, p_raw, ncol=1, nrow=2, heights = c(1,1)) %>%
  ggsave(plot=., filename = "figs/precip_spei.png", height=4, width=7, bg="white")

# getting spei from drake data =================================================
drake_precip_df <- read_csv("/home/a/projects/drake-farm/data/past_data/ages_input/drake58hru/data/reg_precip.csv",
                            skip = 21) %>%
  dplyr::select(2:60) %>%
  pivot_longer(cols = names(.)[2:length(.)],
               names_to = "sensor", values_to = "precip") %>%
  mutate(date = as.Date(time, "%d.%m.%Y"),
         ym = str_sub(time, 4,10)) %>%
  dplyr::select(-time) %>%
  group_by(date, ym) %>%
  summarise(mean_precip = mean(precip)) %>%
  ungroup() %>%
  group_by(ym) %>%
  summarise(sum_precip = sum(mean_precip)) %>%
  ungroup() %>%
  mutate(year = str_sub(ym, 4, 7) %>% as.numeric(),
         month = str_sub(ym, 1,2) %>% as.numeric()) %>%
  arrange(year,month)

drake_tmean_df <- read_csv("/home/a/projects/drake-farm/data/past_data/ages_input/drake58hru/data/reg_tmean.csv",
                        skip=13) %>%
  dplyr::select(2:60) %>%
  pivot_longer(cols = names(.)[2:length(.)],
               names_to = "sensor", values_to = "temp") %>%
  mutate(date = as.Date(time, "%d.%m.%Y"),
         ym = str_sub(time, 4,10)) %>%
  group_by(ym) %>%
  summarise(mean_temp = mean(temp)) %>%
  ungroup() %>%
  mutate(year = str_sub(ym, 4, 7) %>% as.numeric(),
         month = str_sub(ym, 1,2) %>% as.numeric()) %>%
  arrange(year,month); drake_tmean_df

drake_tmean <- drake_tmean_df %>%
  pull(mean_temp);drake_tmean

pet <- SPEI::thornthwaite(Tave = drake_tmean,lat = 40.605054) %>%
  as.vector()
precip <- pull(drake_precip_df, sum_precip)

balance <- precip - pet

spei06 <- SPEI::spei(balance, scale = 6)
spei12 <- SPEI::spei(balance, scale = 12)
spei24 <- SPEI::spei(balance, scale = 24)

fulldf <- drake_tmean_df %>%
  mutate(date = as.Date(paste0(ym, ".01"), "%m.%Y.%d"),
         precip = precip,
         pet = pet,
         spei06 = spei06$fitted %>% as.vector(),
         spei12 = spei12$fitted %>% as.vector(),
         spei24 = spei24$fitted %>% as.vector())

p_agis_spei <- ggplot(fulldf %>% filter(year>2011, year < 2016) %>%
                        pivot_longer(cols = c("spei06", "spei12", "spei24"),
                                     names_to = "scale", values_to = "spei"), 
       aes(x=date, y = spei, color = scale)) +
  geom_line() +
  geom_textvline(xintercept = c(as.Date("2013-05-15")), lty=2, label = "CRP seeding",hjust=.1 ) +
  geom_textvline(xintercept = c(as.Date("2014-05-15")), lty=2, label = "CRP seeding", hjust=.1) +
  ylab("SPEI") +
  theme_light()+
  ggtitle("SPEI from drake AGIS inputs");p_agis_spei

ggsave(p_agis_spei, filename="figs/spei_from_agis_input.png", width=7.5, height =5)  

p_regression <- ggplot(fulldf %>% filter(year>2011, year < 2016) %>%
                         pivot_longer(cols = c("spei06", "spei12", "spei24"),
                                      names_to = "scale", values_to = "spei"),
                       aes(x=spei, y=mean_temp, color = scale)) +
  geom_point()


## temp in situ ==========

library(sf)
library(ggthemes)
cus <- st_read("data/CUS/")
drake <- data.frame(y=40.6, x= -104.84) %>%
  st_as_sf(coords = c("x","y"),crs = 4326)
dmap<- ggplot(cus) + 
  geom_sf(fill="white") +
  geom_sf(data=drake, color="red") +
  theme_classic() +
  theme(panel.background = element_rect(fill=NA, color="black", linewidth=1))

ggsave("figs/drake_map.png", height = 2.5, width=4.5, plot = dmap, bg="white")

# map stuff
annual_precip<- c(183,	289,	285,	334,	181,	238,	293,	456,	418,	396,
                  247,	394,	386,	462,	276,	388,	275,	316,	210,	232)
map <- mean(annual_precip)
sdap <- sd(annual_precip)
