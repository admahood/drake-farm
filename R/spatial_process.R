# spatial models

elev <- elv_df[,c(4:8)]




# spacetime model
smd <- bind_rows(es_drake_n, es_drake_s) %>%
  # filter(probe != "C1") %>%
  mutate(year = lubridate::year(`Date Time`),
         month = lubridate::month(`Date Time`),
         month_group = case_when(month < 3 ~ "pre_jf",
                                 month < 5 & month > 3 ~ "pre_ma",
                                 month > 5 & month < 9 ~ "post_jja",
                                 month > 8 & month <12 ~ "post_son"))%>%
  group_by(year, month, probe, depth)%>%
  summarise(mean_soil_moisture = median(value, na.rm=TRUE),
            n=n()) %>%
  ungroup() %>%
  left_join(sentek_locations, by = c("probe" = "Probe")) %>%
  dplyr::rename(elevation = groundEL_1) %>%
  filter(depth == dpth)


nrow(smd) == length(unique(smd$probe)) * 2 * 12

gpmod <- spT.Gibbs(formula = mean_soil_moisture ~ elevation + fa + twi + TPI + slope,
                   model = 'GP', coords = ~UTME+UTMN,data = smd, 
                   time.data = spT.time(t.series = c(12),2))

predict(object = gpmod, type = "temporal",
        newcoords = dplyr::select(elv_df, UTME=x, UTMN=y) %>% as.matrix(), 
        newdata = elv_df[4:8] %>% dplyr::rename(elevation = dem2018))


# krig
ggplot(smd, aes(x=fa, y=mean_soil_moisture, color = as.factor(month))) +
  geom_point()
smdm <- smd %>% filter(month == 1, year == 2013) %>%
  mutate(mean_soil_moisture = mean_soil_moisture/100)
library(gstat)
vsm <- gstat::variogram(smdm$mean_soil_moisture ~ 1, locations = ~UTMN+UTME, data=smdm%>% as.data.frame())
plot(vsm)

best_vario <- fit.variogram(vsm,  model = vgm(c( "Sph", "Exp")))
best_vario

dists <- smdm %>%
  dplyr::select(x=UTME,y=UTMN) %>% 
  dist() %>%
  as.matrix()

#inverse of distance is weight
dists_weights <- 1/dists
diag(dists_weights) <- 0
library(spdep)
#turn into a weights list
W <- mat2listw(dists_weights)

library(brms)
bor_brms  <- brm(mean_soil_moisture ~ elevation + fa + TPI + slope,
                 family = "beta",
                 autocor = cor_car(W = dists),
                 data = smdm)

library(glmmTMB)

#get things prepped
smdm <- smdm %>% 
  mutate(position = numFactor(UTME, UTMN),
         group = factor(1))

bor_tmb  <- glmmTMB(mean_soil_moisture ~ elevation + fa + TPI + slope+ exp(position + 0 | group),
                    family = beta_family(),
                    data = smdm)

