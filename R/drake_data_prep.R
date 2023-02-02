# data prep for drake
library(tidyverse)
library(sf)
library(fields)
library(terra)
library(raster)
library(topmodel) # for the twi layers
viz<-FALSE
set.seed(331)
# functions ====================================================================

# https://stackoverflow.com/questions/58553966/calculating-twi-in-r
upslope <- function (dem, log = TRUE, atb = FALSE, deg = 0.12, fill.sinks = TRUE) 
{
  if (!all.equal(xres(dem), yres(dem))) {
    stop("Raster has differing x and y cell resolutions. Check that it is in a projected coordinate system (e.g. UTM) and use raster::projectRaster to reproject to one if not. Otherwise consider using raster::resample")
  }
  if (fill.sinks) {
    capture.output(dem <- invisible(raster::setValues(dem, 
                                                      topmodel::sinkfill(raster::as.matrix(dem), 
                                                                         res = xres(dem), 
                                                                         degree = deg))))
  }
  topidx <- topmodel::topidx(raster::as.matrix(dem), res = xres(dem))
  a <- raster::setValues(dem, topidx$area)
  if (log) {
    a <- log(a)
  }
  if (atb) {
    atb <- raster::setValues(dem, topidx$atb)
    a <- addLayer(a, atb)
    names(a) <- c("a", "atb")
  }
  return(a)
}

create_layers <- function (dem, fill.sinks = TRUE, deg = 0.1) 
{
  layers <- stack(dem)
  message("Building upslope areas...")
  a.atb <- upslope(dem, atb = TRUE, fill.sinks = fill.sinks, deg = deg)
  layers <- addLayer(layers, a.atb)
  names(layers) <- c("filled.elevations", "upslope.area", "twi")
  return(layers)
}

# veg data =====================================================================

raw_veg <- read_csv("data/drake_veg_data_2022 - cover_2022(1).csv")

surface_cover <- raw_veg %>%
  filter(str_sub(species_code,1,1)=="_")

plant_cover <- raw_veg %>%
  filter(str_sub(species_code,1,1)!="_") %>%
  mutate(strip_type = str_sub(plot, 7,10)) %>%
  group_by(plot) %>%
  mutate(n_subplots = length(unique(subplot)),
         n_subplots = if_else(n_subplots > 6 , 8, 4))

# topography data ==============================================================

elevation <- terra::rast("data/dem2018.tif")
elevation_r <- raster::raster("data/dem2018.tif", value=TRUE)

twi_stack <- create_layers(elevation_r)

# soil C stuff =================================================================

plots_w_veg <- raw_veg %>%
  transmute(plot = str_sub(plot,3,5)) %>%
  pull(plot) %>%
  unique %>%
  as.numeric()

pipeline_installation <-  c(17, 34, 51, 52, 68, 85, 86, 102, 103, 120, 136, 137,
                            154, 171, 189)

soil_c <- read_csv("data/past_data/Drake_Carbon_2001_2012.csv") %>%
  dplyr::rename(plot = "2012Smp_No",
                strip_type = Mgt_Strip,
                strip_number = Strip_No,
                carbonates_top_15cm_2012 = "12CC_0_6in",
                total_c_top_15cm_2012 = "12C_0_6in",
                total_n_top_15cm_2012 = "12N_0_6in",
                organic_c_top_15cm_2012 = "12SOC_0_6",
                carbonates_15_30cm_2012 = "12CC_6_12",
                total_c_15_30cm_2012 = "12C_6_12in",
                total_n_15_30cm_2012 = "12N_6_12in",
                organic_c_15_30cm_2012 = "12SOC_6_12",
                carbonates_top_30cm_2012 = "12CC_0_1ft",
                total_c_top_30cm_2012 = "12C_0_1ft",
                total_n_top_30cm_2012 = "12N_0_1ft",
                organic_c_top_30cm_2012 = "12SOC_0_1",
                carbonates_top_30cm_2001 = "01CC_0_1ft",
                carbonates_30_60cm_2001 = "01CC_1_2ft",
                carbonates_60_90cm_2001 = "01CC_2_3ft",
                carbonates_top_90cm_2001 = "01CC_0_3ft",
                carbonate_difference = CC_diff,
                carbonate_n3_difference = CC_n3diff,
                utm_e = UTME_2012,
                utm_n = UTMN_2012,
                elv_2001 = ELEV_2001,
                elv_2012 = ELEV_2012,
                elf_difference = ELEVDIFF,
                chg_class_5 = ChgClass_5,
                slope = SLOPE,
                aspect = ASPECT,
                cosign_aspect = COSASP,
                potsolrad = POTSOLRAD,
                curvature = CURV,
                lnscasink = "LNSCASINK",
                wetsink = "WETSINK",
                lnscafill = "LNSCAFILL",
                wetfill = "WETFILL",
                soil_unit = "SOILUNIT",
                soil_unit_name = "UNITNAME") %>%
  mutate(strip_type = ifelse(strip_type == "East", "herb", "shru"),
         seeding_year = ifelse(strip_type == "herb", 2014, 2013),
         pipeline = ifelse(plot %in% pipeline_installation, 
                           "pipeline", "no_pipeline")) %>%
  filter(plot %in% plots_w_veg)
glimpse(soil_c)


# sentek soil temps ============================================================
sentek <- readxl::read_xlsx(
  "data/past_data/sentek probe temperature 2002-2017.xlsx",
  skip=1)

soil_temp_summary <- sentek %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "probe_depth",
               values_to = "soil_temp_c") %>%
  mutate(probe = str_sub(probe_depth, 1,2),
         depth = str_sub(probe_depth, 4,7)) %>%
  filter(depth == "30cm",
         `Date/Time` > as.Date("2012-01-01"),
         `Date/Time` < as.Date("2014-12-01")) %>%
  mutate(year = lubridate::year(`Date/Time`),
         month = lubridate::month(`Date/Time`)) 

t_antecedent_13 <-  soil_temp_summary %>%
  filter(year == 2012) %>%
  mutate(month_group = case_when(month > 8 & month < 12 ~ "pre_seed_son"),
         year = 2013)%>%
  group_by(year, month_group, probe) %>%
  summarise(mean_30cm_soil_temp_c = mean(soil_temp_c, na.rm=TRUE),
            min_30cm_soil_temp_c = min(soil_temp_c, na.rm=TRUE)) %>%
  ungroup()
  
soil_temp_13 <- soil_temp_summary %>%
  filter(year == 2013) %>%
  mutate(month_group = case_when(month < 3 ~ "pre_seed_jf",
                                 month < 5 & month > 3 ~ "pre_seed_ma",
                                 month >5 & month <= 8 ~ "post_seed_jja",
                                 month > 8 & month < 12 ~ "post_seed_son"))%>%
  group_by(year, month_group, probe) %>%
  summarise(mean_30cm_soil_temp_c = mean(soil_temp_c, na.rm=TRUE),
            min_30cm_soil_temp_c = min(soil_temp_c, na.rm=TRUE)) %>%
  ungroup() %>%
  bind_rows(t_antecedent_13)
summary(soil_temp_13); glimpse(soil_temp_13);soil_temp_13

t_antecedent_14 <-  soil_temp_summary %>%
  filter(year == 2013) %>%
  mutate(month_group = case_when(month > 8 & month < 12 ~ "pre_seed_son"),
         year = 2014)%>%
  group_by(year, month_group, probe) %>%
  summarise(mean_30cm_soil_temp_c = mean(soil_temp_c, na.rm=TRUE),
            min_30cm_soil_temp_c = min(soil_temp_c, na.rm=TRUE)) %>%
  ungroup()
soil_temp_14 <- soil_temp_summary %>%
  filter(year == 2014) %>%
  mutate(month_group = case_when(month < 3 ~ "pre_seed_jf",
                                 month < 5 & month > 3 ~ "pre_seed_ma",
                                 month > 5 & month <= 8 ~ "post_seed_jja",
                                 month > 8 & month < 12 ~ "post_seed_son"))%>%
  group_by(year, month_group, probe) %>%
  summarise(mean_30cm_soil_temp_c = mean(soil_temp_c, na.rm=TRUE),
            min_30cm_soil_temp_c = min(soil_temp_c, na.rm=TRUE)) %>%
  ungroup() %>%
  bind_rows(t_antecedent_14)
summary(soil_temp_14); glimpse(soil_temp_14);soil_temp_14


terrain0 <- terra::terrain(elevation, v = c("aspect", "slope", "TPI"))
terrain <- c(terrain0, rast(twi_stack$twi))

sentek_locations <- st_read(
  "data/past_data/Sentek probe locations/sentek_probes.shp") %>%
  bind_cols(terra::extract(terrain, vect(.))) %>%
  mutate(fa = abs(180 - abs(aspect - 225)))

sentek_13 <- sentek_locations %>%
  left_join(soil_temp_13 %>%
              dplyr::select(probe, mean_30cm_soil_temp_c, month_group),
            by=c("Probe"="probe"))%>%
  mutate(year = 2013)

sentek_14<- sentek_locations %>%
  left_join(soil_temp_14 %>%
              dplyr::select(probe, mean_30cm_soil_temp_c, month_group),
            by=c("Probe"="probe"))%>%
  mutate(year = 2014)

sentek_summaries <- bind_rows(sentek_13, sentek_14)

# consider 3dep, or other 1m dem options

elv_df<-c(terrain, elevation)%>%
  as.data.frame(xy=TRUE)%>%
  mutate(fa = abs(180 - abs(aspect - 225))) %>%
  na.omit()
lat <- elv_df$y
lon <- elv_df$x
elev <- elv_df[,4:8]
xps <-cbind(lon, lat)

counter <- 1
result <- list()
spmods <- list()
library(ggpubr)
for(month_grp in (unique(sentek_13$month_group) %>% na.omit())){
  for(sample_year in c(2013:2014)){
    
    # data prep for spatial process model
    x = st_coordinates(sentek_locations)[,c(1,2)]
    z = dplyr::select(sentek_locations, slope, TPI, twi, groundEL_1, fa) %>%
      st_set_geometry(NULL)
    y = sentek_summaries %>%
      filter(year == sample_year, month_group == month_grp) %>%
      pull(mean_30cm_soil_temp_c)
    
    print(paste(month_grp, sample_year))
    spmods[[counter]] <- fields::spatialProcess(x=x, y=y, Z=z, 
                                    profileLambda = TRUE, 
                                    profileARange = TRUE)
    
    # print(summary(spmods[[counter]]))
    ## Predict using the spatial process model fitted
    yp = predict(spmods[[counter]],
                 x=xps,Z=elev)
    
    # convert to a raster
    spmod_rast <- raster::rasterFromXYZ(data.frame(lon=lon, lat=lat,
                                                   prediction = yp),
                                        crs = raster::crs(elevation)) %>%
      rast()
    names(spmod_rast) <- paste0(month_grp, "_", sample_year)
    result[[counter]] <- spmod_rast

    p<- ggarrange(
      spmod_rast %>%
      as.data.frame(xy=TRUE) %>%
        pivot_longer(cols = names(.)[3:length(.)], 
                     values_to = "temp_c", 
                     names_to = "month_grp") %>%
        ggplot() +
        geom_raster(aes(x=x, y=y, fill=temp_c)) + 
        scale_fill_viridis(option="B") +
        facet_wrap(~month_grp,ncol = 1)+
        geom_sf(data = sentek_locations %>%
                  mutate(sentek_t = y), aes(color=sentek_t)) +
        coord_sf(),
    ggarrange(
      ggplot(as.data.frame(y), aes(x=y)) +
        geom_histogram()+
        ggtitle("ST Input", month_grp),
      ggplot(as_tibble(spmod_rast) %>% dplyr::select(x=1), aes(x=x)) +
        geom_histogram()+
        ggtitle("ST predictions",month_grp),
    nrow=2, ncol=1)
    ,nrow=1, ncol=2, widths = c(4,1)
    )
    ggsave(filename = paste0("figs/surfaces/",
                             paste0(month_grp, "_", 
                                    sample_year, "_st.png")),
           plot= p, width=10, height=5, bg="white")
    counter <- counter +1
}}

# model diagnostics

# for(i in 1:length(spmods)){
#   set.panel(2,2);plot(spmods[[i]]);set.panel(1,1)
# }

# sentek soil moisture ================================================================
a_es_drake_e <- readxl::read_xlsx(
  "data/past_data/E.wc.hr.2002-18_Nezat.xlsx",col_types=c("date",rep("numeric",16))) %>%
  filter(timestamp > as.Date("2012-01-01"),
         timestamp < as.Date("2014-01-01")) %>%
  dplyr::select(-numericdate) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,1,2),
         depth = str_sub(variable,4,6) %>% str_c("cm")) %>%
  dplyr::select(-variable) %>%
  filter(value > -9999); a_es_drake_e

# ggplot(a_es_drake_e, aes(x=value)) +
#   geom_histogram()

a_es_drake_n <- readxl::read_xlsx(
  "data/past_data/ESDrakeN_2002_2017_daily.xlsx") %>%
  filter(`Date Time` > as.Date("2012-01-01"),
         `Date Time` < as.Date("2014-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999)%>%
  dplyr::select(-variable, -measure) %>%
  dplyr::rename(timestamp = 1); a_es_drake_n

a_es_drake_s <- readxl::read_xlsx(
  "data/past_data/ESDrakeS_2002_2017_daily.xlsx") %>%
  filter(`Date Time` > as.Date("2012-01-01"),
         `Date Time` < as.Date("2014-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999) %>%
  dplyr::select(-variable, -measure) %>%
  dplyr::rename(timestamp = 1); a_es_drake_s


antecedent_soil_moisture_summaries <- bind_rows(a_es_drake_n, 
                                                a_es_drake_s,
                                                a_es_drake_e) %>%
  mutate(year = lubridate::year(timestamp),
         month = lubridate::month(timestamp),
         month_group = case_when(month > 8 & month < 12 ~ "pre_son"),
         year = year +1)%>%
  group_by(year, month_group, probe, depth)%>%
  summarise(mean_soil_moisture = median(value, na.rm=TRUE)) %>%
  ungroup()

es_drake_e <- readxl::read_xlsx(
  "data/past_data/E.wc.hr.2002-18_Nezat.xlsx",
  col_types=c("date",rep("numeric",16))) %>%
  filter(timestamp > as.Date("2013-01-01"),
         timestamp < as.Date("2015-01-01"))  %>%
  dplyr::select(-numericdate) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,1,2),
         depth = str_sub(variable,4,6) %>% str_c("cm")) %>%
  dplyr::select(-variable) %>%
  filter(value > -9999); es_drake_e

es_drake_n <- readxl::read_xlsx(
  "data/past_data/ESDrakeN_2002_2017_daily.xlsx") %>%
  dplyr::rename(timestamp=1) %>%
  filter(timestamp > as.Date("2013-01-01"),
         timestamp < as.Date("2015-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999) %>%
  dplyr::select(-variable, -measure); es_drake_n

es_drake_s <- readxl::read_xlsx(
  "data/past_data/ESDrakeS_2002_2017_daily.xlsx") %>%
  dplyr::rename(timestamp=1) %>%
  filter(timestamp > as.Date("2013-01-01"),
         timestamp < as.Date("2015-01-01")) %>%
  pivot_longer(cols = names(.)[2:ncol(.)], 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(probe = str_sub(variable,2,3),
         measure = str_split(variable, "\\]", simplify = T)[,2],
         depth = str_split(str_split(variable, "\\]",
                                     simplify = T)[,1],
                           "_", 2, simplify = T)[,2]) %>%
  filter(measure != "Raw", value > -9999) %>%
  dplyr::select(-variable, -measure); es_drake_s


soil_moisture_summaries <- bind_rows(es_drake_n, 
                                     es_drake_s,
                                     es_drake_e) %>%
  mutate(year = lubridate::year(timestamp),
         month = lubridate::month(timestamp),
         month_group = case_when(month < 3 ~ "pre_jf",
                                 month < 5 & month > 3 ~ "pre_ma",
                                 month > 5 & month < 9 ~ "post_jja",
                                 month > 8 & month <12 ~ "post_son"))%>%
  group_by(year, month_group, probe, depth)%>%
    summarise(mean_soil_moisture = median(value, na.rm=TRUE)) %>%
    ungroup() %>%
  bind_rows(antecedent_soil_moisture_summaries)

dpth <- c("30cm")

sm_sentek <- sentek_locations %>%
  left_join(soil_moisture_summaries, by = c("Probe" = "probe")) %>%
  na.omit() %>%
  filter(depth == dpth)

# ggplot(sentek_locations) +
#   geom_text(aes(label=Probe,x=UTME,y=UTMN))

sm_sentek %>%
  group_by(Probe) %>%
  summarise(n=n())

counter <- 1
result_sm <- list()
spmods_sm <- list()
for(month_grp in (unique(sm_sentek$month_group))){
  for(sample_year in c(2013:2014)){
    for(dpt in dpth){
    d <- filter(sm_sentek, 
                year == sample_year, 
                month_grp == month_group,
                depth == dpt)
    # data prep for spatial process model
    x = st_coordinates(d)[,c(1,2)]
    z = dplyr::select(d,  slope, TPI, twi, groundEL_1, fa) %>%
      st_set_geometry(NULL)
    y = sm_sentek %>%
      filter(year == sample_year, 
             month_group == month_grp,
             depth == dpt) %>%
      pull(mean_soil_moisture)
    
    print(paste(month_grp, sample_year,dpt))
    spmods_sm[[counter]] <- fields::spatialProcess(x=x, y=y, Z=z,
                                                   gridN = 10,
                                                   profileGridN = 25,
                                                profileLambda = TRUE,
                                                profileARange = TRUE)
    
    ## Predict using the spatial process model fitted
    yp = predict(spmods_sm[[counter]],
                 x=xps,Z=elev)
    
    # convert to a raster
    spmod_rast_sm <- raster::rasterFromXYZ(data.frame(lon=lon,lat=lat,  
                                                   prediction = yp),
                                        crs = raster::crs(elevation)) %>%
      terra::rast()
    
    names(spmod_rast_sm) <- paste0(month_grp, "_",
                                   dpt, "_", sample_year)
    result_sm[[counter]] <- spmod_rast_sm
    
    p<- ggarrange(
      ggarrange(
        spmod_rast_sm %>%
        as.data.frame(xy=TRUE) %>%
        pivot_longer(cols = names(.)[3:length(.)], 
                     values_to = "soil_moisture_pct", 
                     names_to = "month_grp") %>%
        ggplot() +
        geom_raster(aes(x=x, y=y, fill=soil_moisture_pct)) + 
        scale_fill_viridis(option="A") +
        facet_wrap(~month_grp,ncol = 1)+
        geom_sf(data = d %>%
                  mutate(sentek_m = y), aes(color=sentek_m)) +
        coord_sf()+
          theme_void(),
        spmod_rast_sm %>%
          as.data.frame(xy=TRUE) %>%
          pivot_longer(cols = names(.)[3:length(.)], 
                       values_to = "soil_moisture_pct", 
                       names_to = "month_grp") %>%
          ggplot() +
          geom_raster(aes(x=x, y=y, fill=soil_moisture_pct)) + 
          scale_fill_stepsn(breaks = c(-10,0,10,40,50,1000),
                            colors = c("red", "white", "white",
                                       "grey", "grey", "black"))+
          facet_wrap(~month_grp,ncol = 1)+
          geom_sf(data = d %>%
                    mutate(sentek_m = y), aes(color=sentek_m)) +
          coord_sf() +
          theme_void(), nrow=2, ncol=1),
      ggarrange(
        ggplot(as.data.frame(y), aes(x=y)) +
          geom_histogram()+
          ggtitle("SM Input", month_grp),
        ggplot(as_tibble(spmod_rast_sm) %>% dplyr::select(xx=1), aes(x=xx)) +
          geom_histogram()+
          ggtitle("SM predictions",month_grp),
        nrow=2, ncol=1)
      ,nrow=1, ncol=2, widths = c(4,1)
    )
    ggsave(filename = paste0("figs/surfaces/sm/",
                             paste0(month_grp, "_", 
                                    sample_year, "_sm.png")),
           plot= p, width=15, height=10, bg="white")
    
    counter <- counter +1
  }}}

# spei ===================================================
drake_precip_df <- read_csv("data/past_data/ages_input/drake58hru/data/reg_precip.csv",
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

drake_tmean_df <- read_csv("data/past_data/ages_input/drake58hru/data/reg_tmean.csv",
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

spei_13 <- fulldf %>%
  filter(ym == "03.2013") %>%
  dplyr::select(spei06, spei12, spei24)

spei_13j <- fulldf %>%
  filter(ym == "06.2013") %>%
  dplyr::select(spei06, spei12, spei24)

spei_14 <- fulldf %>%
  filter(ym == "03.2014") %>%
  dplyr::select(spei06, spei12, spei24)

spei_14j <- fulldf %>%
  filter(ym == "06.2014") %>%
  dplyr::select(spei06, spei12, spei24)


# getting sm and st all together ===============================================
soil_temp_rasts_13 <- terra::rast(result)[[c(1,3,5,7,9)]]
names(soil_temp_rasts_13) <- names(soil_temp_rasts_13) %>%
  str_replace_all("post_seed", "soil_temp_post") %>%
  str_replace_all("pre_seed", "soil_temp_pre") %>%
  str_remove_all("_2013")
soil_temp_rasts_14 <- terra::rast(result)[[c(2,4,6,8,10)]]
names(soil_temp_rasts_14) <- names(soil_temp_rasts_14) %>%
  str_replace_all("post_seed", "soil_temp_post") %>%
  str_replace_all("pre_seed", "soil_temp_pre") %>%
  str_remove_all("_2014")

soil_moist_rasts_13 <- terra::rast(result_sm)[[c(1,2,5,6,9,10,13,14,17,18)]]
names(soil_moist_rasts_13) <- names(soil_moist_rasts_13) %>%
  str_remove_all("_2013") %>%
  str_replace_all("pre", "soil_moisture_pre") %>%
  str_replace_all("post", "soil_moisture_post")
soil_moist_rasts_14 <- terra::rast(result_sm)[[c(3,4,7,8,11,12,15,16,19,20)]]
names(soil_moist_rasts_14) <- names(soil_moist_rasts_14) %>%
  str_remove_all("_2014") %>%
  str_replace_all("pre", "soil_moisture_pre") %>%
  str_replace_all("post", "soil_moisture_post")

mc_files <- Sys.glob("data/microclima*.tif")
air_temp_jf <- terra::rast(mc_files[1])
air_temp_mam <- terra::rast(mc_files[2])
air_temp_son_pre <- terra::rast(mc_files[3])

air_temp_rasts_13 <- c(air_temp_jf$jf_13_tmean, air_temp_mam$mam_13_tmean, air_temp_son_pre$son_12_tmean)
names(air_temp_rasts_13) <- names(air_temp_rasts_13) %>%
  str_replace_all("_13_tmean", "_pre")%>%
  str_replace_all("_12_tmean", "_pre") %>%
  str_c("air_temp_",.)
air_temp_rasts_14 <- c(air_temp_jf$jf_14_tmean, air_temp_mam$mam_14_tmean, air_temp_son_pre$son_13_tmean)
names(air_temp_rasts_14) <- names(air_temp_rasts_14)  %>%
  str_replace_all("_14_tmean", "_pre")%>%
  str_replace_all("_13_tmean", "_pre") %>%
  str_c("air_temp_",.)

veg_plot_locations <- st_read("data/sampled_centroids.gpkg") %>% 
  mutate(terra::extract(soil_moist_rasts_13, vect(.), cells=TRUE)) #%>%
 # dplyr::select(cell)

shrub_plots <- soil_c %>% filter(strip_type == "shru") %>% pull(plot)
herb_plots <- soil_c %>% filter(strip_type == "herb") %>% pull(plot)

shrub_clm<- veg_plot_locations %>%
              filter(X2012Flag %in% shrub_plots) %>%
  mutate(terra::extract(c(soil_temp_rasts_13,soil_moist_rasts_13), 
               vect(.)),
         terra::extract(c(air_temp_rasts_13), 
                        vect(.))) %>%
  dplyr::select(-ID) %>%
  mutate(march_spei06 = pull(spei_13, spei06),
         march_spei12 = pull(spei_13, spei12),
         march_spei24 = pull(spei_13, spei24),
         june_spei06 = pull(spei_13j, spei06),
         june_spei12 = pull(spei_13j, spei12),
         june_spei24 = pull(spei_13j, spei24))

herb_clm<- veg_plot_locations %>%
  filter(X2012Flag %in% herb_plots) %>%
  mutate(terra::extract(c(soil_temp_rasts_14,soil_moist_rasts_14), 
                        vect(.)),
         terra::extract(c(air_temp_rasts_14), 
                        vect(.))) %>%
  dplyr::select(-ID) %>%
  mutate(march_spei06 = pull(spei_14, spei06),
         march_spei12 = pull(spei_14, spei12),
         march_spei24 = pull(spei_14, spei24),
         june_spei06 = pull(spei_14j, spei06),
         june_spei12 = pull(spei_14j, spei12),
         june_spei24 = pull(spei_14j, spei24))

# mutate(twi = terra::extract(rast(twi_stack$twi), vect( st_as_sf(.)))[,2]) %>%

# plot(soil_moist_rasts_13[[1]]); plot(veg_plot_locations, add=T)
xdata <- bind_rows(shrub_clm, herb_clm) %>%
  mutate(terra::extract(rast(twi_stack), vect(.))) %>%
  dplyr::select(-ID, -X2001Flag, plot = X2012Flag,
                -cell, -UTME, -UTMN, -n, -dem2018, ) %>%
  st_set_geometry(NULL) %>%
  left_join(soil_c) %>%
  mutate(soil_texture = str_remove_all(soil_unit_name, 
                                       " [:digit:] to [:digit:] percent slope") %>%
           str_remove_all("Kim ") %>%
           str_remove_all("Colby ") %>%
           str_remove_all("Wagonwheel ")) 

