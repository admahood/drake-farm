# microclima for drake
library(tidyverse)
library(microclima)
library(raster)
# devtools::install_github("mrke/NicheMapR")
library(NicheMapR)

drake_dem <- raster("data/dem2018.tif")

# mam ==========================================================================

mam13<-microclima::runauto(r = drake_dem,
                    dstart = "01/03/2013",
                    dfinish = "30/04/2013",
                    hgt = 0.1,
                    l = 1, 
                    x = 1)

mam14<-microclima::runauto(r = drake_dem,
                           dstart = "01/03/2014",
                           dfinish = "30/04/2014",
                           hgt = 0.1,
                           l = 1, 
                           x = 1)


lst <- list()
lst[[1]] <- mam14$tmax
lst[[2]] <- mam14$tmin 
lst[[3]] <- mam14$tmean

lst[[4]] <- mam13$tmax 
lst[[5]] <- mam13$tmin
lst[[6]] <- mam13$tmean 

stk <- lapply(lst,terra::rast) %>% terra::rast()
names(stk) <- c("mam_14_tmax","mam_14_tmin","mam_14_tmean",
                "mam_13_tmax","mam_13_tmin","mam_13_tmean")
terra::writeRaster(stk, "data/microclima_ma.tif")

rm(mam13, mam14)
gc()

# winter temps =================================================================

jf13<-microclima::runauto(r = drake_dem,
                           dstart = "01/01/2013",
                           dfinish = "28/02/2013",
                           hgt = 0.1,
                           l = 1, 
                           x = 1)

jf14<-microclima::runauto(r = drake_dem,
                           dstart = "01/01/2014",
                           dfinish = "28/02/2014",
                           hgt = 0.1,
                           l = 1, 
                           x = 1)

# 1.4 GB -- not sustainable to save the whole thing, only takes 20-ish minutes to create
# save(mam13, mam14, file = "data/microclima_mam.Rda")

lst <- list()
lst[[1]] <- jf13$tmax
lst[[2]] <- jf13$tmin 
lst[[3]] <- jf13$tmean

lst[[4]] <- jf14$tmax 
lst[[5]] <- jf14$tmin
lst[[6]] <- jf14$tmean 

stk <- lapply(lst,terra::rast) %>% terra::rast()
names(stk) <- c("jf_13_tmax","jf_13_tmin","jf_13_tmean",
                "jf_14_tmax","jf_14_tmin","jf_14_tmean")
terra::writeRaster(stk, "data/microclima_jf.tif")

rm(jf13, jf14);gc()

# antecedent fall temps ========================================================

son12<-microclima::runauto(r = drake_dem,
                          dstart = "01/09/2012",
                          dfinish = "30/11/2012",
                          hgt = 0.1,
                          l = 1, 
                          x = 1)

son13<-microclima::runauto(r = drake_dem,
                          dstart = "01/09/2013",
                          dfinish = "30/11/2013",
                          hgt = 0.1,
                          l = 1, 
                          x = 1)

# 1.4 GB -- not sustainable to save the whole thing, only takes 20-ish minutes to create
# save(mam13, mam14, file = "data/microclima_mam.Rda")

lst <- list()
lst[[1]] <- son12$tmax
lst[[2]] <- son12$tmin 
lst[[3]] <- son12$tmean

lst[[4]] <- son13$tmax 
lst[[5]] <- son13$tmin
lst[[6]] <- son13$tmean 

stk <- lapply(lst,terra::rast) %>% terra::rast()
names(stk) <- c("son_12_tmax","son_12_tmin","son_12_tmean",
                "son_13_tmax","son_13_tmin","son_13_tmean")
terra::writeRaster(stk, "data/microclima_son_pre.tif")

rm(son13, son12);gc()

# spei?

#install.packages("SPEI")
library(SPEI)
?SPEI::spei()


# era5 based
# microclima w era5 =============================================================

uid <- "176012"
cds_api_key <- "46a2883c-9c82-4500-a98b-7c629608d48e"

ecmwfr::wf_set_key(user = uid,
                   key = cds_api_key,
                   service = "cds")

req<-mcera5::build_era5_request(ymin = 40.59726526881546,
                                xmin = -104.85857930354454,
                                ymax = 40.614978854837126, 
                                xmax = -104.82755117479745,
                                start_time = as.Date("2020-08-04"),
                                end_time = as.Date("2020-08-14"),
                                outfile_name = "era5_drake")

request_era5(req, uid=uid, out_path = "data/out/era5")

drake_center <- c(40.607243053702604, -104.84138600370632)

point_out <- extract_clim(nc = "data/out/era5/era5_drake_2020.nc", 
                          long = drake_center[2], 
                          lat = drake_center[1],
                          dtr_cor = F,
                          start_time = as.Date("2020-08-04"), 
                          end_time = as.Date("2020-08-14"))

class(point_out) <- "data.frame"
point_out_precip <- extract_precip(nc = "data/out/era5/era5_drake_2020.nc", 
                                   long = drake_center[2], 
                                   lat = drake_center[1],
                                   start_time = as.Date("2020-08-04"), 
                                   end_time = as.Date("2020-08-14"),
                                   convert_daily = TRUE)

library(terra)
r <- microclima::get_dem(lat = drake_center[1], 
                         long = drake_center[2], resolution = 30)
microclima::get_NCEP(drake_center[1], drake_center[2], as.Date("2020-08-08"))->ncep
microe5 <- microclima::runauto(r = dem,
                               dstart = "04/08/2020",
                               dfinish = "14/08/2020",
                               hourlydata = point_out,
                               dailyprecip = point_out_precip,
                               hgt = 2,
                               coastal = FALSE,
                               l = NA, x = NA, weather.elev = "era5",
                               habitat = "Barren or sparsely vegetated")

outfile <- "/project/microclimate/drake_microe5_2020_08.Rda"

save(microe5, file = outfile)
