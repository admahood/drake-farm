# microclima for drake
library(tidyverse)
library(microclima)
library(raster)
# devtools::install_github("mrke/NicheMapR")
library(NicheMapR)

drake_dem <- raster("data/dem2018.tif")

mam13<-microclima::runauto(r = drake_dem,
                    dstart = "01/03/2013",
                    dfinish = "31/05/2013",
                    hgt = 0.1,
                    l = 1, 
                    x = 1)

mam14<-microclima::runauto(r = drake_dem,
                           dstart = "01/03/2014",
                           dfinish = "31/05/2014",
                           hgt = 0.1,
                           l = 1, 
                           x = 1)

# 1.4 GB -- not sustainable to save the whole thing, only takes 20-ish minutes to create
# save(mam13, mam14, file = "data/microclima_mam.Rda")

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
save(stk, file = "data/microclima_mam.Rda")

# winter temps

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
save(stk, file = "data/microclima_jf.Rda")


# antecedent fall temps

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
save(stk, file = "data/microclima_son_pre.Rda")

# spei?

#install.packages("SPEI")
library(SPEI)
?SPEI::spei()
