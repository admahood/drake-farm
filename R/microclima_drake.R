# microclima for drake

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

# spei?

#install.packages("SPEI")
library(SPEI)
?SPEI::spei()
