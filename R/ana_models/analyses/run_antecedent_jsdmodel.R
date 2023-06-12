# Simulated data for Drake farm project
# Ana Miller-ter Kuile
# May 2, 2023

#this script runs a model with simulated data
#for the drake farm jsdm with antecedent climate

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  "R2jags", #jags wrapper
                  "coda",
                  "mcmcplots",
                  "corrplot") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data list ----------------------------------------------------------


n.species <- 20
n.sites <- 5
n.lag <- 3
y <- matrix(rbinom(n.species*n.sites, 1, 0.5), n.species, n.sites)
Tmax <- matrix(rnorm(n.sites*n.lag, mean = 0, sd  =1), n.sites, n.lag)
#Tmax <- rnorm(n.sites)
PPT <- matrix(rnorm(n.sites*n.lag, mean = 0, sd = 1), n.sites, n.lag)
#PPT <- rnorm(n.sites)
n.lv <- (n.species/2)
GroupID <- c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5))
n.groups = length(unique(GroupID))

data_list <- list(n.species = n.species,
                  n.sites = n.sites,
                  n.lag = n.lag,
                  n.groups= n.groups,
                  y = y,
                  Tmax = Tmax,
                  PPT = PPT,
                  n.lv = n.lv,
                  GroupID = GroupID)



# Parameters to track -----------------------------------------------------

params <- c("lv.coef",
            "a0",
            "a1",
            "a2",
            "sig",
            "wA",
            'wB')

# Run Model ---------------------------------------------------------------

model <- here("R", "ana_models", "models",
              "jags_antecedent_jsdmodel_group.R")

mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)


#for convergence of parameters to track only
params <- c("a0",
            "a1",
            "a2",
            "wA",
            'wB')

mod2 <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)

mcmcplot(mod2$samples)



# Correlation matrix ------------------------------------------------------

#calculate the correlation matrix from the latent variables
#pull out the coefficients by species pair
lv.coef<-mod$sims.list$lv.coef
#make this an array of the number of coefficients,
#number of species, number of species
cmall<-array(NA,dim=c(dim(lv.coef)[1],n.species,n.species))
#get the variance explained
eps.res<-apply(lv.coef,c(1,2),function(x)1-sum(x^2))
#get the per iteration value for the correlation
for(i in 1:dim(lv.coef)[1]){ #for each mcmc sample
  cmall[i,,]<-cov2cor(tcrossprod(lv.coef[i,,]) + diag(apply(eps.res,2,mean)))
}
#get the mean of that value to plot in the correlation matrix
cmest<-apply(cmall,c(2,3),mean)

#plot correlation matrix
corplot<-cmest
corrplot(corplot,method="circle",order="hclust",type="lower")
#corrplot(corplot,method="color",order="original",type="lower",outline="black",tl.pos="d",tl.col="black")
#corrplot(corplot,method="number",order="original",type="upper",col="black",addgrid.col="black", add=T,tl.pos="d",tl.col="black", cl.pos="n")
#plot large number of species
#corrplot(corplot,method="color",order="FPC",type="lower",outline="black",tl.pos="ld",tl.col="black",tl.cex=0.6,tl.srt = 45)

