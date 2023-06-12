# Test model data prep for Drake farm
# Ana Miller-ter Kuile
# June 11, 2023

#this script preps data for a test model with real dta
#for the Drake Farm SAM JSDM

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("jagsUI",
                  "coda",
                  "mcmcplots") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

Sys.time()
# Load data ---------------------------------------------------------------

data <- readRDS("/scratch/atm234/drake_farm/inputs/input_data.RDS")


# Compile data ------------------------------------------------------------

data_list <- list(y = data$y,
                  n.species = data$n.species,
                  n.sites = data$n.sites, 
                  n.lv = data$n.lv,
                  n.lag = data$n.lag,
                  n.groups = data$n.groups,
                  GroupID = data$GroupID,
                  SoilM = data$SoilM)

# Run model ---------------------------------------------------------------

# Parameters to track -----------------------------------------------------

params <- c("lv.coef",
            "a0",
            "a1",
            "sig",
            "wA")

# Run Model ---------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file =  "/scratch/atm234/drake_farm/inputs/jags_antecedent_jsdmodel_group_test.R",
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)


saveRDS(mod,
        file = "/scratch/atm234/drake_farm/outputs/model_output_june12.RDS")

Sys.time()


# Convergence check -------------------------------------------------------

mcmcplot(mod$samples,
         dir = '/scratch/atm234/drake_farm/outputs/mcmcplots/')

raf <- raftery.diag(mod$samples)

names <- rownames(raf[[1]]$resmatrix)
ch1 <- raf[[1]]$resmatrix[,2]
ch2 <- raf[[2]]$resmatrix[,2]
ch3 <- raf[[3]]$resmatrix[,2]

raf_all <- as.data.frame(cbind(names, 
                               ch1, ch2, ch3)) %>%
  mutate(ch1 = as.numeric(ch1),
         ch2 = as.numeric(ch2),
         ch3 = as.numeric(ch3)) %>%
  pivot_longer(ch1:ch3,
               names_to = "chain",
               values_to = 'iterations') 

# ggplot(raf_all, aes(x = iterations/3)) +
#   geom_histogram() 

raf_all %>%
  summarise(iterations_90 = quantile(iterations, 
                                     probs = 0.9, 
                                     na.rm = T)/3,
            iterations_95 = quantile(iterations,
                                     probs = 0.95,
                                     na.rm = T)/3,
            max = max(iterations, 
                      na.rm = T)/3)


bu1 <- raf[[1]]$resmatrix[,1]
bu2 <- raf[[2]]$resmatrix[,1]
bu3 <- raf[[3]]$resmatrix[,1]

burn <- as.data.frame(cbind(names, bu1, bu2, bu3)) %>%
  mutate(bu1 = as.numeric(bu1),
         bu2 = as.numeric(bu2),
         bu3 = as.numeric(bu3)) %>%
  filter(!str_detect(names, "z")) %>%
  filter(!str_detect(names, "wA")) %>%
  filter(!str_detect(names, "wB")) %>%
  pivot_longer(bu1:bu3,
               names_to = "chain",
               values_to = 'iterations') 

burn %>%
  summarise(max(iterations, na.rm = T))


