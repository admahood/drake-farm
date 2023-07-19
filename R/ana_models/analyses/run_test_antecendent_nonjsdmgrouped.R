# Test model data prep for Drake farm
# Ana Miller-ter Kuile
# June 11, 2023

#this script preps data for a test model with real dta
#for the Drake Farm SAM JSDM

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

# Load data ---------------------------------------------------------------

data <- read.csv(here("data",
                      "sam_y.csv"))

species <- read.csv(here('data', 
                         'sam_species_groups.csv'))

covs <- read.csv(here('data',
                      'sam_covariates.csv'))

#Questions for Adam:
#is the order of rows in y data and covariate list the same?

#To DO:
#Figure out random effects structure (probably something related to
# subplots and plots, and potentially strip)

# Counts ------------------------------------------------------------------
#number of species
n.species <- ncol(data)

#number of sites
n.sites <- nrow(data)

#number of levels for covariance matrix
#n.lv <- n.species/2

# Species groups ----------------------------------------------------------

species2 <- species %>%
  distinct(group, sam_grp)

sp_group <- data %>%
  rownames_to_column(var = "sampleID") %>%
  pivot_longer(ia_bromus:sece,
               names_to = "species_code",
               values_to = "presence") %>%
  distinct(species_code) %>%
  left_join(species2, by = c("species_code" = "group")) %>%
  mutate(sam_grp = case_when(sam_grp == "AIG" ~ "IAG",
                             sam_grp %in% c("NPF", "NAF") ~ "NF",
                             TRUE ~ sam_grp)) %>%
  distinct(species_code, sam_grp)

#number of groups of species
n.groups <- length(unique(sp_group$sam_grp))

#vector of length of the number of species with Group ID :
GroupID <- as.numeric(as.factor(sp_group$sam_grp))
  
# Climate data prep -------------------------------------------------------

covs2 <- covs %>%
  dplyr::select(strip_type, plot,
                strip_number.x, plot_raw,
                soil_moisture_pre_jf_30cm, 
                soil_moisture_pre_ma_30cm,
                soil_moisture_pre_son_30cm,
                soil_moisture_post_jja_30cm,
                soil_moisture_post_son_30cm) %>%
  rownames_to_column(var = "sampleID") 
  
#matrix of site rows x time periods of n.lags columns
#need to scale all relative to each other, I think?? but maybe check with Kiona
SoilM <- covs2 %>%
  pivot_longer(cols = soil_moisture_pre_jf_30cm:
                 soil_moisture_post_son_30cm,
               names_to = "lag",
               values_to = "soil_moisture") %>%
  mutate(soil_moisture = scale(soil_moisture)) %>%
  pivot_wider(names_from = "lag",
              values_from = "soil_moisture") %>%
  dplyr::select(sampleID, 
                soil_moisture_post_son_30cm,
                soil_moisture_post_jja_30cm,
                soil_moisture_pre_son_30cm,
                soil_moisture_pre_ma_30cm,
                soil_moisture_pre_jf_30cm) %>%
  column_to_rownames(var = "sampleID") %>%
  as.matrix()

#number of climate lags
n.lag <- ncol(SoilM)

# Response data -----------------------------------------------------------

#y is matrix of species (rows) x sites (columns)
y <- as.matrix(t(data))


# Compile data ------------------------------------------------------------

data_list <- list(y = y,
                  n.species = n.species,
                  n.sites = n.sites, 
                  #n.lv = n.lv,
                  n.lag = n.lag,
                  n.groups = n.groups,
                  GroupID = GroupID,
                  SoilM = SoilM)

#saveRDS(data_list, (here("R", "ana_models", 
#                       "monsoon",
#                       "input_data.RDS")))

# Run model ---------------------------------------------------------------

# Parameters to track -----------------------------------------------------

params <- c("a0",
            "a1",
            "wA")

# Run Model ---------------------------------------------------------------

model <- here("R", "ana_models",
              "models",
              "jags_antecedent_nonjsdmodel_group_test.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)
Sys.time()


mcmcplot(mod$samples)
 