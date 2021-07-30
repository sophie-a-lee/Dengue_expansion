###############################################################
####                                                       ####
####     Script to fit spatio-temporal additive models     ####
####                                                       ####
###############################################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, sf, mgcv, ROCR)


#### Load data ####
## Yearly data
df_year <- fread("data/df_model.csv")


#### Create model dataset ####
## Create model variables
df_model <- df_year %>%
  mutate(# Rearrange REGIC to set local centre as reference
         regic07_relevel = -level07_acpnum + 6,
         regic07 = factor(regic07_relevel, levels = 1:5,
                          labels = c("Local centre",
                                     "Zone centre",
                                     "Sub-regional centre",
                                     "Regional capital",
                                     "Metropolis")),
         regic18_relevel = -level18_num + 6,
         regic18 = factor(regic18_relevel, levels = 1:5,
                          labels = c("Local centre",
                                     "Zone centre",
                                     "Sub-regional centre",
                                     "Regional capital",
                                     "Metropolis")),
         regic_comb = factor(ifelse(year %in% 2001:2010, regic07_relevel, 
                                    regic18_relevel),
                             levels = 1:5,
                             labels = c("Local centre",
                                        "Zone centre",
                                        "Sub-regional centre",
                                        "Regional capital",
                                        "Metropolis")),
         # Use urban proportion, 00 for 2001 - 2009 + 10 for 2010 - 2020 (unless mising in 2000)
         urban_prpn = ifelse(!is.na(urban00) & urban00 != 0, 
                             ifelse(year %in% 2001:2009, (urban00/100), 
                                    (urban10/100)), (urban10/100)),
         prior_outbreak = factor(prior_outbreak, levels = 0:1,
                                 labels =  c("No", "Yes")),
         outbreak_prevyr = factor(outbreak_lag1, levels = 0:1,
                                  labels =  c("No", "Yes"))) %>%
  dplyr::select(municip_code_ibge, year, lon, lat, outbreak, outbreak100,
                regic_comb, urban_prpn, months_suitable.aeg, months_suitable.both,
                prior_outbreak, outbreak_prevyr)
   

#### Fit baseline model (with only spatio-temporal smooths) ####
model_base <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                    ti(lon, lat, year, d = c(2, 1), k = c(100, 8)),
                  data = df_model, family = binomial, method = "REML")

# summary(model_base)
# gam.check(model_base)
# saveRDS(model_base, file = "output/model_base.rds")


#### Fit full model (with urbanisation, REGIC, months suitable and prior outbreak) ####
model_full <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                    ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                    urban_prpn + region_comb + months_suitable.both + prior_outbreak,
                  data = df_model, family = binomial, method = "REML")

# summary(model_full)
# gam.check(model_full)
# saveRDS(model_full, file = "output/model_full.rds")


#### Fit full model using aegypti cut-off (sensitivity analysis) ####
model_aeg <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                    ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                    urban_prpn + region_comb + months_suitable.aeg + prior_outbreak,
                  data = df_model, family = binomial, method = "REML")

# summary(model_aeg)
# gam.check(model_aeg)
# saveRDS(model_aeg, file = "output/model_aeg.rds")


#### Fit full model using prior outbreak from previous year only (sensitivity analysis) ####
model_outbreak <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                   ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                   urban_prpn + region_comb + months_suitable.both + outbreak_prevyr,
                 data = df_model, family = binomial, method = "REML")

# summary(model_outbreak)
# gam.check(model_outbreak)
# saveRDS(model_outbreak, file = "output/model_outbreak.rds")

#### Fit outbreak100 model (sensitivity analysis) ####
model_full100 <- bam(outbreak100 ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                       ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                       urban_prpn + region_comb + months_suitable.both + prior_outbreak,
                  data = df_model, family = binomial, method = "REML")

# summary(model_full100)
# gam.check(model_full100)
# saveRDS(model_full100, file = "output/model_full100.rds")


#### Remove each covariate in turn (sensitivity analysis) ####
## Months suitable
model_clim <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                       ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                       urban_prpn + region_comb + prior_outbreak,
                     data = df_model, family = binomial, method = "REML")
# saveRDS(model_clim, file = "output/model_clim.rds")


## Prior outbreak
model_prior <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                     ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                     urban_prpn + region_comb + months_suitable.both,
                   data = df_model, family = binomial, method = "REML")

# saveRDS(model_prior, file = "output/model_prior.rds")


## Urbanisation
model_urb <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                   ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                   region_comb + months_suitable.both + prior_outbreak,
                 data = df_model, family = binomial, method = "REML")

# saveRDS(model_urb, file = "output/model_urb.rds")


model_regic <- bam(outbreak ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                     ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                     urban_prpn + months_suitable.both + prior_outbreak,
                    data = df_model, family = binomial, method = "REML")

# saveRDS(model_regic, file = "output/model_regic.rds")



#### Compare full and baseline models with AIC, AUC and Brier score ####
## Table S2
model_checks <- function(model, outbreak_obs = df_model$outbreak) {
  
  # Print AIC
  print(paste0("AIC = ", AIC(model)))
  
  # Get the predicted probabilities for each sample
  model.pred <- predict(model, type = "response")
  
  # Brier score 
  brier <- mean((model.pred - outbreak_obs)^2)
  print(paste0("Brier score = ", brier))
  
  # Create ROC object for ROCR package
  rp <- prediction(as.numeric(model.pred), as.numeric(outbreak_obs))
  
  # Calculate area under the curve
  auc <- performance(rp, "auc")@y.values[[1]]
  print(paste0("AUC = ", auc))
  
  # plot ROC curve
  plot(performance(rp, "tpr", "fpr") )
  
}

## Baseline model checks
model_checks(model_base)


## Full model checks
model_checks(model_full)



