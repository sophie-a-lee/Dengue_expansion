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
  mutate(# Recode REGIC to set local centre as reference
         regic_comb = factor(ifelse(year %in% 2001:2010, regic07_relevel, 
                                    regic18_relevel),
                             levels = 1:5,
                             labels = c("Local centre",
                                        "Zone centre",
                                        "Sub-regional centre",
                                        "Regional capital",
                                        "Metropolis")))
   

#### Fit baseline model (with only spatio-temporal smooths) ####
model_base <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                    ti(lon, lat, year, d = c(2, 1), k = c(100, 8)),
                  data = df_model, family = binomial, method = "REML")

# summary(model_base)
# gam.check(model_base)
# saveRDS(model_base, file = "output/model_base.rds")


#### Fit full model (with urbanisation, REGIC, months suitable and prior outbreak) ####
model_full <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                    ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                    urban_prpn + regic_comb + months_suitable.both + prior_outbreak_fix,
                  data = df_model, family = binomial, method = "REML")

# summary(model_full)
# gam.check(model_full)
# saveRDS(model_full, file = "output/model_full.rds")


#### Add the number of extremely wet months ####
model_wet <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                 ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                 urban_prpn + regic_comb + months_suitable.both + 
                 prior_outbreak_fix + months_wet,
                 data = df_model, family = binomial, method = "REML")

# summary(model_wet)
# gam.check(model_wet)
# saveRDS(model_wet, file = "output/model_wet.rds")


#### Fit full model using aegypti cut-off (sensitivity analysis) ####
model_aeg <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                   ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                   urban_prpn + regic_comb + months_suitable.aeg + prior_outbreak_fix,
                 data = df_model, family = binomial, method = "REML")

# summary(model_aeg)
# gam.check(model_aeg)
# saveRDS(model_aeg, file = "output/model_aeg.rds")


#### Fit outbreak100 model (sensitivity analysis) ####
model_full100 <- bam(outbreak_fix100 ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                       ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                       urban_prpn + regic_comb + months_suitable.both + prior_outbreak_fix100,
                  data = df_model, family = binomial, method = "REML")

# summary(model_full100)
# gam.check(model_full100)
# saveRDS(model_full100, file = "output/model_full100.rds")


#### Fit 75th percentile model (sensitivity analysis) ####
model_full_perc75 <- bam(outbreak_perc75 ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                         ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                         urban_prpn + regic_comb + months_suitable.both + prior_outbreak_perc75,
                       data = df_model, family = binomial, method = "REML")

# summary(model_full_perc75)
# gam.check(model_full_perc75)
# saveRDS(model_full_perc75, file = "output/model_full_perc75.rds")


#### Add each covariate in turn (sensitivity analysis) ####
## Months suitable
model_temp <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                    ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                    months_suitable.both,
                  data = df_model, family = binomial, method = "REML")
# saveRDS(model_clim, file = "output/model_clim.rds")


## Prior outbreak
model_prior <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                     ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                     prior_outbreak_fix,
                   data = df_model, family = binomial, method = "REML")

# saveRDS(model_prior, file = "output/model_prior.rds")


## Urbanisation
model_urb <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                   ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                   urban_prpn,
                 data = df_model, family = binomial, method = "REML")

# saveRDS(model_urb, file = "output/model_urb.rds")


model_regic <- bam(outbreak_fix ~ s(lon, lat, k = 400) + s(year, k = 10) + 
                     ti(lon, lat,  year, d = c(2, 1), k = c(100, 8)) +
                     regic_comb,
                    data = df_model, family = binomial, method = "REML")

# saveRDS(model_regic, file = "output/model_regic.rds")



