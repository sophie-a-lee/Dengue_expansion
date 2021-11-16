##################################################
####                                          ####
####     Script to load packages and data     ####
####                                          ####
##################################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, sf, geofacet, cowplot, geobr,
               ggthemes, gghalves, geobr, ggrepel, devtools, mgcv,
               mgcViz, pROC)


# Overcome issues with ggdist and sf packages
devtools::install_github("mjskay/ggdist")
library(ggdist)
sf_use_s2(FALSE)


#### Load datasets ####
## Monthly epidemiological data for each region
df_month <- fread("data/dir_month_region.csv")


## Monthly climate data aggregated by state
df_climate <- fread("data/climate_month_state.csv")


## Yearly epi, climate + socioeconomic data for each municipality
df_year <- fread("data/df_model.csv")


## Create dataset used to  fit models
df_model <- df_year %>%
  # Relevel REGIC categories to set zone centres as the reference
  mutate(regic07_relevel = -level07_acpnum + 6,
         regic18_relevel = -level18_num + 6,
         # Use 2007 levels 2001 - 2009 and 2018 level 2010 - 2020
         regic_comb = factor(ifelse(year %in% 2001:2009, 
                                    regic07_relevel, regic18_relevel),
                             levels = 1:5,
                             labels = c("Local centre",
                                        "Zone centre",
                                        "Sub-regional centre",
                                        "Regional capital",
                                        "Metropolis")),
         # Use 2000 urban data for 2001 - 2009 and 2010 data for 2010 - 2020
         urban = ifelse(urban00 != 0 & !is.na(urban00) & year %in% 2001:2009, 
                        urban00, urban10),
         # Convert % to proportion for model
         urban_prpn = urban/100)


## Load municipality shapefile (from IBGE, 2010)
shp <- read_municipality()

## File to convert raw  shapefile to 'parent municipality' shapefile
# Combines municipalities founded since 2001 with parent municipalities
parent_conv <- fread("data/parent_municip_conv.csv")


## Aggregate shapefile to parent municipalities
shp_parent <- left_join(shp, parent_conv,
                        by = c("code_muni" = "municip_code_ibge")) %>%
  # Remove lakes included in shape but not data
  filter(!is.na(municip_parent_name)) %>%
  group_by(municip_parent_code, municip_parent_name) %>%
  summarise() %>%
  ungroup() %>%
  # Rename variables with municipality code and name
  rename(municip_code_ibge = municip_parent_code,
         municip_name = municip_parent_name) %>%
  mutate(municip_code = as.numeric(substr(municip_code_ibge, 1, 6)))


