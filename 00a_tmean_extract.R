##################################################################
####                                                          ####
####     Code to extract tmean variables from ERA5 raster     ####
####                                                          ####
##################################################################

#### Load in packages ####
pacman::p_load(tidyverse, data.table, sf, raster, ncdf4, exactextractr,
               lubridate, geobr)


#### Load data ####
### 1. Shapefile 
shp <- read_municipality()


### 2. Convert to parent municipalities to avoid issues with boundary changes
parent_conv <- fread("data/parent_municip_conv.csv")

shp_clean <- left_join(shp, parent_conv, 
                       by = c("code_muni" = "municip_code_ibge")) %>%
  # Remove lakes included in shape but not data
  filter(!is.na(municip_parent_name)) %>%
  group_by(municip_parent_code, municip_parent_name) %>%
  # Join new and parent municipalities
  summarise() %>%
  rename(municip_code_ibge = municip_parent_code,
         municip_name = municip_parent_name) %>%
  st_as_sf() %>%
  # Ensure CRS matches ERA5 
  st_transform(., 4326)

# Number variables in shapefile (for later function)
nshp <- ncol(shp_clean)


### 3. Load ERA5 tmean raster file
temp2m <- raster::brick("data/era5_tmean_200120.nc")


#### Extract monthly tmean for each municipality ####
## Recode layer labels (or Z)
dates <- seq(ymd("20010101"), ymd("20201201"), by = "months")
temp2m <- setZ(temp2m, dates)


## Crop raster file to Brazil 
ras_temp2m <- crop(temp2m, extent(shp_clean), snap = "out")


## Update shape to include small regions partially covered by tiles
ras_shp <- rasterize(shp_clean, ras_temp2m, getCover = T)


## Crop the raster to polygons
ras_temp2m <- mask(ras_temp2m, ras_shp)


## Convert Kelvin to celcius
k_to_c_fun <- function(x) {x - 273.15}
ras_temp2mc <- calc(ras_temp2m, k_to_c_fun)

## Extract mean of each municipality and add to dataset
# Number layers in raster
n1 <- nlayers(ras_temp2mc)

shp_clean[ , (nshp + 1):(nshp + n1)] <- exact_extract(ras_temp2mc, shp_clean, 'mean')


## Convert into data frame
temp2m <- gather(shp_clean, colname, tmean, 
                 names(shp_clean[, nshp + 1])[[1]]:names(shp_clean[, nshp + n1])[[1]]) 

# Add time
temp2m$year <- rep(2001:2020, each = 12*nrow(shp_clean))
temp2m$month <- rep(rep(1:12, each = nrow(shp_clean)), times = 20)


temp2m <- ungroup(temp2m) %>%
  dplyr::select(municip_code_ibge, year, month, tmean)

shp_clean <- shp_clean[ , -c((nshp + 1):(nshp + n1))]


#### Convert into months suitable per year ####
tmean_suitable <- temp2m %>%
  mutate(suitable = ifelse(tmean >= 17.8 & tmean <= 34.5, 1, 0)) %>%
  group_by(municip_code_ibge, year) %>%
  summarise(months_suitable = sum(suitable)) %>%
  ungroup()


## Save data 
fwrite(tmean_suitable, "data/temp_suitable.csv")


