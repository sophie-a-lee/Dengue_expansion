##########################################
####                                  ####
####     Script to visualise data     ####
####                                  ####
##########################################


## Load packages and data 
source("00_load_packages_data.R")


## Load region & state shapefiles for maps
shp_region <- read_region() %>%
  # Add English region names
  mutate(region_name = recode(code_region,
                              `5` = "Centre-West",
                              `1` = "North",
                              `2` = "Northeast",
                              `3` = "Southeast",
                              `4` = "South"))

shp_state <- read_state()


## State grid to produce faceted plots
grid_br <- fread("data/br_states_grid.csv", encoding = "Latin-1")


## Breaks and labels for plots 
year_brk <- (seq(from = 0, to = length(unique(df_month$year))-2, 
                 by = 2)*12) + 1

year_lab <- as.character(seq(from = 2001, 
                             to = max(df_month$year)-1, by = 2))


## Set up colour palette for regions and REGIC levels
region_col <- c("#ef476f", "#FFBA08", "#06d6a0", "#118ab2", "#073b4c")

connect_cols <- c("#eaac8b", "#e56b6f", "#b56576", "#6d597a", "#355070")


#### Maps showing regions and state (Figure S1) ####
region_map <- ggplot(data = shp_region) +
  geom_sf(aes(fill = region_name), lwd = .05) +
  scale_fill_manual(values = region_col) +
  theme_void() +
  theme(legend.position = "none")

ggsave(region_map + 
         # Add region labels
         geom_sf_label(aes(label = region_name)), 
       filename = "output/region_map.png")


state_map <- ggplot(data = shp_state) +
  geom_sf(lwd = .1, fill = "#cbc0d3") + 
  # Adds state labels, jitters to avoid overlap
  ggrepel::geom_label_repel(
    aes(label = abbrev_state, geometry = geom),
    stat = "sf_coordinates") +
  theme_void() 

ggsave(state_map, filename = "output/state_map.png")


#### Climate plots ####
## State monthly tmean geofacet  (Figure S2)
# Combine climate data with state grid
df_climate_grid <- df_climate %>%
  full_join(., grid_br, by = c("state_code" = "code_num")) 


tmean_heat <- ggplot(data = df_climate_grid, 
                     aes(x = month, y = year, fill = tmean)) +
  geom_raster() +
  ylab("Year") + 
  xlab("Month") + 
  scale_fill_distiller(name = "Mean temp (°C)", palette = "YlOrRd",
                       direction = 1) + 
  scale_x_continuous(breaks = c(1, 4, 7, 10), 
                     labels = c("Jan", "Apr", "Jul", "Oct"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(2000, 2020, by = 5), expand = c(0, 0)) +
  theme_bw() + 
  # Plot separate raster per state in the shape of Brazil (using grid_br)
  facet_geo(~name, grid = grid_br) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15))

ggsave(tmean_heat, filename = "output/tmean_heat.png",
       height = 20, width = 20)


## Number of months suitable for transmission (Figure S3)
# Create df with mean number of months with suitable temp by decade
df_suitable <- df_year %>%
  # Create factor for decade 1 (2001 - 2010) + decade 2 (2011 - 2020)
  mutate(decade = cut(year, breaks = c(-Inf, 2010, Inf),
                      labels = c("dec1", "dec2"))) %>%
  group_by(municip_code_ibge, decade) %>%
  # Return mean number of months per decade per municipality
  summarise(months_suitable = mean(months_suitable.both)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


# Plot average number of months with temp suitable per decade
suitable_map <- ggplot(data = df_suitable) +
  geom_sf(aes(fill = months_suitable), lwd = 0) +
  scale_fill_viridis_c(name = "Months suitable \nper year",
                       direction = -1, option = "B") +
  expand_limits(fill = c(0, 12)) +
  facet_wrap(~decade) +
  coord_sf(datum = NA) +
  theme_void() 

ggsave(suitable_map, filename = "output/era_suit.png")


## Difference in months suitable (Figure 2)
# Convert decade df into wide format (variable per decade)
df_suitable_diff <- spread(st_drop_geometry(df_suitable), 
                           decade, months_suitable) %>%
  left_join(., shp_parent, by = "municip_code_ibge") 


# Find number difference in mean suitable months per decade
df_suitable_diff <- mutate(df_suitable_diff,
                           suitable_diff = dec2 - dec1,
                           # For plot (no difference = NA, then set NAs to white)
                           suitable_diff_na = 
                             ifelse(suitable_diff == 0, NA, suitable_diff)) %>%
  # Remove FdN to avoid extra wide map
  filter(municip_code_ibge != 2605459) %>%
  st_as_sf()


suitable_diff <- ggplot(data = df_suitable_diff) +
  geom_sf(aes(fill = suitable_diff_na), lwd = .05) +
  scale_fill_distiller(name = "Change in \nmonths suitable",
                       palette = "PiYG", na.value = "white") +
  expand_limits(fill = c(1, -1)) +
  coord_sf(datum = NA) +
  theme_void()

ggsave(suitable_diff, filename = "output/era_suit_diff.png")


## State monthly scPDSI set out in shape of Brazil (using grid_br)  (Figure M1)
pdsi_heat <- ggplot(data = df_climate_grid, 
                     aes(x = month, y = year, fill = pdsi)) +
  geom_raster() +
  ylab("Year") + 
  xlab("Month") + 
  scale_fill_distiller(name = "scPDSI", palette = "BrBG",
                       direction = 1) + 
  scale_x_continuous(breaks = c(1, 4, 7, 10), 
                     labels = c("Jan", "Apr", "Jul", "Oct"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(2000, 2020, by = 5), expand = c(0, 0)) +
  expand_limits(fill = c(-5,5)) +
  theme_bw() + 
  # Plot separate raster per state in the shape of Brazil (using grid_br)
  facet_geo(~name, grid = grid_br) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15))

ggsave(pdsi_heat, filename = "output/pdsi_heat.png",
       height = 20, width = 20)


## Create df with average number of months 'extremely wet' per decade per municipality
df_wet <- df_year %>%
  # Create factor for decade 1 (2001 - 2010) + decade 2 (2011 - 2020)
  mutate(decade = cut(year, breaks = c(-Inf, 2010, Inf), 
                      labels = c("dec1", "dec2"))) %>%
  group_by(municip_code_ibge, decade) %>%
  summarise(months_wet = mean(months_wet)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


## Difference in months extremely wet (Figure M2)
# Convert decade df into wide format (variable per decade)
df_wet_diff <- spread(st_drop_geometry(df_wet), 
                      decade, months_wet) %>%
  left_join(., shp_parent, by = "municip_code_ibge") 


# Find number difference in mean number of wet months per decade
df_wet_diff <- mutate(df_wet_diff,
                      wet_diff = dec2 - dec1,
                      # Set no difference to NA for plot
                      wet_diff_na = ifelse(wet_diff == 0, NA, wet_diff)) %>%
  # Remove FdN to avoid extra wide map
  filter(municip_code_ibge != 2605459) %>%
  st_as_sf()


wet_diff <- ggplot(data = df_wet_diff) +
  geom_sf(aes(fill = wet_diff_na), lwd = .05) +
  scale_fill_distiller(name = "Change in \nmonths wet",
                       palette = "PiYG", na.value = "white") +
  expand_limits(fill = c(1, -1)) +
  coord_sf(datum = NA) +
  theme_void()

ggsave(wet_diff, filename = "output/wet_diff.png")



#### Census plots ####
## Urbanisation maps (Figure S4)
# Create df with census and REGIC variables (select 2010 only to avoid duplicates)
df_census <- df_year %>%
  filter(year == 2010) %>% 
  left_join(., shp_parent, by = c("municip_code_ibge")) %>%
  mutate(regic07 = factor(level07_acpnum, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre")),
         regic18 = factor(level18_num, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre"))) %>%
  st_as_sf()


# Plot % urbanisation from 2000 census
urban00_map <- ggplot(data = df_census) +
  geom_sf(aes(fill = urban00), lwd = 0) +
  scale_fill_continuous_tableau(name = "% urbanisation", 
                                palette = "Blue-Teal") +
  # Set legend to 0-100% (to match 2010 map)
  expand_limits(fill = c(0, 100)) +
  theme_void()


# Plot % urbanisation from 2010 census
urban10_map <- ggplot(data = df_census) +
  geom_sf(aes(fill = urban10), lwd = 0) +
  scale_fill_continuous_tableau(name = "% urbanisation", 
                                palette = "Blue-Teal") +
  # Set legend to 0-100% (to match 2000 map)
  expand_limits(fill = c(0, 100)) +
  theme_void()


# Combine urbanisation maps in same plot
urban_maps <- plot_grid(urban00_map + theme(legend.position = "none"),
                        urban10_map + theme(legend.position = "none"),
                        labels = c("2000 census", "2010 census"))


# Extract legend from one of the maps to add to plot
urban_leg <- get_legend(urban00_map)

# Combine maps and legend on same plot
urban_maps <- plot_grid(urban_maps, urban_leg, rel_widths = c(3, .4))

ggsave(urban_maps, filename = "output/urban_decade.png", 
       height = 5, width = 10)


## Urban vs basic services scatterplots (Figure S5)
# Urbanisation vs. piped water
urban_water <- ggplot(data = df_census) +
  geom_point(aes(x = urban10, y = water_network, colour = region_name)) + 
  scale_x_continuous(name = "% urbanisation", expand = c(0, 5)) +
  scale_y_continuous(name = "% with access to piped water") +
  scale_colour_manual(name = "Region", values = region_col) +
  theme_light() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))


ggsave(urban_water, filename = "output/urban_water.png", 
       height = 5, width = 10)

# Return correlation and p-value between urbanisation and piped water access
cor.test(df_census$urban10, df_census$water_network) 
# r = 0.656 [0.641, 0.671] p. < 0.001


# Urbanisation vs. refuse collection
urban_refuse <- ggplot(data = df_census) +
  geom_point(aes(x = urban10, y = total_collected, colour = region_name)) + 
  scale_x_continuous(name = "% urbanisation", expand = c(0, 5)) +
  scale_y_continuous(name = "% with refuse collection") +
  scale_colour_manual(name = "Region", values = region_col) +
  theme_light() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))


ggsave(urban_refuse, filename = "output/urban_refuse.png", 
       height = 5, width = 10)

# Return correlation and p-value between urbanisation and refuse collection
cor.test(df_census$urban10, df_census$total_collected) 
# r = 0.794 [0.784, 0.804] p. < 0.001


#### REGIC plots ####
## Map showing REGIC categories from 2018 (Figure 3)
regic18_map <- ggplot(data = df_census) +
  geom_sf(aes(fill = regic18), lwd = 0) +
  scale_fill_manual(name = "Level of \ninfluence", values = connect_cols) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
                      legend.text = element_text(size = 10))


ggsave(regic18_map, filename = "output/regic18_map.png")


## Histograms of % municipalities in each level per region (Figure S6)
# Create df with percentage of municipalities in each region in each REGIC category (2007)
df_regic07_perc <- st_drop_geometry(df_census) %>%
  group_by(region_name, regic07) %>%
  # Return number of municipalities per region in each category
  summarise(regic07_n = n()) %>%
  # Convert number into a percentage per region
  mutate(regic07_perc = regic07_n/sum(regic07_n) * 100)


regic07_hist <- ggplot(data = df_regic07_perc) +
  geom_bar(aes(x = region_name, y = regic07_perc, fill = regic07), 
           stat = "identity") +
  scale_y_continuous(name = "% of municipalities", expand = c(0, 0)) +
  scale_x_discrete(name = "Region", expand = c(0, 0)) +
  scale_fill_manual(name = "Level of influence",  values = connect_cols) +
  theme_light()


# Create df with percentage of municipalities in each region in each REGIC category (2018)
df_regic18_perc <- st_drop_geometry(df_census) %>%
  group_by(region_name, regic18) %>%
  # Return number of municipalities per region in each category
  summarise(regic18_n = n()) %>%
  # Convert number into a percentage per region
  mutate(regic18_perc = regic18_n/sum(regic18_n) * 100)


regic18_hist <- ggplot(data = df_regic18_perc) +
  geom_bar(aes(x = region_name, y = regic18_perc, fill = regic18), stat = "identity") +
  scale_y_continuous(name = "% of municipalities", 
                     expand = c(0, 0)) +
  scale_x_discrete(name = "Region", expand = c(0, 0)) +
  scale_fill_manual(name = "Level of influence",  values = connect_cols) +
  theme_light()


# Combine histograms from 2007 & 2018 onto same plot
regic_hist <- plot_grid(regic07_hist + theme(legend.position = "none"),
                        regic18_hist + theme(legend.position = "none"))

# Extract legend (same for both histograms)
regic_leg <- get_legend(regic07_hist)

# Add legend to the combined histogram plot
regic_hist <- plot_grid(regic_hist, regic_leg, rel_widths = c(3, .4))

ggsave(regic_hist, filename = "output/regic_hist.png",
       height = 5, width = 15)


## Raincloud plots REGIC vs census variables (Figure S7)
# Urbanisation
regic_urb <- ggplot(df_census, aes(x = regic18, y = urban10, fill = regic18)) +
  # Add density plot, % urbanisation per REGIC category
  stat_halfeye(adjust = .5,
               .width = 0,
               justification = -.2, scale = .6) + 
  # Add boxplot, remove outliers as these are shown in the points below
  geom_boxplot(width = .15, outlier.shape = NA, position = "dodge") +
  # Adds points, jittered as there is a lot of overlap
  geom_half_point(side = "l", range_scale = .4,
                  alpha = .3, position = "dodge") +
  labs(x = "Level of influence", y = "% urbanisation") +
  scale_fill_manual(values = connect_cols) +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))


ggsave(regic_urb, filename = "output/regic_urb_rain.png",
       width = 20, height = 10)


# Access to water network
regic_water <- ggplot(df_census, aes(x = regic18, y = water_network, 
                                     fill = regic18)) +
  # Add density plot, % access to water network per REGIC category
  stat_halfeye(adjust = .5,
               .width = 0,
               justification = -.2,
               scale = .6) +
  # Add boxplot, remove outliers as these are shown in the points below
  geom_boxplot(width = .15, outlier.shape = NA) +
  # Adds points, jittered as there is a lot of overlap
  geom_half_point(side = "l", range_scale = .4,
                  alpha = .3) +
  labs(x = "Level of influence", y = "% with access to water network") +
  scale_fill_manual(values = connect_cols) +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))


ggsave(regic_water, filename = "output/regic_water_rain.png",
       width = 20, height = 10)


# Refuse collection
regic_refuse <- ggplot(df_census, aes(x = regic18, y = total_collected, 
                                      fill = regic18)) +
  # Add density plot, % refuse collection per REGIC category
  stat_halfeye(adjust = .5,
               .width = 0,
               justification = -.2,
               scale = .6) +
  # Add boxplot, remove outliers as these are shown in the points below
  geom_boxplot(width = .15, outlier.shape = NA) +
  # Adds points, jittered as there is a lot of overlap
  geom_half_point(side = "l", range_scale = .4,
                  alpha = .3) +
  labs(x = "Level of influence", y = "% with refuse collection") +
  scale_fill_manual(values = connect_cols) +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))


ggsave(regic_refuse, filename = "output/regic_refuse_rain.png",
       width = 15, height = 10)


#### REGIC tables ####
## Number of municipalities per region in each REGIC category (Table S1)
df_regic07 <- df_year %>%
  filter(year == 2007) %>% 
  mutate(regic07 = factor(level07_acpnum, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre")))

addmargins(table(df_regic07$region_name, df_regic07$regic07))


df_regic18 <- df_year %>%
  filter(year == 2018) %>% 
  mutate(regic18 = factor(level18_num, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre")))

addmargins(table(df_regic18$region_name, df_regic18$regic18))



#### Dengue case exploratory plots ####
## Regional DIR timeplot (Figure 4)
# Create df with DIR per region per month
df_region_dir <- df_month %>%
  group_by(region_name, time) %>%
  summarise(dengue_total = sum(dengue_cases), 
            population_total = sum(population)) %>%
  mutate(DIR = (dengue_total/population_total) * 10^5) %>%
  ungroup() 


# Time series with DIR per region per month
region_dir_time <-  ggplot(data = df_region_dir) +
  geom_line(aes(x = time, y = DIR, colour = region_name,  group = region_name)) +
  scale_x_continuous(name = "Year", 
                     breaks = year_brk, labels = year_lab, expand = c(0,0)) +
  scale_y_continuous(name = "Dengue incidence rate (per 100,000 residents)", 
                     labels = scales::comma, expand = c(0, 0)) +
  scale_color_manual(name = "Region", values = region_col) +
  theme_light() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


# Extract legend from time series
region_leg <- get_legend(region_dir_time)

# Combine time series with region map
region_line_inset <- ggdraw() +
  draw_plot(region_dir_time + theme(legend.position = "none")) +
  draw_plot(region_map + theme(legend.position = "none"), 
            x = .01, y = .5, width = .5, height = .5) 

# Add legend to same plot
region_line_comp <- plot_grid(region_line_inset, region_leg,
                              rel_widths = c(3, .4))

ggsave(region_line_comp, filename = "output/region_dir_inset.png",
       height = 5, width = 10)


#### Outbreak plots ####
## Number of years with outbreak (Figure 1)
# Create df with number of year where DIR > 300
df_permanence <- df_year %>%
  mutate(DIR_year = (dengue_year/population)*10^5,
         outbreak = ifelse(DIR_year >= 300, 1, 0)) %>%
  group_by(municip_code_ibge) %>%
  summarise(number_outbreaks = sum(outbreak)) %>%
  ungroup() %>%
  # Set regions with no prior outbreak to NA to change colour in the map
  mutate(number_outbreaks_na = ifelse(number_outbreaks == 0, NA, 
                                      number_outbreaks)) %>%
  full_join(., shp_parent, by = "municip_code_ibge") %>%
  # Remove  FdN to avoid extra whitespace
  filter(municip_code_ibge !=2605459) %>%
  st_as_sf()


# Plot map with number of years DIR > 300
number_outbreaks <- ggplot(data = df_permanence) +
  geom_sf(aes(fill = number_outbreaks_na, colour = ""), lwd = 0) +
  scale_fill_gradient_tableau("Purple", 
                              name = "Number of years \nwith outbreaks", 
                              # Change colour for municipalities with no outbreaks
                              na.value = "grey") +
  # Add regions with no outbreaks to the legend
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No outbreaks")) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


ggsave(number_outbreaks, filename = "output/number_outbreaks.png",
       height = 10, width = 10)


## Percentage of municipalities per year experiencing outbreaks (Figure S8)
# DF with number of  municipalities per year
municip_region <- st_drop_geometry(df_census) %>%
  group_by(region_name) %>%
  summarise(n_municip = n()) 


# DF with  percentage of municipalities per region with outbreaks each year
df_region_outbreaks <- df_year %>%
  mutate(DIR_year = (dengue_year/population)*10^5,
         outbreak = ifelse(DIR_year >= 300, 1, 0)) %>%
  group_by(region_name, year) %>%
  # Number of municipalities per region with outbreaks each year
  summarise(total_outbreak = sum(outbreak)) %>%
  ungroup() %>%
  # Combine with total number of municipalities per region
  full_join(., municip_region, by = "region_name") %>%
  mutate(perc_outbreak = (total_outbreak/n_municip) * 100)


# Time series with % municipalities in each region with outbreaks per year
perc_region_outbreak <- ggplot(data = df_region_outbreaks) +
  geom_line(aes(x = year, y = perc_outbreak, colour = region_name,  
                group = region_name)) +
  scale_x_continuous(name = "Year", expand = c(0,0)) +
  scale_y_continuous(name = "% of municipalities with outbreaks", 
                     labels = scales::comma, expand = c(0, 0)) +
  scale_color_manual(name = "Region", values = region_col) +
  theme_light() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


# Add region map to time series
perc_region_inset <- ggdraw() +
  draw_plot(perc_region_outbreak + theme(legend.position = "none")) +
  draw_plot(region_map + theme(legend.position = "none"), 
            x = .01, y = .5, width = .5, height = .5) 

# Combine with legend (from earlier)
perc_region_inset <- plot_grid(perc_region_inset, region_leg,
                               rel_widths = c(3, .4))

ggsave(perc_region_inset, filename = "output/region_perc_inset.png",
       height = 5, width = 10)


## Year of first outbreak (Figure 5)
# DF with first year muncipality had DIR > 300
df_outbreak <- df_year %>%
  mutate(DIR = (dengue_year/population)*10^5) %>%
  # Remove non-outbreak years
  filter(DIR >= 300) %>%
  arrange(municip_code_ibge, year) %>%
  group_by(municip_code_ibge) %>%
  # Return earliest year  each municipality had an outbreak
  mutate(first_fix = first(year)) %>%
  # Remove years after the first outbreak
  filter(year == first_fix) %>%
  ungroup()

# Add year of first outbreak to sf data, add NA for any municipality yet to experience an outbreak
df_outbreak <- df_outbreak[!duplicated(df_outbreak$municip_code_ibge), 
                           c("municip_code_ibge", "first_fix")] %>%
  full_join(., shp_parent, by = "municip_code_ibge") %>%
  # Remove FdN to avoid excessive whitespace
  filter(municip_code_ibge !=2605459) %>%
  # Add variable to show year of first outbreak prior to 2010
  mutate(first_outbreak10 = ifelse(first_fix > 2010, NA, first_fix)) %>%
  st_as_sf()


# Map of first outbreak 2001 - 2020 
first_outbreak20 <- ggplot(data = df_outbreak) + 
  geom_sf(aes(fill = first_fix, colour = ""), lwd = .05) +
  scale_fill_gradient_tableau("Purple", 
                              # Change colour for municipalities with no outbreaks
                              na.value = "grey",
                              name = "Year of \nfirst outbreak",
                              trans = "reverse") +
  expand_limits(fill = c(2001, 2020)) +
  # Add municipalities without outbreaks to legend
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No outbreaks")) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


# Map with year of first outbreak (between 2001 - 2010)
first_outbreak10 <- ggplot(data = df_outbreak) + 
  geom_sf(aes(fill = first_outbreak10, colour = ""), lwd = .05) +
  scale_fill_gradient_tableau("Purple", 
                              # Change colour for municipalities with no outbreaks by 2010
                              na.value = "grey",
                              name = "Year of \nfirst outbreak",
                              trans = "reverse") +
  # Ensure legends are on same scale
  expand_limits(fill = c(2001, 2020)) +
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No outbreaks")) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


# Extract legend from first outbreak map (same for both)
outbreak_leg <- get_legend(first_outbreak20)


# Combine maps from 2001 - 2010 and 2001 - 2020
first_outbreak_decade <- plot_grid(first_outbreak10 + theme(legend.position = "none"),
                                   first_outbreak20 + theme(legend.position = "none"),
                                   labels = c("2001 - 2010", "2001 - 2020"))


# Add legend
first_outbreak_decade <- plot_grid(first_outbreak_decade, outbreak_leg, 
                                   rel_widths = c(3, .4))


ggsave(first_outbreak_decade, filename = "output/first_outbreak_decade.png",
       height = 5, width = 10)


#### Plot 75th percentile outbreak threshold (Figure M) ####
df_cutoff <- df_year %>% 
  group_by(municip_code_ibge) %>% 
  # Calculate the 75th percentile of DIR per municipality
  summarise(perc_75 = quantile(DIR, .75),
            pop_mean = mean(population)) %>% 
  ungroup() %>% 
  # Add minimum cutoff (equivalent to 5 cases) 
  mutate(DIR_5 = (5/pop_mean)*10^5,
         perc75_cutoff = ifelse(perc_75 < DIR_5, DIR_5, perc_75)) %>% 
  left_join(., shp_parent, by = "municip_code_ibge") %>% 
  st_as_sf()


# Plot cutooff threshold for each municipality
perc75_cutoff <- ggplot(data = df_cutoff) +
  geom_sf(aes(fill = perc75_cutoff), lwd = .05) +
  scale_fill_gradient2_tableau("Gold-Purple Diverging", 
                               name = "Outbreak threshold",
                               trans = "log1p",
                               breaks = c(0, 30, 100, 300, 1000, 3000)) +
  theme_void() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

ggsave(perc75_cutoff, filename = "output/perc75_cutoff.png")
