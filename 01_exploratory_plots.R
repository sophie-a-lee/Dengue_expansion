##########################################
####                                  ####
####     Script to visualise data     ####
####                                  ####
##########################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, sf, geofacet, cowplot, geobr,
               ggthemes, ggdist, gghalves, geobr, ggrepel)
sf_use_s2(FALSE)


#### Load data ####
## Monthly regional data
df_month <- fread("data/dir_month_region.csv")


## Yearly data 
df_year <- fread("data/df_model.csv")


## Monthly climate state data
df_climate <- fread("data/climate_month_state.csv", encoding = "Latin-1")


## Shapefile
# Load raw municipality data
shp <- read_municipality()

# Convert into parent municipalities
parent_conv <- fread("data/parent_municip_conv.csv")

# Aggregate shapefile
shp_parent <- left_join(shp, parent_conv,
                        by = c("code_muni" = "municip_code_ibge")) %>%
  # Remove lakes included in shape but not data
  filter(!is.na(municip_parent_name)) %>%
  group_by(municip_parent_code, municip_parent_name) %>%
  summarise() %>%
  ungroup() %>%
  rename(municip_code_ibge = municip_parent_code,
         municip_name = municip_parent_name) %>%
  mutate(municip_code = as.numeric(substr(municip_code_ibge, 1, 6)),
         area_sqkm = as.numeric(st_area(.)/10^6))


## Load region & state shapefiles for maps
shp_region <- read_region() %>%
  mutate(region_name = recode(code_region,
                              `5` = "Centre-West",
                              `1` = "North",
                              `2` = "Northeast",
                              `3` = "Southeast",
                              `4` = "South"))

shp_state <- read_state()


## Grid for faceted plots
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

ggsave(region_map + geom_sf_label(aes(label = region_name)), 
       filename = "output/region_map.png")


state_map <- ggplot(data = shp_state) +
  geom_sf(lwd = .1, fill = "#cbc0d3") + 
  ggrepel::geom_label_repel(
    aes(label = abbrev_state, geometry = geom),
    stat = "sf_coordinates") +
  theme_void() 

ggsave(state_map, filename = "output/state_map.png")


#### Climate plots ####
## State monthly tmean geofacet  (Figure S2)
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
  # using name intead of State
  facet_geo(~name, grid = grid_br) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15))

ggsave(tmean_heat, filename = "output/tmean_heat.png",
       height = 20, width = 20)


## Number of months suitable for transmission (Figure S3)
df_suitable <- df_year %>%
  mutate(decade = factor(ifelse(year %in% 2001:2010, "2001 - 2010", 
                                ifelse(year %in% 2011:2020, "2011 - 2020", NA)))) %>%
  group_by(municip_code_ibge, decade) %>%
  summarise(months_suitable = mean(months_suitable.both)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


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
df_suitable_diff <- spread(st_drop_geometry(df_suitable), decade, months_suitable) %>%
  left_join(., shp_parent, by = "municip_code_ibge") 

names(df_suitable_diff)[names(df_suitable_diff)=="2001 - 2010"] <- "dec1"
names(df_suitable_diff)[names(df_suitable_diff)=="2011 - 2020"] <- "dec2"

df_suitable_diff <- mutate(df_suitable_diff,
                           suitable_diff = dec2 - dec1,
                           # For plot (no difference = white)
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


## State monthly scPDSI geofacet  (Figure M1)
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
  # using name intead of State
  facet_geo(~name, grid = grid_br) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15))

ggsave(pdsi_heat, filename = "output/pdsi_heat.png",
       height = 20, width = 20)


## Number of months wet (scPDSI > 4) 
df_wet <- df_year %>%
  mutate(decade = factor(ifelse(year %in% 2001:2010, "2001 - 2010", 
                                ifelse(year %in% 2011:2020, "2011 - 2020", NA)))) %>%
  group_by(municip_code_ibge, decade) %>%
  summarise(months_wet = mean(months_wet)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()



## Difference in months extremely wet (Figure M2)
df_wet_diff <- spread(st_drop_geometry(df_wet), decade, months_wet) %>%
  left_join(., shp_parent, by = "municip_code_ibge") 

names(df_wet_diff)[names(df_wet_diff)=="2001 - 2010"] <- "dec1"
names(df_wet_diff)[names(df_wet_diff)=="2011 - 2020"] <- "dec2"

df_wet_diff <- mutate(df_wet_diff,
                      wet_diff = dec2 - dec1,
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
df_census <- df_year[df_year$year == 2010,] %>%
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


# 2000 census
urban00_map <- ggplot(data = df_census) +
  geom_sf(aes(fill = urban00), lwd = 0) +
  scale_fill_continuous_tableau(name = "% urbanisation", palette = "Blue-Teal") +
  theme_void()


# 2010 census
urban10_map <- ggplot(data = df_census) +
  geom_sf(aes(fill = urban10), lwd = 0) +
  scale_fill_continuous_tableau(name = "% urbanisation", palette = "Blue-Teal") +
  theme_void()


urban_maps <- plot_grid(urban00_map + theme(legend.position = "none"),
                        urban10_map + theme(legend.position = "none"),
                        labels = c("2000 census", "2010 census"))

urban_leg <- get_legend(urban00_map)

urban_maps <- plot_grid(urban_maps, urban_leg, rel_widths = c(3, .4))

ggsave(urban_maps, filename = "output/urban_decade.png", 
       height = 5, width = 10)


## Urban vs basic services scatterplots (Figure S5)
# Urbanisation vs. piped water
urban_water <- ggplot(data = df_census) +
  geom_point(aes(x = urban10, y = water_network, colour = region_name)) + 
  # geom_smooth(aes(x = urban10, y = water_network)) +
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

cor.test(df_census$urban10, df_census$water_network) # r = 0.656 [0.641, 0.671] p. < 0.001


# Urbanisation vs. refuse collection
urban_refuse <- ggplot(data = df_census) +
  geom_point(aes(x = urban10, y = total_collected, colour = region_name)) + 
  # geom_smooth(aes(x = urban10, y = total_collected)) +
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

cor.test(df_census$urban10, df_census$total_collected) # r = 0.794 [0.784, 0.804] p. < 0.001


#### REGIC plots ####
## REGIC 18 map (Figure 3)
regic18_map <- ggplot(data = df_census) +
  geom_sf(aes(fill = regic18), lwd = 0) +
  scale_fill_manual(name = "Level of \ninfluence", values = connect_cols) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
                      legend.text = element_text(size = 10))


ggsave(regic18_map, filename = "output/regic18_map.png")


## Histograms of % municipalities in each level per region (Figure S6)
df_regic07_perc <- st_drop_geometry(df_census) %>%
  group_by(region_name, regic07) %>%
  summarise(regic07_n = n()) %>%
  drop_na() %>%
  mutate(regic07_perc = regic07_n/sum(regic07_n) * 100)


regic07_hist <- ggplot(data = df_regic07_perc) +
  geom_bar(aes(x = region_name, y = regic07_perc, fill = regic07), stat = "identity") +
  scale_y_continuous(name = "% of municipalities", expand = c(0, 0)) +
  scale_x_discrete(name = "Region", expand = c(0, 0)) +
  scale_fill_manual(name = "Level of influence",  values = connect_cols) +
  theme_light()


df_regic18_perc <- st_drop_geometry(df_census) %>%
  group_by(region_name, regic18) %>%
  summarise(regic18_n = n()) %>%
  drop_na() %>%
  mutate(regic18_perc = regic18_n/sum(regic18_n) * 100)


regic18_hist <- ggplot(data = df_regic18_perc) +
  geom_bar(aes(x = region_name, y = regic18_perc, fill = regic18), stat = "identity") +
  scale_y_continuous(name = "% of municipalities", expand = c(0, 0)) +
  scale_x_discrete(name = "Region", expand = c(0, 0)) +
  scale_fill_manual(name = "Level of influence",  values = connect_cols) +
  theme_light()


regic_hist <- plot_grid(regic07_hist + theme(legend.position = "none"),
                        regic18_hist + theme(legend.position = "none"))

regic_leg <- get_legend(regic07_hist)

regic_hist <- plot_grid(regic_hist, regic_leg, rel_widths = c(3, .4))

ggsave(regic_hist, filename = "output/regic_hist.png",
       height = 5, width = 15)




## Raincloud plots REGIC vs census variables (Figure S7)
# Urbanisation
regic_urb <- ggplot(df_census, aes(x = regic18, y = urban10, fill = regic18)) +
  stat_halfeye(adjust = .5,
               .width = 0,
               justification = -.2, scale = .6)+ 
  geom_boxplot(width = .15, outlier.shape = NA, position = "dodge") +
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


# Piped water
regic_water <- ggplot(df_census, aes(x = regic18, y = water_network, fill = regic18)) +
  stat_halfeye(adjust = .5,
               .width = 0,
               justification = -.2,
               scale = .6) +
  geom_boxplot(width = .15, outlier.shape = NA) +
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


regic_refuse <- ggplot(df_census, aes(x = regic18, y = total_collected, fill = regic18)) +
  stat_halfeye(adjust = .5,
               .width = 0,
               justification = -.2,
               scale = .6) +
  geom_boxplot(width = .15, outlier.shape = NA) +
  geom_half_point(side = "l", range_scale = .4,
                  alpha = .3) +
  labs(x = "Level of influence", y = "% with refusse collection") +
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
## Basic summaries of each level (Table S1)
df_regic07 <- df_year[df_year$year == 2007,] %>%
  mutate(regic07 = factor(level07_acpnum, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre")))


df_regic18 <- df_year[df_year$year == 2018,] %>%
  mutate(regic18 = factor(level18_num, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre")))


# By region
addmargins(table(df_regic07$region_name, df_regic07$regic07))


addmargins(table(df_regic18$region_name, df_regic18$regic18))



#### Dengue case exploratory plots ####
## Regional DIR timeplot (Figure 4)
df_region_dir <- df_month %>%
  group_by(region_name, time) %>%
  summarise(dengue_total = sum(dengue_cases), 
            population_total = sum(population)) %>%
  mutate(DIR = (dengue_total/population_total) * 10^5) %>%
  ungroup() 


## DIR time series per month
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


## Combine timeseries with inset region map
# Use map legend
region_leg <- get_legend(region_dir_time)

region_line_inset <- ggdraw() +
  draw_plot(region_dir_time + theme(legend.position = "none")) +
  draw_plot(region_map + theme(legend.position = "none"), 
            x = .01, y = .5, width = .5, height = .5) 


region_line_comp <- plot_grid(region_line_inset, region_leg,
                              rel_widths = c(3, .4))

ggsave(region_line_comp, filename = "output/region_dir_inset.png",
       height = 5, width = 10)


#### Outbreak plots ####
## Number of years with outbreak (Figure 1)
df_permanence <- df_year %>%
  mutate(DIR_year = (dengue_year/population)*10^5,
         outbreak = ifelse(DIR_year >= 300, 1, 0)) %>%
  group_by(municip_code_ibge) %>%
  summarise(number_outbreaks = sum(outbreak)) %>%
  ungroup() %>%
  # To make 0 outbreaks grey
  mutate(number_outbreaks_na = ifelse(number_outbreaks == 0, NA, 
                                      number_outbreaks)) %>%
  full_join(., shp_parent, by = "municip_code_ibge") %>%
  filter(municip_code_ibge !=2605459) %>%
  st_as_sf()


number_outbreaks <- ggplot(data = df_permanence) +
  geom_sf(aes(fill = number_outbreaks_na, colour = ""), lwd = 0) +
  scale_fill_gradient_tableau("Purple", name = "Number of years \nwith outbreaks", 
                              na.value = "grey") +
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No outbreaks")) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


ggsave(number_outbreaks, filename = "output/number_outbreaks.png",
       height = 10, width = 10)


## Percentage of municipalities per year experiencing outbreaks (Figure S8)
municip_region <- st_drop_geometry(df_census) %>%
  group_by(region_name) %>%
  summarise(n_municip = n()) 


df_region_outbreaks <- df_year %>%
  mutate(DIR_year = (dengue_year/population)*10^5,
         outbreak = ifelse(DIR_year >= 300, 1, 0)) %>%
  group_by(region_name, year) %>%
  summarise(total_outbreak = sum(outbreak)) %>%
  ungroup() %>%
  full_join(., municip_region, by = "region_name") %>%
  mutate(perc_outbreak = (total_outbreak/n_municip) * 100)


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


perc_region_inset <- ggdraw() +
  draw_plot(perc_region_outbreak + theme(legend.position = "none")) +
  draw_plot(region_map + theme(legend.position = "none"), 
            x = .01, y = .5, width = .5, height = .5) 


perc_region_inset <- plot_grid(perc_region_inset, region_leg,
                               rel_widths = c(3, .4))

ggsave(perc_region_inset, filename = "output/region_perc_inset.png",
       height = 5, width = 10)


## Year of first outbreak (Figure 5)
df_outbreak <- df_year %>%
  mutate(DIR = (dengue_year/population)*10^5) %>%
  filter(DIR >= 300) %>%
  arrange(municip_code_ibge, year) %>%
  group_by(municip_code_ibge) %>%
  mutate(first_fix = first(year)) %>%
  filter(year == first_fix) %>%
  ungroup()

df_outbreak <- df_outbreak[!duplicated(df_outbreak$municip_code_ibge), 
                           c("municip_code_ibge", "first_fix")] %>%
  full_join(., shp_parent, by = "municip_code_ibge") %>%
  filter(municip_code_ibge !=2605459) %>%
  mutate(first_outbreak10 = ifelse(first_fix > 2009, NA, first_fix)) %>%
  st_as_sf()


# First outbreak 2001 - 2020 
first_outbreak20 <- ggplot(data = df_outbreak) + 
  geom_sf(aes(fill = first_fix, colour = ""), lwd = .05) +
  scale_fill_gradient_tableau("Purple", na.value = "grey",
                              name = "Year of \nfirst outbreak",
                              trans = "reverse") +
  expand_limits(fill = c(2001, 2020)) +
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No outbreaks")) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


# First outbreak 2001 - 2010
first_outbreak10 <- ggplot(data = df_outbreak) + 
  geom_sf(aes(fill = first_outbreak10, colour = ""), lwd = .05) +
  scale_fill_gradient_tableau("Purple", na.value = "grey",
                              name = "Year of \nfirst outbreak",
                              trans = "reverse") +
  expand_limits(fill = c(2001, 2020)) +
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No outbreaks")) +
  theme_void() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))

outbreak_leg <- get_legend(first_outbreak20)


first_outbreak_decade <- plot_grid(first_outbreak10 + theme(legend.position = "none"),
                                   first_outbreak20 + theme(legend.position = "none"),
                                   labels = c("2001 - 2010", "2001 - 2020"))

first_outbreak_decade <- plot_grid(first_outbreak_decade, outbreak_leg, 
                                   rel_widths = c(3, .4))


ggsave(first_outbreak_decade, filename = "output/first_outbreak_decade.png",
       height = 5, width = 10)


#### Plot 75th percentile outbreak threshold (Figure M) ####
df_cutoff <- df_year %>% 
  group_by(municip_code_ibge) %>% 
  summarise(perc_75 = quantile(DIR, .75),
            pop_mean = mean(population)) %>% 
  ungroup() %>% 
  # Add minimum cutoff (5 cases) 
  mutate(DIR_5 = (5/pop_mean)*10^5,
         perc75_cutoff = ifelse(perc_75 < DIR_5, DIR_5, perc_75)) %>% 
  left_join(., shp_parent, by = "municip_code_ibge") %>% 
  st_as_sf()

perc75_cutoff <- ggplot(data = df_cutoff) +
  geom_sf(aes(fill = perc75_cutoff), lwd = .05) +
  scale_fill_gradient2_tableau("Gold-Purple Diverging", 
                               name = "Outbreak threshold",
                               trans = "log1p",
                               breaks = c(0, 30, 100, 300, 1000, 3000)) +
  # expand_limits(fill = c(0, 10000))  +
  theme_void() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

ggsave(perc75_cutoff, filename = "output/perc75_cutoff.png")
