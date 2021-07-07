#################################################
####                                         ####
####     Script to explore model results     ####
####                                         ####
#################################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, sf, mgcv, mgcViz, cowplot)


#### Load data ####
## Model data
df_year <- fread("data/dengue_year.csv")


df_model <- df_year %>%
  mutate(DIR_year = (dengue_year/population)*10^5,
         outbreak = ifelse(DIR_year >= 300, 1, 0),
         outbreak100 = ifelse(DIR_year >= 100, 1, 0),
         # Rearrange REGIC to set local centre as reference
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
         regic_comb = factor(ifelse(year %in% 2001:2009, regic07_relevel, regic18_relevel),
                             levels = 1:5,
                             labels = c("Local centre",
                                        "Zone centre",
                                        "Sub-regional centre",
                                        "Regional capital",
                                        "Metropolis")),
         # Use urban proportion, 00 for 2001 - 2009 + 10 for 2010 - 2020 (unless mising in 2000)
         urban_prpn = ifelse(!is.na(urban00), 
                             ifelse(year %in% 2001:2009, (urban00/100), 
                                    (urban10/100)), (urban10/100))) %>%
  dplyr::select(municip_code_ibge, year, lon, lat, outbreak, regic_comb, 
                urban_prpn, months_suitable.era, prior_outbreak)


## Shapefile
# Load raw municipality data
shp <- read_municipality()

# Convert into parent municipalities
parent_conv <- fread("data/parent_municip.csv")

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
  mutate(municip_code = as.numeric(substr(municip_code_ibge, 1, 6)))


#### Load model objects ####
model_base <- readRDS("output/model_base.rds")
model_full <- readRDS("output/model_full.rds")

model_full100 <- readRDS("output/model_full100.rds")


set.seed(123)


#### Investigate fixed model parameters ####
## Simulate from posterior distribution of parameters
# Function to obtain MVN random variates
rmvn <- function(n, mu, sig) { 
  L <- mroot(sig); m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}


# Set number of simulations
n.sims <- 1000


# Draw random samples from beta posterior distribution
betas_base <- rmvn(n.sims, coef(model_base), model_base$Vp)
betas_full <- rmvn(n.sims, coef(model_full), model_full$Vp) 

betas_full100 <- rmvn(n.sims, coef(model_full100), model_full100$Vp) 



## Obtain estimate + 95% credible interval for each 
# Full model
beta_ci <- data.frame(betas_full)
names(beta_ci) <- gsub("[[:punct:][:blank:]]+","", names(model_full$coefficients))


beta_ci <- beta_ci %>%
  summarise(across(2:8, list(mean = mean, lower = ~quantile(.x, 0.025), 
                             upper = ~quantile(.x, 0.975)),
                   names = "{.col}.fn{.fn}")) %>%
  # Convert onto aOR scale
  mutate(across(1:21, ~exp(.x)))


# Convert into form that can be plotted
beta_mean <- beta_ci %>%
  pivot_longer(cols = ends_with("_mean"),
               names_to = "Covariate",
               values_to = "Mean") %>%
  transmute(Covariate = gsub("\\_.*", "", Covariate),
            Mean = Mean)

beta_lci <- beta_ci %>%
  pivot_longer(cols = ends_with("_lower"),
               names_to = "Covariate",
               values_to = "LCI") %>%
  transmute(Covariate = gsub("\\_.*", "", Covariate),
            LCI = LCI)

beta_uci <- beta_ci %>%
  pivot_longer(cols = ends_with("upper"),
               names_to = "Covariate",
               values_to = "UCI") %>%
  transmute(Covariate = gsub("\\_.*", "", Covariate),
            UCI = UCI)


## Mean & 95% CI of parameter estimates (Table 1)
beta_estimates <- full_join(beta_mean, beta_lci, by = "Covariate") %>%
  full_join(., beta_uci, by = "Covariate") 


# Plot mean & 95% CI (Figure 7)
beta_ci_plot <- ggplot(data = beta_estimates, aes(Covariate)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
  geom_point(aes(x = Covariate, y = Mean, colour = Covariate), shape = 16, 
             size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Covariate), lwd = 1,
                 position = position_dodge(width = 0.5)) +
  scale_x_discrete(name = "Model covariate", 
                   labels = c("Months suitable", "Prior outbreak: Yes", 
                              "Metropolis", "Regional capital",
                              "Sub-regional centre", "Zone centre", 
                              "Urbanisation")) + 
  scale_colour_manual(values = c("#ff0054", "#390099", rep("#9e0059", 4), "#197278")) +
  expand_limits(y = c(1, 4)) +
  coord_flip() + 
  scale_y_continuous(name = "Coefficient estimate (aOR)",
                     position = "right") +
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = "None")  

ggsave(beta_ci_plot, filename = "output/beta_ci.png")


## Sensitivity analysis - set DIR > 100 as outbreak
beta_ci100 <- data.frame(betas_full100)
names(beta_ci100) <- gsub("[[:punct:][:blank:]]+","", names(model_full100$coefficients))


beta_ci100 <- beta_ci100 %>%
  summarise(across(2:8, list(mean = mean, lower = ~quantile(.x, 0.025), upper = ~quantile(.x, 0.975)),
                   names = "{.col}.fn{.fn}")) %>%
  # Convert onto aOR scale
  mutate(across(1:21, ~exp(.x)))


# Convert into form that can be plotted
beta_mean100 <- beta_ci100 %>%
  pivot_longer(cols = ends_with("_mean"),
               names_to = "Covariate",
               values_to = "Mean") %>%
  transmute(Covariate = gsub("\\_.*", "", Covariate),
            Mean = Mean)

beta_lci100 <- beta_ci100 %>%
  pivot_longer(cols = ends_with("_lower"),
               names_to = "Covariate",
               values_to = "LCI") %>%
  transmute(Covariate = gsub("\\_.*", "", Covariate),
            LCI = LCI)

beta_uci100 <- beta_ci100 %>%
  pivot_longer(cols = ends_with("upper"),
               names_to = "Covariate",
               values_to = "UCI") %>%
  transmute(Covariate = gsub("\\_.*", "", Covariate),
            UCI = UCI)


## Mean & 95% CI of parameter estimates (Table S1)
beta_estimates100 <- full_join(beta_mean100, beta_lci100, by = "Covariate") %>%
  full_join(., beta_uci100, by = "Covariate") 


# Plot estimates + 95% CI (Figure S7)
beta_ci_plot100 <- ggplot(data = beta_estimates100, aes(Covariate)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
  geom_point(aes(x = Covariate, y = Mean, colour = Covariate), shape = 16, 
             size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Covariate), lwd = 1,
                 position = position_dodge(width = 0.5)) +
  scale_x_discrete(name = "Model covariate", 
                   labels = c("Months suitable", "Prior outbreak: Yes", 
                              "Metropolis", "Regional capital",
                              "Sub-regional centre", "Zone centre", 
                              "Urbanisation")) + 
  scale_colour_manual(values = c("#ff0054", "#390099", rep("#9e0059", 4), "#197278")) +
  expand_limits(y = c(1, 4)) +
  coord_flip() + 
  scale_y_continuous(name = "Coefficient estimate (aOR)",
                     position = "right") +
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = "None")  

ggsave(beta_ci_plot100, filename = "output/beta_ci100.png")



#### Investigate smooth functions ####
## Use mgcViz to plot temporal smooth function 
gam_full_viz <-getViz(model_full)


## Plot temporal smooth function (Figure 8a)
temp_smooth_plot.full <- plot(sm(gam_full_viz, 2), trans = exp) +
  labs(x = "Year", y = "s(Year)") +
  l_fitLine( ) +
  l_ciLine() + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))


## Return matrix of linear predictors (used to generate simulations)
model_mat.full <- predict(object = model_full, newdata = df_model, 
                          type = "lpmatrix")


# Select only spatial parameters
spat_mat.full <- model_mat.full[ ,substr(names(coefficients(model_full)), 1, 5) == "s(lon"]
betas_spat.full<- betas_full[ , substr(names(coefficients(model_full)), 1, 5) == "s(lon"]


# Estimate spatial random effects + transform from log scale
spat_est.full <- exp(spat_mat.full %*% t(betas_spat.full))
# As all years are equal, just keep 2001
spat_est01.full <- spat_est.full[df_model$year == 2001, ] 

# Take average of simulations
mean_spat.full <- apply(spat_est01.full, 1, mean)


# Join random effects to dataset 
spat_re.full <- data.table(municip_code_ibge = df_model[df_model$year == 2001,]$municip_code_ibge,
                           spat_re = mean_spat.full,
                           # Make areas with no change in risk from expected white
                           spat_re_na = ifelse(mean_spat.full == 1, NA, mean_spat.full)) %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


## Map with spatial random field (Figure 8b)
spat_re_plot <- ggplot(data = spat_re.full) +
  geom_sf(aes(fill = spat_re), lwd = 0) +
  scale_fill_gradient2(name = "Spatial random\nfield", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 1, 
                       trans = "log1p", breaks = c(0, 1, 5, 10, 15)) +
  theme_void()

ggsave(spat_re_plot, filename = "output/spatial_re_map_full.png",
       height = 10, width = 10)


#### Obtain predictions of the probability of an outbreak ####
## Use linear prediction matrix and parameter simulations to estimate posterior mean
mean_prob <- exp(model_mat.full %*% t(betas_full))/(1 + exp(model_mat.full %*% t(betas_full)))
  
## Simulate from the posterior distribution
preds <- apply(mean_prob, 1, function(x){rbinom(length(x), size = 1, prob = x)})
  
## Obtain mean probability
mean_pred <- apply(preds, 2, mean)
  
## Convert into a dataset to plot  
predictions <- data.table(municip_code_ibge = df_model$municip_code_ibge,
                          year = df_model$year,
                          outbreak_obs = df_model$outbreak,
                          prob_pred = mean_pred)  %>%
  full_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


## Plot predicted probability per year (Figure S8)
pred_year_map <- ggplot(data = predictions) +
  geom_sf(aes(fill = prob_pred), lwd = 0) +
  facet_wrap(~ year) +
  scale_fill_viridis_c(name = "Probability of\nan outbreak", direction = 1) +
  coord_sf(datum = NA) +
  expand_limits(fill = c(0, 1)) +
  theme_void()


ggsave(pred_year_map, filename = "output/pred_year_maps.png")


## Plot average probability 2001 - 2010 vs 2011 - 2020 
predictions_decade <- st_drop_geometry(predictions) %>%
  mutate(decade = factor(ifelse(year %in% 2001:2010, 0, 
                                ifelse(year %in% 2011:2020, 1, NA)),
                         levels = c(0, 1),
                         labels = c("2001 - 2010", "2011 - 2020"))) %>%
  group_by(municip_code_ibge, decade) %>%
  summarise(mean_prob = mean(prob_pred, na.rm = T)) %>%
  ungroup() %>%
  # Add different risk levels to identify barriers
  mutate(no_risk = factor(ifelse(mean_prob == 0, 0, 1), levels = c(0, 1), 
                          labels = c("0%", "> 0%")),
         risk_5 = factor(ifelse(mean_prob <= .05, 0, 1), levels = c(0, 1),
                         labels = c("\u2264 5%", "> 5%")),
         risk_10 = factor(ifelse(mean_prob <= .1, 0, 1), levels = c(0, 1),
                          labels = c("\u2264 10%", "> 10%")),
         risk_15 = factor(ifelse(mean_prob <= .15, 0, 1), levels = c(0, 1),
                          labels = c("\u2264 15%", "> 15%")),
         risk_20 = factor(ifelse(mean_prob <= .2, 0, 1), levels = c(0, 1),
                          labels = c("\u2264 20%", "> 20%")),
         risk_25 = factor(ifelse(mean_prob <= .25, 0, 1), levels = c(0, 1),
                          labels = c("\u2264 25%", "> 25%"))) %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


## Plot mean probability per decade, Brazil (Figure 9a)
predictions_decade_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = mean_prob), col = NA) +
  facet_wrap(~ decade) +
  scale_fill_viridis_c(name = "Probability of\noutbreak", direction = 1) +
  coord_sf(datum = NA) +
  expand_limits(fill = c(0, 0.8)) +
  theme_void()

ggsave(predictions_decade_map, filename = "output/predictions_decade_map.png",
       height = 10, width = 15)

## Plot mean probability per decade, South Brazil (Figure 9b)
predictions_decade_sth <- ggplot(data = predictions_decade[substr(predictions_decade$municip_code_ibge, 1, 1) == "4",]) +
  geom_sf(aes(fill = mean_prob), col = NA) +
  facet_wrap(~ decade) +
  scale_fill_viridis_c(name = "Probability of\noutbreak", direction = 1) +
  coord_sf(datum = NA) +
  expand_limits(fill = c(0, 0.8)) +
  theme_void()

ggsave(predictions_decade_sth, filename = "output/predictions_decade_sth.png",
       height = 10, width = 15)

## Plot mean probability per decade, Acre + Amazonas (Figure 9c)
predictions_decade_am <- ggplot(data = predictions_decade[substr(predictions_decade$municip_code_ibge, 1, 2) %in% c("13","12"),]) +
  geom_sf(aes(fill = mean_prob), col = NA) +
  facet_wrap(~ decade) +
  scale_fill_viridis_c(name = "Probability of\noutbreak", direction = 1) +
  coord_sf(datum = NA) +
  expand_limits(fill = c(0, 0.8)) +
  theme_void()

ggsave(predictions_decade_am, filename = "output/predictions_decade_am.png",
       height = 10, width = 15)


#### Maps showing barriers to dengue transmission (where P(outbreak) < threshold) ####
## Plot showing mean(P(outbreak)) = 0%
no_risk_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = no_risk), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "Status", labels = c("Protected", "Not protected"), 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()


## Plot showing mean(P(outbreak)) = 5%
risk_5_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = risk_5), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "P(outbreak)", 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()


## Plot showing mean(P(outbreak)) = 10% (Figure 10)
risk_10_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = risk_10), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "P(outbreak)", 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()

ggsave(risk_10_map, filename = "output/risk_10_map.png",
       height = 10, width = 10)

## Plot showing mean(P(outbreak)) = 15%
risk_15_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = risk_15), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "P(outbreak)", 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()


## Plot different threshold barriers (Figure S9)
barrier_maps <- plot_grid(no_risk_map + theme(legend.position = "none"), 
                          risk_5_map + theme(legend.position = "none"), 
                          risk_10_map + theme(legend.position = "none"),
                          risk_15_map + theme(legend.position = "none"),
                          labels = c("Threshold: 0%", "Threshold: 5%", 
                                     "Threshold: 10%", "Threshold: 15%"),
                          nrow = 2
                          )

protect_leg <- get_legend(no_risk_map)

barrier_maps <- plot_grid(barrier_maps, protect_leg, rel_widths = c(4, .3))

ggsave(barrier_maps, filename = "output/barrier_comp.png",
       height = 10, width = 20)
