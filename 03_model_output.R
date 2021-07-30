#################################################
####                                         ####
####     Script to explore model results     ####
####                                         ####
#################################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, sf, geobr, mgcv, mgcViz, cowplot)

# Run models
# source("02_model_fitting.R")

#### Load data ####
## Model data
df_year <- fread("data/df_model.csv")


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
                prior_outbreak, outbreak_prevyr) %>%
  drop_na()


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
  mutate(municip_code = as.numeric(substr(municip_code_ibge, 1, 6)))


#### Load model objects ####
model_base <- readRDS("output/model_base.rds")
model_full <- readRDS("output/model_full.rds")

## Sensitivity analysis - different variable options
model_full100 <- readRDS("output/model_full100.rds")
model_outbreak <- readRDS("output/model_outbreak.rds")
model_aeg <- readRDS("output/model_aeg.rds")

## Sensistivity analysis - remove each variable in turn
model_clim <- readRDS("output/model_clim.rds")
model_prior <- readRDS("output/model_prior.rds")
model_urb <- readRDS("output/model_urb.rds")
model_regic <- readRDS("output/model_regic.rds")


set.seed(123)


#### Bayesian estimation of the beta parameters ####
# Function to obtain MVN random variates
rmvn <- function(n, mu, sig) { 
  L <- mroot(sig); m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}


# Number of simulations
n.sims <- 1000

betas_full <- rmvn(n.sims, coef(model_full), model_full$Vp)
betas_full100 <- rmvn(n.sims, coef(model_full100), model_full100$Vp)
betas_aeg <- rmvn(n.sims, coef(model_aeg), model_aeg$Vp)
betas_outbreak <- rmvn(n.sims, coef(model_outbreak), model_outbreak$Vp)


## Plot betas + 95% CI
beta_ci <- function(model, betas) {
  beta_ci <- as_tibble(betas) 
  names(beta_ci) <- gsub("[[:punct:][:blank:]]+","", names(model$coefficients))
  
  beta_ci <- beta_ci %>%
    summarise(across(2:8, list(mean = mean, lower = ~quantile(.x, 0.025), upper = ~quantile(.x, 0.975)),
                     names = "{.col}.fn{.fn}")) %>%
    mutate(across(everything(), ~exp(.x)))
  
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
  
  
  beta_estimates <- full_join(beta_mean, beta_lci, by = "Covariate") %>%
    full_join(., beta_uci, by = "Covariate") 
  
  
  return(beta_estimates)
}

# Full model (Figure 6)
beta_ci_full <- beta_ci(model_full, betas_full)

beta_plot_full <- ggplot(data = beta_ci_full, aes(Covariate)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
  geom_point(aes(x = Covariate, y = Mean, colour = Covariate), shape = 16, 
             size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Covariate), lwd = 1,
                 position = position_dodge(width = 0.5)) +
  scale_x_discrete(name = "Model covariate",
                   labels = c("Months suitable", "Prior outbreak: Yes", "Metropolis", "Regional capital",
                              "Sub-regional centre", "Zone centre", "Urbanisation")) +
  scale_colour_manual(values = c("#ff0054", "#390099", rep("#9e0059", 4), "#197278")) +
  coord_flip() + 
  scale_y_continuous(name = "Coefficient estimate (aOR)",
                     position = "right",
                     expand = c(0,0), limits = c(1, 4)) +
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = "None")  

ggsave(beta_plot_full, filename = "output/beta_ci_full.png")


# Outbreak > 100
beta_ci_full100 <- beta_ci(model_full100, betas_full100)

# beta_plot_full100 <- ggplot(data = beta_ci_full100, aes(Covariate)) +
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
#   geom_point(aes(x = Covariate, y = Mean, colour = Covariate), shape = 16, 
#              size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Covariate), lwd = 1,
#                  position = position_dodge(width = 0.5)) +
#   scale_x_discrete(name = "Model covariate",
#                    labels = c("Months suitable", "Prior outbreak: Yes", "Metropolis", "Regional capital",
#                               "Sub-regional centre", "Zone centre", "Urbanisation")) +
#   scale_colour_manual(values = c("#ff0054", "#390099", rep("#9e0059", 4), "#197278")) +
#   coord_flip() + 
#   scale_y_continuous(name = "Coefficient estimate (aOR)",
#                      position = "right",
#                      expand = c(0,0), limits = c(1, 4)) +
#   theme_classic() +
#   theme(plot.margin = unit(c(1,1,1,1),"cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15),
#         legend.position = "None")  
# 
# ggsave(beta_plot_full100, filename = "output/beta_ci_full100.png")


# Using aegypti boundary
beta_ci_aeg <- beta_ci(model_aeg, betas_aeg)

# beta_plot_aeg <- ggplot(data = beta_ci_aeg, aes(Covariate)) +
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
#   geom_point(aes(x = Covariate, y = Mean, colour = Covariate), shape = 16, 
#              size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Covariate), lwd = 1,
#                  position = position_dodge(width = 0.5)) +
#   scale_x_discrete(name = "Model covariate",
#                    labels = c("Months suitable", "Prior outbreak: Yes", "Metropolis", "Regional capital",
#                               "Sub-regional centre", "Zone centre", "Urbanisation")) +
#   scale_colour_manual(values = c("#ff0054", "#390099", rep("#9e0059", 4), "#197278")) +
#   coord_flip() + 
#   scale_y_continuous(name = "Coefficient estimate (aOR)",
#                      position = "right",
#                      expand = c(0,0), limits = c(1, 4)) +
#   theme_classic() +
#   theme(plot.margin = unit(c(1,1,1,1),"cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15),
#         legend.position = "None")  
# 
# ggsave(beta_plot_aeg, filename = "output/beta_ci_aeg.png")


# Using prior outbreak form previous year
beta_ci_outbreak <- beta_ci(model_outbreak, betas_outbreak)

# beta_plot_outbreak <- ggplot(data = beta_ci_outbreak, aes(Covariate)) +
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
#   geom_point(aes(x = Covariate, y = Mean, colour = Covariate), shape = 16, 
#              size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Covariate), lwd = 1,
#                  position = position_dodge(width = 0.5)) +
#   scale_x_discrete(name = "Model covariate",
#                    labels = c("Months suitable", "Prior outbreak: Yes", "Metropolis", "Regional capital",
#                               "Sub-regional centre", "Zone centre", "Urbanisation")) +
#   scale_colour_manual(values = c("#ff0054", "#390099", rep("#9e0059", 4), "#197278")) +
#   coord_flip() + 
#   scale_y_continuous(name = "Coefficient estimate (aOR)",
#                      position = "right",
#                      expand = c(0,0), limits = c(1, 4.5)) +
#   theme_classic() +
#   theme(plot.margin = unit(c(1,1,1,1),"cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15),
#         legend.position = "None")  
# 
# ggsave(beta_plot_outbreak, filename = "output/beta_ci_outbreak.png")



#### Plot smooth functions ####
## Temporal smooth using mgcViz (Figure 7a)
model_full_viz <- getViz(model_full)

temp_plot.full <- plot(sm(model_full_viz, 2), trans = exp) +
  labs(x = "Year", y = "s(Year)") +
  l_fitLine( ) +
  l_ciLine() +
  geom_hline(yintercept = 1, linetype = "dashed") 


## Spatial smooth using simulations from the posterior
model_mat.full <- predict(object = model_full, newdata = df_model, type = "lpmatrix")

# Select only spatial parameters
spat_mat.full <- model_mat.full[ ,substr(names(coefficients(model_full)), 1, 5) == "s(lon"]
betas_spat.full<- betas_full[ , substr(names(coefficients(model_full)), 1, 5) == "s(lon"]

# Estimate spatial random effects + transform from log scale
spat_est.full <- exp(spat_mat.full %*% t(betas_spat.full))
# As all years are equal, just keep 2001
spat_est01.full <- spat_est.full[df_model$year == 2001, ] 

# Take average of simulations
mean_spat.full <- apply(spat_est01.full, 1, mean)

# Create 95%. credible interval
lci_spat.full <- apply(spat_est01.full, 1, quantile, probs=0.025)
uci_spat.full <- apply(spat_est01.full, 1, quantile, probs=0.975)


# Join random effects to dataset + plot on a map
spat_smooth.full <- data.table(municip_code_ibge = df_model[df_model$year == 2001,]$municip_code_ibge,
                               spat_smooth = mean_spat.full,
                               spat_re_uci = ifelse(uci_spat.full == 1, NA, uci_spat.full),
                               spat_re_lci =  ifelse(lci_spat.full == 1, NA, lci_spat.full)) %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


# Plot spatial smooth field (Figure 7b)
spat_smooth_plot <- ggplot(data = spat_smooth.full) +
  geom_sf(aes(fill = spat_smooth), lwd = 0) +
  scale_fill_gradient2(name = "Spatial random\nfield", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 1, 
                       trans = "log1p", breaks = c(0, 1, 5, 10, 15)) +
  theme_void()

ggsave(spat_smooth_plot, filename = "output/spat_smooth_plot.png",
       height = 10, width = 10)


#### Compare smooth terms between models ####
## Function to obtain smooth functions
smooth_estimates <- function(model, data) {
  # Obtain simulations from beta posterior distributions (only keep smooth functions)
  betas <- rmvn(1000, coef(model), model$Vp) 
  betas_smooth <- betas[ , substr(names(coefficients(model)), 1, 2) %in%
                           c("s(", "ti")]
  
  # Obtain linear predictor matrix for estimates (only keep smooth functions)
  model_mat <- predict(object = model, newdata = data, type = "lpmatrix")
  model_mat_smooth <- model_mat[ , substr(names(coefficients(model)), 1, 2) %in%
                                   c("s(", "ti")]
  
  # Estimate smooth function for each simulation
  smooth_est <- model_mat_smooth %*% t(betas_smooth)
  
  
  return(smooth_est)
  
}

## Baseline model smooths
smooth_base <- smooth_estimates(model_base, df_model)

## Final model smooths
smooth_final <- smooth_estimates(model_full, df_model)

base_full_diff <- abs(smooth_final) - abs(smooth_base)

base_full_diff.med <- apply(base_full_diff, 1, median)
base_full_diff.lq <- apply(base_full_diff, 1, quantile, .25)
base_full_diff.uq <- apply(base_full_diff, 1, quantile, .75)



## Compare baseline and final models, plot absolute differences (Figure 8)
# (if final model is better, would expect reduction in absolute values/shrinkage to 0)
# df_model$smooth_absdiff <- abs(df_model$smooth_final) - abs(df_model$smooth_base)
df_model$smooth_absdiff <- base_full_diff.med

## Join to map and plot differences per year
base_full_comp <- left_join(df_model, shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


## Plot each year inidividually
# for(i in min(base_full_comp$year):max(base_full_comp$year)) {
#   base_full_map <- ggplot(data = base_full_comp[base_full_comp$year == i,]) +
#     geom_sf(aes(fill = smooth_absdiff), lwd = .05) +
#     scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
#                          mid = "white", high = "#c51b7d", midpoint = 0) +
#     expand_limits(fill = c(-4, 4)) +
#     theme_void() +
#     facet_wrap(~year) 
#   
#   ggsave(base_full_map, filename = paste0("output/base_full_comp", i, ".png"))
# }


# Take average over the time period and plot
beta_full_comp_av <- group_by(df_model, municip_code_ibge) %>%
  summarise(diff_med = median(smooth_absdiff)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


base_full_med <- ggplot(data = beta_full_comp_av) +
  geom_sf(aes(fill = diff_med), lwd = .05) +
  scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void() 

ggsave(base_full_med, file = "output/base_full_med.png")


#### Compare smooth functions after removing covariates in turn ####
# Urbanisation
smooth_urb <- smooth_estimates(model_urb, df_model)
# df_model$urb_full.diff <- abs(smooth_final) - abs(smooth_urb)
urb_full_diff <- abs(smooth_final) - abs(smooth_urb)

urb_full_diff.med <- apply(urb_full_diff, 1, median)
urb_full_diff.lq <- apply(urb_full_diff, 1, quantile, .25)
urb_full_diff.uq <- apply(urb_full_diff, 1, quantile, .75)

df_model$urb_full.diff <- urb_full_diff.med

# Prior outbreak
smooth_prior <- smooth_estimates(model_prior, df_model)

prior_full_diff <- abs(smooth_final) - abs(smooth_prior)

prior_full_diff.med <- apply(prior_full_diff, 1, median)
prior_full_diff.lq <- apply(prior_full_diff, 1, quantile, .25)
prior_full_diff.uq <- apply(prior_full_diff, 1, quantile, .75)

df_model$prior_full.diff <- prior_full_diff.med

# Connectivity
smooth_regic <- smooth_estimates(model_regic, df_model)

regic_full_diff <- abs(smooth_final) - abs(smooth_regic)

regic_full_diff.med <- apply(regic_full_diff, 1, median)
regic_full_diff.lq <- apply(regic_full_diff, 1, quantile, .25)
regic_full_diff.uq <- apply(regic_full_diff, 1, quantile, .75)

df_model$regic_full.diff <- regic_full_diff.med

# Number of months suitable
smooth_clim <- smooth_estimates(model_clim, df_model)

clim_full_diff <- abs(smooth_final) - abs(smooth_clim)

clim_full_diff.med <- apply(clim_full_diff, 1, median)
clim_full_diff.lq <- apply(clim_full_diff, 1, quantile, .25)
clim_full_diff.uq <- apply(clim_full_diff, 1, quantile, .75)

df_model$clim_full.diff <- clim_full_diff.med


## Join the differences to map and plot differences per year
base_full_comp <- left_join(df_model, shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


# for(i in min(base_full_comp$year):max(base_full_comp$year)) {
#   base_full_map <- ggplot(data = base_full_comp[base_full_comp$year == i,]) +
#     geom_sf(aes(fill = urb_full.diff), lwd = .05) +
#     scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
#                          mid = "white", high = "#c51b7d", midpoint = 0) +
#     expand_limits(fill = c(-4, 4)) +
#     theme_void() +
#     facet_wrap(~year) 
#   
#   ggsave(base_full_map, 
#          filename = paste0("output/full_urb_comp", i, ".png"))
# }
# 
# 
# for(i in min(base_full_comp$year):max(base_full_comp$year)) {
#   base_full_map <- ggplot(data = base_full_comp[base_full_comp$year == i,]) +
#     geom_sf(aes(fill = clim_full.diff), lwd = .05) +
#     scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
#                          mid = "white", high = "#c51b7d", midpoint = 0) +
#     expand_limits(fill = c(-4, 4)) +
#     theme_void() +
#     facet_wrap(~year) 
#   
#   ggsave(base_full_map, 
#          filename = paste0("output/full_clim_comp", i, ".png"))
# }

## Plot 2020 difference for climate (Figure S9)
clim_diff_20 <- ggplot(data = base_full_comp[base_full_comp$year == 2020,]) +
      geom_sf(aes(fill = clim_full.diff), lwd = .05) +
      scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221",
                           mid = "white", high = "#c51b7d", midpoint = 0) +
      expand_limits(fill = c(-3, 3)) +
      theme_void()


ggsave(clim_diff_20, filename = "output/clim_diff20.png")
# 
# for(i in min(base_full_comp$year):max(base_full_comp$year)) {
#   base_full_map <- ggplot(data = base_full_comp[base_full_comp$year == i,]) +
#     geom_sf(aes(fill = regic_full.diff), lwd = .05) +
#     scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
#                          mid = "white", high = "#c51b7d", midpoint = 0) +
#     expand_limits(fill = c(-4, 4)) +
#     theme_void() +
#     facet_wrap(~year) 
#   
#   ggsave(base_full_map, 
#          filename = paste0("output/full_regic_comp", i, ".png"))
# }
# 
# for(i in min(base_full_comp$year):max(base_full_comp$year)) {
#   base_full_map <- ggplot(data = base_full_comp[base_full_comp$year == i,]) +
#     geom_sf(aes(fill = prior_full.diff), lwd = .05) +
#     scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
#                          mid = "white", high = "#c51b7d", midpoint = 0) +
#     expand_limits(fill = c(-4, 4)) +
#     theme_void() +
#     facet_wrap(~year) 
#   
#   ggsave(base_full_map, 
#          filename = paste0("output/full_prior_comp", i, ".png"))
# }

## Plot median over time period 
full_cov_comp_av <- group_by(df_model, municip_code_ibge) %>%
  summarise(clim_med = median(clim_full.diff),
            urb_med = median(urb_full.diff),
            regic_med = median(regic_full.diff),
            prior_med = median(prior_full.diff)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


# Climate (Figure 9a)
clim_diff_med_map <- ggplot(data = full_cov_comp_av) +
  geom_sf(aes(fill = clim_med), lwd = .05) +
  scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void() 

ggsave(clim_diff_med_map, filename = "output/clim_diff_av.png")


# Prior outbreak (Figure 9b)
prior_diff_med_map <- ggplot(data = full_cov_comp_av) +
  geom_sf(aes(fill = prior_med), lwd = .05) +
  scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void() 

ggsave(prior_diff_med_map, filename = "output/prior_diff_av.png")


# Urbanisation  (Figure  9c)
urb_diff_med_map <- ggplot(data = full_cov_comp_av) +
  geom_sf(aes(fill = urb_med), lwd = .05) +
  scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void() 

ggsave(urb_diff_med_map, filename = "output/urb_diff_av.png")


# REGIC (Figure 9d)
regic_diff_med_map <- ggplot(data = full_cov_comp_av) +
  geom_sf(aes(fill = regic_med), lwd = .05) +
  scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void() 

ggsave(regic_diff_med_map, filename = "output/regic_diff_av.png")



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
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


## Plot predicted probability per year (Figure S10)
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


## Plot mean probability per decade, Brazil (Figure 10a)
predictions_decade_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = mean_prob), col = NA) +
  facet_wrap(~ decade) +
  scale_fill_viridis_c(name = "Probability of\noutbreak", direction = 1) +
  coord_sf(datum = NA) +
  expand_limits(fill = c(0, 0.8)) +
  theme_void()

ggsave(predictions_decade_map, filename = "output/predictions_decade_map.png",
       height = 10, width = 15)

## Plot mean probability per decade, South Brazil (Figure 10b)
predictions_decade_sth <- ggplot(data = predictions_decade[substr(predictions_decade$municip_code_ibge, 1, 1) == "4",]) +
  geom_sf(aes(fill = mean_prob), col = NA) +
  facet_wrap(~ decade) +
  scale_fill_viridis_c(name = "Probability of\noutbreak", direction = 1) +
  coord_sf(datum = NA) +
  expand_limits(fill = c(0, 0.8)) +
  theme_void()

ggsave(predictions_decade_sth, filename = "output/predictions_decade_sth.png",
       height = 10, width = 15)

## Plot mean probability per decade, Acre + Amazonas (Figure 10c)
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
## Plot showing mean(P(outbreak)) = 0% (Figure S11a)
no_risk_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = no_risk), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "Status", labels = c("Protected", "Not protected"), 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()

ggsave(no_risk_map, filename = "output/risk0_map.png")


## Plot showing mean(P(outbreak)) = 5% (Figure S11b)
risk_5_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = risk_5), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "P(outbreak)", 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()

ggsave(risk_5_map, filename = "output/risk5_map.png")


## Plot showing mean(P(outbreak)) = 10% (Figure 11/Figure S11c)
risk_10_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = risk_10), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "P(outbreak)", 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()

ggsave(risk_10_map, filename = "output/risk_10_map.png",
       height = 10, width = 10)

## Plot showing mean(P(outbreak)) = 15% (Figure S11d)
risk_15_map <- ggplot(data = predictions_decade) +
  geom_sf(aes(fill = risk_15), lwd = .05) +
  facet_wrap(~ decade) +
  scale_fill_manual(name = "P(outbreak)", 
                    values = c("#68AA9C", "White")) +
  coord_sf(datum = NA) +
  theme_void()

ggsave(risk_15_map, filename = "output/risk15_map.png")






