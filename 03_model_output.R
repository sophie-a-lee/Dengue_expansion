#################################################
####                                         ####
####     Script to explore model results     ####
####                                         ####
#################################################


## Load packages and data 
source("00_load_packages_data.R")


#### Load model objects ####
## NOTE: This code will only work if 02_model_fitting.R has been run
model_base <- readRDS("output/model_base.rds")
model_full <- readRDS("output/model_full.rds")


## Sensitivity analysis - different variable options
model_full100 <- readRDS("output/model_full100.rds")
model_full_perc75 <- readRDS("output/model_full_perc75.rds")
model_aeg <- readRDS("output/model_aeg.rds")
model_wet <- readRDS("output/model_wet.rds")


## Sensitivity analysis - remove each variable in turn
model_temp <- readRDS("output/model_temp.rds")
model_prior <- readRDS("output/model_prior.rds")
model_urb <- readRDS("output/model_urb.rds")
model_regic <- readRDS("output/model_regic.rds")



#### Bayesian estimation of the beta parameters ####
# Function to obtain random numbers from a multivariate normal distribution
# n = number of number of random  draws
# mu = matrix of means
# sig = matrix of standard deviations
rmvn <- function(n, mu, sig) { 
  L <- mroot(sig); m <- ncol(L);
  t(mu + L %*% matrix(rnorm(m*n), m, n)) 
}


# Number of simulations
n.sims <- 1000


# Take 1000 random samples from the joint posterior of model coefficients (betas)
betas_full <- rmvn(n.sims, coef(model_full), model_full$Vp)
betas_full100 <- rmvn(n.sims, coef(model_full100), model_full100$Vp)
betas_perc75 <- rmvn(n.sims, coef(model_full_perc75), model_full_perc75$Vp)
betas_aeg <- rmvn(n.sims, coef(model_aeg), model_aeg$Vp)
betas_wet <- rmvn(n.sims, coef(model_wet), model_wet$Vp)



## Estimate beta mean and 95% credible interval (Table 1 & Table S2)
beta_ci <- function(model, betas) {
  
  # Extract number of fixed covariates
  ncov <- min(which(substr(names(coef(model)), 1, 2) == "s(")) - 1
  
  beta_ci <- as_tibble(betas) 
  # rename columns after coefficient (remove punctuation and spaces)
  names(beta_ci) <- gsub("[[:punct:][:blank:]]+","", names(model$coefficients))
  
  # Return the mean, 2.5th percentile and 97.5th percentile of the random draws from the beta posterior
  beta_ci <- beta_ci %>%
    summarise(across(2:ncov, list(mean = mean, lower = ~quantile(.x, 0.025), 
                               upper = ~quantile(.x, 0.975)),
                     names = "{.col}.fn{.fn}")) %>%
    # Transform the estimates onto the adjusted odds ratio scale
    mutate(across(everything(), ~exp(.x)))
  
  # Manipulate data into format that can be plotted
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

# Plot aOR and 95% CI
beta_plot_full <- ggplot(data = beta_ci_full, aes(Covariate)) +
  # Add reference line: aOR = 1 = no difference
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
  # Add aOR estimates
  geom_point(aes(x = Covariate, y = Mean, colour = Covariate), shape = 16, 
             size = 2, position = position_dodge(width = 0.5)) +
  # Add aOR 95% CIs
  geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Covariate), lwd = 1,
                 position = position_dodge(width = 0.5)) +
  # Add coefficient labels
  scale_x_discrete(name = "Model covariate",
                   labels = c("Months suitable", "Prior outbreak: Yes", 
                              "Metropolis", "Regional capital", 
                              "Sub-regional centre",  "Zone centre", 
                              "Urbanisation")) +
  # Add coour scheme (keep REGIC categorical variables same colour to aid interpretation)
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


# Medium risk model (outbreak > 100)
beta_ci_full100 <- beta_ci(model_full100, betas_full100)

# 75th percentile model (outbreak > 75th percentile)
beta_ci_perc75<- beta_ci(model_full_perc75, betas_perc75)


## Plot betas from alternative thresholds (Figure S8)
# Create df with aOR and 95% CI for each model, add Threshold variable 
beta_thresholds <- rbind(mutate(beta_ci_full, Covariate = paste0(Covariate, "300"),
                                Threshold = "DIR > 300"),
                         mutate(beta_ci_full100, Covariate = paste0(Covariate, "100"),
                                Threshold = "DIR > 100"),
                         mutate(beta_ci_perc75, Covariate = paste0(Covariate, "75"),
                                Threshold = "75th percentile"))


beta_comp_plot <- ggplot(data = beta_thresholds, aes(Covariate)) +
  # Add reference line: aOR = 1 = no difference
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5) +
  # Add aOR estimate, set colour different for each threshold
  geom_point(aes(x = Covariate, y = Mean, colour = Threshold), shape = 16, 
             size = 2, position = position_dodge(width = 0.5)) +
  # Add aOR 95% CI, set colour different for each threshold
  geom_linerange(aes(ymin = LCI, ymax = UCI, colour = Threshold), lwd = 1,
                 position = position_dodge(width = 0.5)) +
  # Add covariate labels (only one per covariate to avoid overcrowding on axis)
  scale_x_discrete(name = "Model covariate",
                   labels = c("", "Months suitable", "", 
                              "", "Prior outbreak: Yes", "",
                              "", "Metropolis", "",
                              "", "Regional capital", "",
                              "", "Sub-regional centre","",
                              "", "Zone centre", "", 
                              "", "Urbanisation", "")) +
  # Set colour per threshold model
  scale_colour_manual(values = c("#390099", "#9e0059", "#197278")) +
  coord_flip() + 
  scale_y_continuous(name = "Coefficient estimate (aOR)",
                     position = "right",
                     expand = c(0,0), limits = c(.8, 4)) +
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))  

# Extract legend
threshold_leg <- get_legend(beta_comp_plot)

# Add legend to the pllot  area (to avoid overlapping on the original plot)
beta_comp_plot <- plot_grid(beta_comp_plot + theme(legend.position = "none"),
                            threshold_leg,
                            rel_widths = c(3, .4))

ggsave(beta_comp_plot, filename = "output/beta_comp_plot.png",
       width = 10, height = 5)


#### Plot smooth functions ####
## Temporal smooth using mgcViz (Figure 7a)
model_full_viz <- getViz(model_full)

temp_plot.full <- plot(sm(model_full_viz, 2), trans = exp) +
  labs(x = "Year", y = "s(Year)") +
  l_fitLine( ) +
  l_ciLine() +
  # Add reference lline for aOR = 1 = no difference
  geom_hline(yintercept = 1, linetype = "dashed") 


## Plot the spatial smooth function
# Return a matrix of linear predictors (including for smooth functions)
model_mat.full <- predict(object = model_full, newdata = df_model, 
                          type = "lpmatrix")

# Select linear predictors for only spatial parameters (named s(lon, lat))
spat_mat.full <- model_mat.full[ , substr(names(coefficients(model_full)), 1, 5) == "s(lon"]
# Select coefficient estimates for only spatial parameters 
betas_spat.full<- betas_full[ , substr(names(coefficients(model_full)), 1, 5) == "s(lon"]


# Estimate spatial random effects (multiply coefficient values by linear predictors) + transform from log scale
spat_est.full <- exp(spat_mat.full %*% t(betas_spat.full))
# As all years are equal, just keep 2001
spat_est01.full <- spat_est.full[df_model$year == 2001, ] 

# Take average of simulations
mean_spat.full <- apply(spat_est01.full, 1, mean)

# Create 95% credible interval by taking 2.5th and 97.5th percentiles
lci_spat.full <- apply(spat_est01.full, 1, quantile, probs=0.025)
uci_spat.full <- apply(spat_est01.full, 1, quantile, probs=0.975)


# Join random effects to model df + plot on a map
spat_smooth.full <- data.table(municip_code_ibge = df_model[df_model$year == 2001,]$municip_code_ibge,
                               spat_smooth = mean_spat.full,
                               spat_re_uci = uci_spat.full,
                               spat_re_lci =  lci_spat.full) %>%
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
  
  # Obtain simulations from beta posterior distributions 
  betas <- rmvn(1000, coef(model), model$Vp) 
  # Only keep estimates from smooth functions (beginning 's(' or 'ti(')
  betas_smooth <- betas[ , substr(names(coefficients(model)), 1, 2) %in%
                           c("s(", "ti")]
  
  # Obtain linear predictor matrix for estimates
  model_mat <- predict(object = model, newdata = data, type = "lpmatrix")
  # Only keep smooth functions
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

# Find the absolute difference between smooth functions (shrinkage towards zero)
base_full_diff <- abs(smooth_final) - abs(smooth_base)

# Keep average 
base_full_diff.med <- apply(base_full_diff, 1, median)



## Compare baseline and final models, plot absolute differences (Figure 8)
# (if final model is better, would expect reduction in absolute values/shrinkage to 0)
base_full_comp <- df_model %>% 
  mutate(smooth_absdiff = base_full_diff.med) %>% 
  left_join(df_model, shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


# df with average absolute difference per municipality
beta_full_comp_av <- group_by(df_model, municip_code_ibge) %>%
  summarise(diff_med = median(smooth_absdiff)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


# Plot average difference, green = final model performs better
base_full_med <- ggplot(data = beta_full_comp_av) +
  geom_sf(aes(fill = diff_med), lwd = .05) +
  scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void() 

ggsave(base_full_med, file = "output/base_full_med.png")


#### Compare smooth functions after adding covariates in turn ####
# Urbanisation
smooth_urb <- smooth_estimates(model_urb, df_model)

# Estimate absolute difference 
urb_full_diff <- abs(smooth_urb) - abs(smooth_base)

# Take average of simulations and add to df
urb_full_diff.med <- apply(urb_full_diff, 1, median)
df_model$urb_full.diff <- urb_full_diff.med


# Prior outbreak
smooth_prior <- smooth_estimates(model_prior, df_model)

# Estimate absolute difference 
prior_full_diff <- abs(smooth_prior) - abs(smooth_base)

# Take average of simulations and add to df
prior_full_diff.med <- apply(prior_full_diff, 1, median)
df_model$prior_full.diff <- prior_full_diff.med


# Connectivity
smooth_regic <- smooth_estimates(model_regic, df_model)

# Estimate absolute difference 
regic_full_diff <- abs(smooth_regic) - abs(smooth_base)

# Take average of simulations and add to df
regic_full_diff.med <- apply(regic_full_diff, 1, median)
df_model$regic_full.diff <- regic_full_diff.med


# Number of months suitable
smooth_temp <- smooth_estimates(model_temp, df_model)

# Estimate absolute difference 
temp_full_diff <- abs(smooth_temp) - abs(smooth_base)

# Take average of simulations and add to df
temp_full_diff.med <- apply(temp_full_diff, 1, median)
df_model$temp_full.diff <- temp_full_diff.med


## Find average difference for each municipality for each covariate
full_cov_comp_av <- group_by(df_model, municip_code_ibge) %>%
  summarise(temp_med = median(temp_full.diff),
            urb_med = median(urb_full.diff),
            regic_med = median(regic_full.diff),
            prior_med = median(prior_full.diff)) %>%
  ungroup() %>%
  left_join(., shp_parent, by = "municip_code_ibge") %>%
  st_as_sf()


## Plot average difference in smooth functions for each covariate
# Temperature (Figure 9a)
temp_diff_med_map <- ggplot(data = full_cov_comp_av) +
  geom_sf(aes(fill = temp_med), lwd = .05) +
  scale_fill_gradient2(name = "Difference in\nsmooth terms", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void() 

ggsave(temp_diff_med_map, filename = "output/temp_diff_av.png")


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


## Calculate % of municipalities with shrinkage to 0(improvement compared to base model)
# Full vs base model
sum(beta_full_comp_av$diff_med <0)/
  length(unique(beta_full_comp_av$municip_code_ibge)) * 100

# Temperature vs base
sum(full_cov_comp_av$temp_med < 0)/
  length(unique(full_cov_comp_av$municip_code_ibge)) * 100 # 91.16%

#  Prior vs base
sum(full_cov_comp_av$prior_med < 0)/
  length(unique(full_cov_comp_av$municip_code_ibge)) * 100 # 94.28

# Urbanisation vs base
sum(full_cov_comp_av$urb_med < 0)/
  length(unique(full_cov_comp_av$municip_code_ibge)) * 100 # 57.50%

# REGIC vs base
sum(full_cov_comp_av$regic_med < 0)/
  length(unique(full_cov_comp_av$municip_code_ibge)) * 100 # 45.08%


#### Compare model fit between threshold models (Figure S7) ####
## Function to obtain predictions of the probability of an outbreak 
model_predict <- function(model, betas) {
  
  # Extract matrix of linear predictors from model
  model_mat <- predict(object = model, newdata = df_model, type = "lpmatrix")
  
  ## Use linear prediction matrix and parameter simulations to estimate posterior mean
  # Convert using the probit function to estimate the probability of an outbreak
  mean_prob <- exp(model_mat %*% t(betas))/(1 + exp(model_mat %*% t(betas)))
  
  ## Simulate from the posterior distribution (using mean probability, draw from a bernoulli distribution)
  preds <- apply(mean_prob, 1, function(x){rbinom(length(x), size = 1, prob = x)})
  
  ## Obtain mean probability (number of outbreaks/total number of draws)
  mean_pred <- apply(preds, 2, mean)
  
  return(mean_pred)
}


# Predict the probability of an outbreak for each outbreak threshold model
model_full_pred <- model_predict(model_full, betas_full)
model_full100_pred <- model_predict(model_full100, betas_full100)
model_perc75_pred <- model_predict(model_full_perc75, betas_perc75)

  
## Convert into a dataset to plot  
predictions <- data.table(municip_code_ibge = df_model$municip_code_ibge,
                          year = df_model$year,
                          outbreak_obs300 = df_model$outbreak_fix,
                          outbreak_obs100 = df_model$outbreak_fix100,
                          outbreak_obs_perc75 =  df_model$outbreak_perc75,
                          prob_pred300 = model_full_pred,
                          prob_pred100 = model_full100_pred,
                          prob_pred_perc75 = model_perc75_pred)  %>%
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


#### Model checking using the ROC curve ####
# Create ROC object based on observed outbreak status and predicted predictions
# This also returns area under the ROC curve with 95% CI
roc_obj300 <- roc(predictions$outbreak_obs300, predictions$prob_pred300,
               auc = T, ci = T, plot = T)
roc_obj100 <- roc(predictions$outbreak_obs100, predictions$prob_pred100,
               auc = T, ci = T, plot = T)
roc_obj75 <- roc(predictions$outbreak_obs_perc75, predictions$prob_pred_perc75,
               auc = T, ci = T, plot = T)

# Combine ROC data into one df (for plotting)
roc_obj <- data.table(roc300_sens = roc_obj300$sensitivities,
                      roc300_spec = 1 - roc_obj300$specificities,
                      roc100_sens = roc_obj100$sensitivities,
                      roc100_spec = 1 - roc_obj100$specificities,
                      roc75_sens = roc_obj75$sensitivities,
                      roc75_spec = 1 - roc_obj75$specificities)


## Plot ROC curve (Figure S7) ##
roc_curve <- ggplot(data = roc_obj) +
  # Add a line per threshold
  geom_line(aes(x = roc300_spec,  y = roc300_sens)) +
  geom_line(aes(x = roc100_spec,  y = roc100_sens), linetype = "dashed",
            col = "red") +
  geom_line(aes(x = roc75_spec,  y = roc75_sens), linetype = "dotdash",
            col = "blue") +
  # Add reference line (y = x line represents chance)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "True negative rate", y = "True positive rate") +
  theme_light()

ggsave(roc_curve, filename = "output/ROC_curve.png")


## Calculate Brier scores (Table S3) 
mean((predictions$prob_pred300 - predictions$outbreak_obs300)^2)
mean((predictions$prob_pred100 - predictions$outbreak_obs100)^2)
mean((predictions$prob_pred_perc75 - predictions$outbreak_obs_perc75)^2)

## Calculate model diagnostics for  extremely wet model
# Return predicted probabilities
prob_pred_wet <- model_predict(model_wet, betas_wet)

# Return area under the ROC curve (with 95% CI)
roc(predictions$outbreak_obs300, prob_pred_wet,
    auc = T, ci = T, plot = T)
# Calculllate Brier score
mean((prob_pred_wet - predictions$outbreak_obs300)^2)


## Calculate average probability of an outbreak 2001 - 2010 vs 2011 - 2020 
predictions_decade <- st_drop_geometry(predictions) %>%
  # Create 'decade' factor
  mutate(decade = cut(year, 
                      breaks = c(-Inf, 2010, Inf), 
                      labels = c("2001 - 2010", "2011 - 2020"))) %>%
  group_by(municip_code_ibge, decade) %>%
  # Estimate average probability per municipality per deacde
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






