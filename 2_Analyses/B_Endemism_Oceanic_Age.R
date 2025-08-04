# This script:
# 1. models and plots results for drivers of endemic richness
# 2. models and plots results for drivers of non-endemic occurrence and diversification (where present) associated with speciation and coexistence (without mutualism data)
# 3. models and plots results for drivers of non-endemic occurrence and diversification (where present) associated with speciation and coexistence for each mutualism type (4)
# Outputs for all models: model averaging and best model

############################
###### LOAD PACKAGES #######
############################

library(mgcv); library(gridExtra); library(betareg); library(MASS); library(lme4); library(lmerTest); library(lsmeans); library(ggeffects); library(spdep); library(ggplot2); library(ncf); library(ape); library(sjPlot); library(gridExtra); library(MuMIn);library(maps); library(sf); library(car);library(viridis);library(tidyverse);library(glmmTMB);library(DHARMa);library(purrr)
options(na.action = "na.fail")

############################
#### PLOTTING FUNCTION #####
############################

# set function for plotting single vars
predict_over_variable <- function(var, model, data, n = 1000) {
  vars <- all.vars(formula(model))[-1]  # get predictor names
  base_means <- summarise(na.omit(data[, vars]), across(everything(), ~mean(.x)))
  
  # Gendnonrate sequence for the target variable
  seq_var <- seq(
    max(data[[var]], na.rm = TRUE),
    min(data[[var]], na.rm = TRUE),
    length.out = n
  )
  
  # Repeat the mean row and replace the target column with the sequence
  new_data <- base_means[rep(1, n), ]
  new_data[[var]] <- seq_var
  
  # Predict
  pred <- predict(model, newdata = new_data, se.fit = TRUE, type = "response")
  
  return(list(pred = pred, new_data = new_data, variable = var, sequence = seq_var))
}

####################################
####################################
####### NON-MUTUALIST MODELS #######
####################################
####################################

# set manual color
colScale <- scale_colour_manual(values =c ("rosybrown2","lightcoral", "darkgrey"))
fillScale <- scale_fill_manual(values =c ("rosybrown2","lightcoral",  "darkgrey"))

# end/non dat
dat.endnon <- readRDS("data/endnonend_data_2025.RDS") %>%
  select(c("entity_ID","age_Ma"))

# endemisp data quality
dat.check <- readRDS("data/endnonend_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# clean data
debt <- readRDS("data/deficit_latitude_data_2025.RDS") %>%
  filter(entity_class2 == "Oceanic") %>%
  left_join(dat.check, by = "entity_ID") %>%
  left_join(dat.endnon, by = "entity_ID") %>%
  mutate(pend = endemic/(endemic + nonendemic)) %>%
  mutate(logitend = log(endemic + 0.1/(endemic+ 0.1 + nonendemic+ 0.1)/(nonendemic+ 0.1/(endemic+ 0.1 + nonendemic+ 0.1)))) %>%
  mutate(log_endemic = log10(endemic)) %>%
  mutate(end_exists = as.integer(ifelse(endemic>0, 1,0))) %>%
  mutate(
    us_age_Ma = age_Ma,
    us_area = area, 
    us_dist = dist,
    us_temp = temp,
    us_prec = prec,
    us_elev_range = elev_range) %>%
  mutate(
    age_Ma = as.vector(scale(log10((age_Ma)+.01))),
    Tdebt = as.vector(scale(log10((Tdebt + 1)))), 
    sc_nonendemic = as.vector(scale(log10((nonendemic)+.01))),
    area = as.vector(scale(log10((area)+.01))), 
    dist = as.vector(scale(dist)),
    temp = as.vector(scale(log10((temp)+.01))),
    prec = as.vector(scale(log10((prec)+.01))),
    elev_range = as.vector(scale(log10((elev_range)+.01))),
    abslatitude = as.vector(scale(log10((abslatitude)+.01))))  %>%
    drop_na() %>%
    filter(!entity_ID == "984")

# Colinearity check
debt.co <- debt %>%                                                                                       
  dplyr::select(c("Tdebt","sc_nonendemic","age_Ma","area","dist","elev_range","abslatitude", "temp", "prec"))                               
debt.cor <- as.matrix(cor(debt.co))  

#########################
## NON-ENDEMIC COUNTS ###
#########################

dat <- debt
theta_start <- theta.ml(dat$endemic, mean(dat$endemic))
mod <- glm.nb(nonendemic ~ age_Ma + dist + area + elev_range + temp + prec,
                weights = pnata, data = dat, init.theta = theta_start) 
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Broad_NE_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Broad_NE_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model)

#########################
##### PRED END OCC ######
#### ONLY NONENDEMIC ####
#########################

dat <- debt
mod <- glm(end_exists ~ sc_nonendemic + age_Ma + dist + area + elev_range +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range,
             family = binomial, weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Broad_YN_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Broad_YN_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model) 

#########################
##### PRED END OCC ######
####### ONLY TDEBT ######
#########################

#dat <- debt 
#mod <- glm(end_exists  ~ Tdebt + age_Ma + dist + area + elev_range +
#             Tdebt:age_Ma + Tdebt:dist + Tdebt:area + Tdebt:elev_range,
#             family = binomial, weights = pnata, data = dat) 
#summary(mod)

#dredged_models <- dredge(mod, rank = AICc)
#averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
#saveRDS(averaged_model_results, "models/Broad_YN_ENDTdebt_averaged_model_results.RDS")
#averaged_model_results <- readRDS("models/Broad_YN_ENDTdebt_averaged_model_results.RDS")
#summary(averaged_model_results) 

#best_model <- get.models(averaged_model_results, subset = 1)[[1]]
#summary(best_model) 

#########################
##### PRED END OCC ######
### NONENDEMIC & TDEBT ##
#########################

#dat <- debt 
#mod <- glm(end_exists  ~ sc_nonendemic + age_Ma + dist + area + elev_range +
#             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range +
#             Tdebt:age_Ma + Tdebt:dist + Tdebt:area + Tdebt:elev_range,
#             family = binomial, weights = pnata, data = dat) 
#summary(mod)

#dredged_models <- dredge(mod, rank = AICc)
#averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
#saveRDS(averaged_model_results, "models/Broad_YN_ENDNonEndTdebt_averaged_model_results.RDS")
#averaged_model_results <- readRDS("models/Broad_YN_ENDNonEndTdebt_averaged_model_results.RDS")
#summary(averaged_model_results) 

#best_model <- get.models(averaged_model_results, subset = 1)[[1]]
#summary(best_model) 

#########################
#### PRED END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

dat <- debt %>% filter(end_exists == 1) 
mod <- lm(logitend ~ sc_nonendemic + age_Ma + dist + area + elev_range +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range,
             weights = pnata, data = dat) 
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Broad_DIVERS_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Broad_DIVERS_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model) 

#########################
#### PRED END DIVERS ####
####### ONLY TDEBT ######
#########################

#dat <- debt %>% filter(end_exists == 1) 
#mod <- lm(logitend ~ Tdebt + age_Ma + dist + area + elev_range +
#           Tdebt:age_Ma + Tdebt:dist + Tdebt:area + Tdebt:elev_range,
#           weights = pnata, data = dat) 
#summary(mod)

#dredged_models <- dredge(mod, rank = AICc)
#averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
#saveRDS(averaged_model_results, "models/Broad_DIVERS_ENDTdebt_averaged_model_results.RDS")
#averaged_model_results <- readRDS("models/Broad_DIVERS_ENDTdebt_averaged_model_results.RDS")
#summary(averaged_model_results) 

#best_model <- get.models(averaged_model_results, subset = 1)[[1]]
#summary(best_model) 

#########################
#### PRED END DIVERS ####
### NONENDEMIC & TDEBT ##
#########################

#dat <- debt %>% filter(end_exists == 1) 
#mod <- lm(logitend ~ sc_nonendemic + age_Ma + dist + area + elev_range +
#             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range +
#             Tdebt:age_Ma + Tdebt:dist + Tdebt:area + Tdebt:elev_range,
#             weights = pnata, data = dat) 
#summary(mod)

#dredged_models <- dredge(mod, rank = AICc)
#averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
#saveRDS(averaged_model_results, "models/Broad_DIVERS_ENDNonEndTdebt_averaged_model_results.RDS")
#averaged_model_results <- readRDS("models/Broad_DIVERS_ENDNonEndTdebt_averaged_model_results.RDS")
#summary(averaged_model_results) 

#best_model <- get.models(averaged_model_results, subset = 1)[[1]]
#summary(best_model) 

#########################
######### PLOT ##########
#########################

#########################
######## PLOT NE ########
#########################

########## AGE ##########
result <- predict_over_variable("age_Ma", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$age_Ma), max(new_data$age_Ma), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_age_Ma + 0.01), na.rm = TRUE) + mean(log10(dat$us_age_Ma + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = age_Ma, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = age_Ma, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = age_Ma, y = nonendemic), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Colonist richness") +
  xlab("Age") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_NE_age_Ma.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## AREA ##########
result <- predict_over_variable("area", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$area), max(new_data$area), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_area + 0.01), na.rm = TRUE) + mean(log10(dat$us_area + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = area, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = area, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = area, y = nonendemic), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Colonist richness") +
  xlab("Area") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_NE_area.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## DIST ##########
result <- predict_over_variable("dist", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

min_scaled_for_positive <- -mean(dat$us_dist, na.rm = TRUE) / sd(dat$us_dist, na.rm = TRUE)
x_breaks <- seq(min_scaled_for_positive + 0.1, max(new_data$dist), length.out = 5)
x_labels <- x_breaks * sd(dat$us_dist, na.rm = TRUE) + mean(dat$us_dist, na.rm = TRUE)

plot <- ggplot() +
  geom_line(data = new_data, aes(x = dist, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = dist, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = dist, y = nonendemic), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Colonist richness") +
  xlab("Distance") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_NE_dist.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## ELEV RANGE ##########
result <- predict_over_variable("elev_range", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$elev_range), max(new_data$elev_range), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_elev_range + 0.01), na.rm = TRUE) + mean(log10(dat$us_elev_range + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = elev_range, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = elev_range, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = elev_range, y = nonendemic), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Colonist richness") +
  xlab("Elevation range") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_NE_elev_range.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## PREC ##########
result <- predict_over_variable("prec", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$prec), max(new_data$prec), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_prec + 0.01), na.rm = TRUE) + mean(log10(dat$us_prec + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = prec, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = prec, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = prec, y = nonendemic), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Colonist richness") +
  xlab("Precipitation") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_NE_prec.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## TEMP ##########
result <- predict_over_variable("temp", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$temp), max(new_data$temp), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_temp + 0.01), na.rm = TRUE) + mean(log10(dat$us_temp + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = temp, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = temp, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = temp, y = nonendemic), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Colonist richness") +
  xlab("Temperature") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_NE_temp.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

#########################
##### PLOT END OCC ######
#### ONLY NONENDEMIC ####
#########################

########## AREA #########
result <- predict_over_variable("area", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$area), max(new_data$area), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_area + 0.01), na.rm = TRUE) + mean(log10(dat$us_area + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = area, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = area, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = area, y = end_exists), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic occurence") +
  xlab("Area") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_YN_ENDNonEnd_area.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## NONENDEMIC ##########
result <- predict_over_variable("sc_nonendemic", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$sc_nonendemic), max(new_data$sc_nonendemic), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$nonendemic + 0.01), na.rm = TRUE) + mean(log10(dat$nonendemic + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = sc_nonendemic, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = sc_nonendemic, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = sc_nonendemic, y = end_exists), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic occurence") +
  xlab("Colonist richness") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_YN_ENDNonEnd_sc_nonendemic.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## AGE ######### 
###### NS AVG MOD ######
result <- predict_over_variable("age_Ma", averaged_model_results, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$age_Ma), max(new_data$age_Ma), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_age_Ma + 0.01), na.rm = TRUE) + mean(log10(dat$us_age_Ma + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = age_Ma, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = age_Ma, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.2) +  
  geom_point(data = dat, aes(x = age_Ma, y = end_exists), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic occurence") +
  xlab("Age") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_YN_ENDNonEnd_age_Ma.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## DIST ########
###### NS AVG MOD ######
result <- predict_over_variable("dist", averaged_model_results, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

min_scaled_for_positive <- -mean(dat$us_dist, na.rm = TRUE) / sd(dat$us_dist, na.rm = TRUE)
x_breaks <- seq(min_scaled_for_positive + 0.1, max(new_data$dist), length.out = 5)
x_labels <- x_breaks * sd(dat$us_dist, na.rm = TRUE) + mean(dat$us_dist, na.rm = TRUE)

plot <- ggplot() +
  geom_line(data = new_data, aes(x = dist, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = dist, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.2) +  
  geom_point(data = dat, aes(x = dist, y = end_exists), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic occurence") +
  xlab("Distance") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_YN_ENDNonEnd_dist.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

###### ELEV RANGE ######
###### NS AVG MOD ######
result <- predict_over_variable("elev_range", averaged_model_results, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$elev_range), max(new_data$elev_range), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_elev_range + 0.01), na.rm = TRUE) + mean(log10(dat$us_elev_range + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = elev_range, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = elev_range, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.2) +  
  geom_point(data = dat, aes(x = elev_range, y = end_exists), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic occurence") +
  xlab("Elevation range") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_YN_ENDNonEnd_elev_range.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

#########################
#### PLOT END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

########## AGE ##########
result <- predict_over_variable("age_Ma", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$age_Ma), max(new_data$age_Ma), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_age_Ma + 0.01), na.rm = TRUE) + mean(log10(dat$us_age_Ma + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = age_Ma, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = age_Ma, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = age_Ma, y = logitend), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") +
  xlab("Age") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_DIVERS_ENDNonEnd_age_Ma.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## NONENDEMIC ##########
result <- predict_over_variable("sc_nonendemic", best_model, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$sc_nonendemic), max(new_data$sc_nonendemic), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$nonendemic + 0.01), na.rm = TRUE) + mean(log10(dat$nonendemic + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = sc_nonendemic, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = sc_nonendemic, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.5) +  
  geom_point(data = dat, aes(x = sc_nonendemic, y = logitend), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") +
  xlab("Colonist richness") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_DIVERS_ENDNonEnd_sc_nonendemic.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## AREA*NONENDEMIC ##########
sc_nonendemic.range <- seq(from = quantile(dat$sc_nonendemic, 0.001), to = quantile(dat$sc_nonendemic, 0.999), length.out = 100)
area.minmax <- data.frame(area = quantile(dat$area, probs = c(0.25, 0.75)), level = c("low", "high"))
ex.grid <- expand.grid(sc_nonendemic = sc_nonendemic.range, area = area.minmax$area)

predictor_vars <- c("temp", "prec", "dist", "pnata", "age_Ma", "elev_range")
pred.dat <- cbind(ex.grid, t(colMeans(na.omit(dat[, predictor_vars]), na.rm = TRUE)))
pred <- predict(best_model, newdata = pred.dat, se.fit = TRUE, type = "response")

pdat <- cbind(pred.dat, fit = pred$fit, se.fit = pred$se.fit) %>%
  mutate(conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit) %>%
  left_join(area.minmax, by = "area")
pdat$level <- factor(pdat$level, levels = c("low", "high"))

x_breaks <- seq(min(pdat$sc_nonendemic), max(pdat$sc_nonendemic), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$nonendemic + 0.01), na.rm = TRUE) + mean(log10(dat$nonendemic + 0.01), na.rm = TRUE)) - 0.01

dat <- dat %>%
  mutate(level = case_when(
    area <= area.minmax[1,1] ~ "low",
    area >= area.minmax[2,1] ~ "high",
    TRUE ~ "other"))

plot <- 
  ggplot(pdat) + 
  geom_line(aes(x = sc_nonendemic, y = fit, color = level)) + 
  geom_ribbon(aes(x = sc_nonendemic, ymin = fit - se.fit, ymax = fit + se.fit, fill = level), alpha = 0.8) +  
  geom_point(data = dat, aes(x = sc_nonendemic, y = logitend, color = level)) +  
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") +
  xlab("Colonist richness") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") +
  colScale + fillScale

png("figures/Broad_DIVERS_ENDNonEnd_areabysc_nonendemic.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## DIST*NONENDEMIC ##########
sc_nonendemic.range <- seq(from = quantile(dat$sc_nonendemic, 0.001), to = quantile(dat$sc_nonendemic, 0.999), length.out = 100)
dist.minmax <- data.frame(dist = quantile(dat$dist, probs = c(0.25, 0.75)), level = c("low", "high"))
ex.grid <- expand.grid(sc_nonendemic = sc_nonendemic.range, dist = dist.minmax$dist)

predictor_vars <- c("temp", "prec", "area", "pnata", "age_Ma", "elev_range")
pred.dat <- cbind(ex.grid, t(colMeans(na.omit(dat[, predictor_vars]), na.rm = TRUE)))
pred <- predict(best_model, newdata = pred.dat, se.fit = TRUE, type = "response")

pdat <- cbind(pred.dat, fit = pred$fit, se.fit = pred$se.fit) %>%
  mutate(conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit) %>%
  left_join(dist.minmax, by = "dist")
pdat$level <- factor(pdat$level, levels = c("low", "high"))

x_breaks <- seq(min(pdat$sc_nonendemic), max(pdat$sc_nonendemic), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$nonendemic + 0.01), na.rm = TRUE) + mean(log10(dat$nonendemic + 0.01), na.rm = TRUE)) - 0.01

dat <- dat %>%
  mutate(level = case_when(
    dist <= dist.minmax[1,1] ~ "low",
    dist >= dist.minmax[2,1] ~ "high",
    TRUE ~ "other"))

plot <- 
  ggplot(pdat) + 
  geom_line(aes(x = sc_nonendemic, y = fit, color = level)) + 
  geom_ribbon(aes(x = sc_nonendemic, ymin = fit - se.fit, ymax = fit + se.fit, fill = level), alpha = 0.8) +  
  geom_point(data = dat, aes(x = sc_nonendemic, y = logitend, color = level)) +  
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") +
  xlab("Colonist richness") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") +
  colScale + fillScale

png("figures/Broad_DIVERS_ENDNonEnd_distbysc_nonendemic.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

###### ELEV RANGE ######
###### NS AVG MOD ######

result <- predict_over_variable("elev_range", averaged_model_results, dat)
new_data <- result$new_data %>%
  mutate(fit = result$pred$fit, se.fit = result$pred$se.fit, conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(new_data$elev_range), max(new_data$elev_range), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$us_elev_range + 0.01), na.rm = TRUE) + mean(log10(dat$us_elev_range + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = new_data, aes(x = elev_range, y = fit), color = "darkgrey", size = 1) +  
  geom_ribbon(data = new_data, aes(x = elev_range, ymin = conf.low, ymax = conf.high), 
              fill = "darkgrey", alpha = 0.2) +  
  geom_point(data = dat, aes(x = elev_range, y = logitend), color = "darkgrey", size = 3, shape = 16, alpha = 0.5) +  
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") +
  xlab("Elevation range") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position = "none") 

png("figures/Broad_DIVERS_ENDNonEnd_elev_range.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

####################################
####################################
############### POLL ###############
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "polltype", values = c("deeppink3","darkgrey"))
colFill <- scale_fill_manual(name = "polltype", values = c("deeppink3","darkgrey"))

colScale.pt <- scale_colour_manual(name = "polltype", values = c("darkgrey","deeppink3"))
colFill.pt <- scale_fill_manual(name = "polltype", values = c("darkgrey", "deeppink3"))

# end/non dat
p.endnon <- readRDS("data/merge_poll_data_2025.RDS") %>%
  select(c("entity_ID","flora","poll.b","poll.ab","age_Ma")) %>%
  filter(!flora == "native") %>%
  pivot_longer(cols = c(poll.b, poll.ab), names_to = "polltype", values_to = "spcount") %>%
  pivot_wider(names_from = flora, values_from = spcount) 

# exp rich dat
p.expsprich <- readRDS("data/deficit_poll_latitude_data_2025.RDS") %>%
  select(c("entity_ID","poll.b_exp","poll.ab_exp")) %>%
  pivot_longer(cols = c(poll.b_exp, poll.ab_exp), names_to = "polltype", values_to = "exp_sp") %>%
  mutate(polltype = ifelse(polltype == "poll.b_exp", "poll.b", "poll.ab"))

# endemism data quality
p.check <- readRDS("data/poll_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# clean data
p.debt <- readRDS("data/deficit_poll_latitude_data_2025.RDS") %>%
  select(-c("poll.b_exp", "poll.ab_exp", "propb_exp", "poll.b.diff", "poll.ab.diff", "sprichdiff", "Tdebt","C_poll.b.debt", "C_poll.ab.debt", "poll.b", "poll.ab")) %>%
  filter(entity_class2 == "Oceanic") %>%
  left_join(p.endnon, by ="entity_ID") %>%
  left_join(p.check, by = "entity_ID") %>%
  left_join(p.expsprich, by = c("entity_ID", "polltype")) %>%
  mutate(pend = endemic/(endemic + nonendemic)) %>%
  mutate(logitend = log(endemic + 0.1/(endemic+ 0.1 + nonendemic+ 0.1)/(nonendemic+ 0.1/(endemic+ 0.1 + nonendemic+ 0.1)))) %>%
  mutate(end_exists = as.integer(ifelse(endemic>0, 1,0))) %>%
  mutate(
    age_Ma = as.vector(scale(log10((age_Ma)+.01))),
    sc_nonendemic = as.vector(scale(log10((nonendemic)+.01))),
    area = as.vector(scale(log10((area)+.01))), 
    dist = as.vector(scale(dist)),
    temp = as.vector(scale(log10((temp)+.01))),
    prec = as.vector(scale(log10((prec)+.01))),
    elev_range = as.vector(scale(log10((elev_range)+.01))),
    abslatitude = as.vector(scale(log10((abslatitude)+.01)))) %>%
  drop_na() 

#########################
##### PRED END OCC ######
#### ONLY NONENDEMIC ####
#########################

dat <- p.debt
mod <- glm(end_exists ~ polltype + sc_nonendemic + age_Ma + dist + area + elev_range + #(1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:polltype + 
             polltype:age_Ma + polltype:dist + polltype:area + polltype:elev_range,
             family = binomial, weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Poll_YN_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Poll_YN_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model)

#########################
#### PRED END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

dat <- p.debt %>% filter(end_exists == 1) 
mod <- lmer(logitend ~ polltype + sc_nonendemic + age_Ma + dist + area + elev_range + (1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:polltype + 
             polltype:age_Ma + polltype:dist + polltype:area + polltype:elev_range,
             weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Poll_DIVERS_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Poll_DIVERS_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model) 

# find simple slope by muttype
emtrends(best_model, ~ polltype, var="sc_nonendemic")

#########################
######### PLOT ##########
#########################

#########################
##### PLOT END OCC ######
#### ONLY NONENDEMIC ####
#########################

####### POLLTYPE ########
ref <- lsmeans(best_model, pairwise ~ polltype, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)

plot <- 
  ggplot(ref.table, aes(reorder(polltype, lsmean), lsmean, color = polltype, fill = polltype)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(polltype, end_exists), position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Endemic occurence") + 
  xlab("Mutualist type") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  scale_x_discrete(labels = c("Abiotic", "Biotic")) +
  colFill.pt + colScale.pt 

png("figures/Poll_YN_ENDNonEnd_polltype.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

#########################
#### PLOT END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

####### POLLTYPE ########
ref <- emmeans(best_model, pairwise ~ polltype, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$emmeans)

plot <- 
  ggplot(ref.table, aes(polltype, emmean, color = polltype, fill = polltype)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(polltype, logitend), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Endemic prevalence") + 
  xlab("Mutualist type") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  scale_x_discrete(labels = c("Abiotic", "Biotic")) +
  colFill.pt + colScale.pt +
  ylim(-0.5,3)

png("figures/Poll_DIVERS_ENDNonEnd_polltype.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

##### POLLTYPE*NONENDEMIC ######
sc_nonendemic.range <- seq(max(dat$sc_nonendemic),min(dat$sc_nonendemic),length.out=1000)
polltype.levels <- unique(dat$polltype)
ex.grid <- expand.grid(sc_nonendemic = sc_nonendemic.range, polltype = polltype.levels)

predictor_vars <- c("temp", "prec", "area", "pnata", "age_Ma", "elev_range", "dist")
pred.dat <- cbind(ex.grid, t(colMeans(na.omit(dat[, predictor_vars]), na.rm = TRUE)))
pred <- predict(best_model, newdata = pred.dat, se.fit = TRUE, type = "response", re.form = NA)

pred.nonend <- cbind(pred.dat, fit = pred$fit, se.fit = pred$se.fit) %>%
  mutate(conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit) 

x_breaks <- seq(min(pred.nonend$sc_nonendemic), max(pred.nonend$sc_nonendemic), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$nonendemic + 0.01), na.rm = TRUE) + mean(log10(dat$nonendemic + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = pred.nonend, mapping = aes(x = sc_nonendemic, y = fit, color = polltype, fill = polltype), size = 1) +
  geom_ribbon(data = pred.nonend, aes(x = sc_nonendemic, ymin = conf.low, ymax = conf.high, color = polltype, fill = polltype),  alpha = 0.5) +
  geom_point(data = dat,aes(x = sc_nonendemic, y = logitend, color = polltype, fill = polltype),size = 3,shape = 16, alpha = 0.5) +
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") + 
  xlab("Colonist richness") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position="none") +
  colFill + colScale +
  ylim(-2.5,6)

png("figures/Poll_DIVERS_ENDNonEnd_polltypebysc_nonendemic.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

####################################
####################################
############### DISP ###############
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "disptype" , values = c("brown4","darkgrey"))
colFill <- scale_fill_manual(name = "disptype", values = c("brown4","darkgrey"))

colScale.pt <- scale_colour_manual(name = "disptype" , values = c("darkgrey", "brown4"))
colFill.pt <- scale_fill_manual(name = "disptype", values = c("darkgrey", "brown4"))

# end/non dat
d.endnon <- readRDS("data/merge_disp_data_2025.RDS") %>%
  select(c("entity_ID","flora","disp.b","disp.ab","age_Ma")) %>%
  filter(!flora == "native") %>%
  pivot_longer(cols = c(disp.b, disp.ab), names_to = "disptype", values_to = "spcount") %>%
  pivot_wider(names_from = flora, values_from = spcount) 

# exp rich dat
d.expsprich <- readRDS("data/deficit_disp_latitude_data_2025.RDS") %>%
  select(c("entity_ID","disp.b_exp","disp.ab_exp")) %>%
  pivot_longer(cols = c(disp.b_exp, disp.ab_exp), names_to = "disptype", values_to = "exp_sp") %>%
  mutate(disptype = ifelse(disptype == "disp.b_exp", "disp.b", "disp.ab"))

# endemism data quality
d.check <- readRDS("data/disp_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# clean data: vars are scaled
d.debt <- readRDS("data/deficit_disp_latitude_data_2025.RDS") %>%
  select(-c("disp.b_exp", "disp.ab_exp", "propb_exp", "disp.b.diff", "disp.ab.diff", "sprichdiff", "Tdebt","C_disp.b.debt", "C_disp.ab.debt", "disp.b", "disp.ab")) %>%
  filter(entity_class2 == "Oceanic") %>%
  left_join(d.endnon, by ="entity_ID") %>%
  left_join(d.check, by = "entity_ID") %>%
  left_join(d.expsprich, by = c("entity_ID", "disptype")) %>%
  mutate(pend = endemic/(endemic + nonendemic)) %>%
  mutate(logitend = log(endemic + 0.01/(endemic+ 0.01 + nonendemic+ 0.01)/(nonendemic+ 0.01/(endemic+ 0.01 + nonendemic+ 0.01)))) %>%
  mutate(end_exists = as.integer(ifelse(endemic>0, 1,0))) %>%
  mutate(
    age_Ma = as.vector(scale(log10((age_Ma)+.01))),
    sc_nonendemic = as.vector(scale(log10((nonendemic)+.01))),
    area = as.vector(scale(log10((area)+.01))), 
    dist = as.vector(scale(dist)),
    temp = as.vector(scale(log10((temp)+.01))),
    prec = as.vector(scale(log10((prec)+.01))),
    elev_range = as.vector(scale(log10((elev_range)+.01))),
    abslatitude = as.vector(scale(log10((abslatitude)+.01)))) %>%
  drop_na() 

#########################
##### PRED END OCC ######
#### ONLY NONENDEMIC ####
#########################

table(dat$end_exists, dat$entity_ID)
# Some entity_ID have all same values of end_exists

dat <- d.debt
mod <- glm(end_exists ~ disptype + sc_nonendemic + age_Ma + dist + area + elev_range + #(1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:disptype + 
             disptype:age_Ma + disptype:dist + disptype:area + disptype:elev_range,
             family = binomial, weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/DispYN_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/DispYN_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model)

#########################
#### PRED END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

dat <- d.debt %>% filter(end_exists == 1) 
mod <- lmer(logitend ~ disptype + sc_nonendemic + age_Ma + dist + area + elev_range + (1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:disptype + 
             disptype:age_Ma + disptype:dist + disptype:area + disptype:elev_range,
             weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/DispDIVERS_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/DispDIVERS_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model) 

#########################
######### PLOT ##########
#########################

#########################
##### PLOT END OCC ######
#### ONLY NONENDEMIC ####
#########################

########## DISPTYPE ##########
ref <- lsmeans(best_model, pairwise ~ disptype, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)

plot <- 
  ggplot(ref.table, aes(reorder(disptype, lsmean), lsmean, color = disptype, fill = disptype)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(disptype, end_exists), position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Endemic occurence") + 
  xlab("Mutualist type") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  scale_x_discrete(labels = c("Abiotic", "Biotic")) +
  colFill.pt + colScale.pt 

png("figures/Disp_YN_ENDNonEnd_disptype.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## DISPTYPE*DIST ##########
dist.range <- seq(max(dat$dist),min(dat$dist),length.out=1000)
disptype.levels <- unique(dat$disptype)
ex.grid <- expand.grid(dist = dist.range, disptype = disptype.levels)

predictor_vars <- c("temp", "prec", "area", "pnata", "age_Ma", "elev_range", "sc_nonendemic")
pred.dat <- cbind(ex.grid, t(colMeans(na.omit(dat[, predictor_vars]), na.rm = TRUE)))
pred <- predict(best_model, newdata = pred.dat, se.fit = TRUE, type = "response", re.form = NA)

pred.nonend <- cbind(pred.dat, fit = pred$fit, se.fit = pred$se.fit) %>%
  mutate(conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

plot <- ggplot() +
  geom_line(data = pred.nonend, mapping = aes(x = dist, y = fit, color = disptype, fill = disptype), size = 1) +
  geom_ribbon(data = pred.nonend, aes(x = dist, ymin = conf.low, ymax = conf.high, color = disptype, fill = disptype),  alpha = 0.5) +
  geom_point(data = dat,aes(x = dist, y = end_exists, color = disptype, fill = disptype),size = 3,shape = 16, alpha = 0.5) +
  theme_classic(base_size = 40) +
  ylab("Endemic occurence") + 
  xlab("Distance") +
  theme(legend.position="none") +
  colFill + colScale

png("figures/Disp_YN_ENDNonEnd_disptypebydist.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

#########################
#### PLOT END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

####################################
####################################
################ MYC ###############
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "myctype", values = c("royalblue4","darkgrey"))
colFill <- scale_fill_manual(name = "myctype", values = c("royalblue4","darkgrey"))

colScale.pt <- scale_colour_manual(name = "myctype", values = c("darkgrey", "royalblue4"))
colFill.pt <- scale_fill_manual(name = "myctype", values = c("darkgrey", "royalblue4"))

# end/non dat
m.endnon <- readRDS("data/merge_myc_data_2025.RDS") %>%
  select(c("entity_ID","flora","AM","NM","age_Ma")) %>%
  filter(!flora == "native") %>%
  pivot_longer(cols = c(AM, NM), names_to = "myctype", values_to = "spcount") %>%
  pivot_wider(names_from = flora, values_from = spcount) 

# exp rich dat
m.expsprich <- readRDS("data/deficit_myc_latitude_data_2025.RDS") %>%
  select(c("entity_ID","AM_exp","NM_exp")) %>%
  pivot_longer(cols = c(AM_exp, NM_exp), names_to = "myctype", values_to = "exp_sp") %>%
  mutate(myctype = ifelse(myctype == "AM_exp", "AM", "NM"))

# endemisp data quality
m.check <- readRDS("data/myc_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# clean data: vars are scaled
m.debt <- readRDS("data/deficit_myc_latitude_data_2025.RDS") %>%
  select(-c("AM_exp", "NM_exp",  "ORC_exp", "EM_exp", "propAM_exp","propEM_exp", "propORC_exp", "AMdiff","EMdiff","ORCdiff", "NMdiff", "sprichdiff","EMdebt","ORCdebt", "Tdebt", "C_AMdebt", "C_EMdebt" ,"C_ORCdebt" ,"C_NMdebt")) %>%
  filter(entity_class2 == "Oceanic") %>%
  left_join(m.endnon, by ="entity_ID") %>%
  left_join(m.check, by = "entity_ID") %>%
  left_join(m.expsprich, by = c("entity_ID", "myctype")) %>%
  mutate(pend = endemic/(endemic + nonendemic)) %>%
  mutate(logitend = log(endemic + 0.01/(endemic+ 0.01 + nonendemic+ 0.01)/(nonendemic+ 0.01/(endemic+ 0.01 + nonendemic+ 0.01)))) %>%
  mutate(end_exists = as.integer(ifelse(endemic>0, 1,0))) %>%
  mutate(
    age_Ma = as.vector(scale(log10((age_Ma)+.01))),
    sc_nonendemic = as.vector(scale(log10((nonendemic)+.01))),
    area = as.vector(scale(log10((area)+.01))), 
    dist = as.vector(scale(dist)),
    temp = as.vector(scale(log10((temp)+.01))),
    prec = as.vector(scale(log10((prec)+.01))),
    elev_range = as.vector(scale(log10((elev_range)+.01))),
    abslatitude = as.vector(scale(log10((abslatitude)+.01)))) %>%
  drop_na()  

#########################
##### PRED END OCC ######
#### ONLY NONENDEMIC ####
#########################

table(dat$end_exists, dat$entity_ID)
# Some entity_ID have all same values of end_exists

dat <- m.debt
mod <- glm(end_exists ~ myctype + sc_nonendemic + age_Ma + dist + area + elev_range + #(1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:myctype + 
             myctype:age_Ma + myctype:dist + myctype:area + myctype:elev_range,
             family = binomial, weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Myc_YN_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Myc_YN_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model)

#########################
#### PRED END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

dat <- m.debt %>% filter(end_exists == 1) 
mod <- lmer(logitend ~ myctype + sc_nonendemic + age_Ma + dist + area + elev_range + (1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:myctype + 
             myctype:age_Ma + myctype:dist + myctype:area + myctype:elev_range,
             weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Myc_DIVERS_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Myc_DIVERS_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model) 

# find simple slope by muttype
emtrends(best_model, ~ myctype, var="sc_nonendemic")

#########################
######### PLOT ##########
#########################

#########################
##### PLOT END OCC ######
#### ONLY NONENDEMIC ####
#########################

#########################
#### PLOT END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

########## MYCTYPE ##########
ref <- emmeans(best_model, pairwise ~ myctype, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$emmeans)

plot <- 
  ggplot(ref.table, aes(reorder(myctype, emmean), emmean, color = myctype, fill = myctype)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(myctype, logitend), position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Endemic prevalence") + 
  xlab("Mutualist type") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  scale_x_discrete(labels = c("Abiotic", "Biotic")) +
  colFill + colScale +
  ylim(-0.5,3)

png("figures/Myc_DIVERS_ENDNonEnd_myctype.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## MYCTYPE*NONENDEMIC ##########
sc_nonendemic.range <- seq(max(dat$sc_nonendemic),min(dat$sc_nonendemic),length.out=1000)
myctype.levels <- unique(dat$myctype)
ex.grid <- expand.grid(sc_nonendemic = sc_nonendemic.range, myctype = myctype.levels)

predictor_vars <- c("temp", "prec", "area", "pnata", "age_Ma", "elev_range", "dist")
pred.dat <- cbind(ex.grid, t(colMeans(na.omit(dat[, predictor_vars]), na.rm = TRUE)))
pred <- predict(best_model, newdata = pred.dat, se.fit = TRUE, type = "response", re.form = NA)

pred.nonend <- cbind(pred.dat, fit = pred$fit, se.fit = pred$se.fit) %>%
  mutate(conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(pred.nonend$sc_nonendemic), max(pred.nonend$sc_nonendemic), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$nonendemic + 0.01), na.rm = TRUE) + mean(log10(dat$nonendemic + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = pred.nonend, mapping = aes(x = sc_nonendemic, y = fit, color = myctype, fill = myctype), size = 1) +
  geom_ribbon(data = pred.nonend, aes(x = sc_nonendemic, ymin = conf.low, ymax = conf.high, color = myctype, fill = myctype),  alpha = 0.5) +
  geom_point(data = dat,aes(x = sc_nonendemic, y = logitend, color = myctype, fill = myctype),size = 3,shape = 16, alpha = 0.5) +
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") + 
  xlab("Colonist richness") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position="none") +
  colFill + colScale +
  ylim(-2.5,6)

png("figures/Myc_DIVERS_ENDNonEnd_myctypebysc_nonendemic.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

####################################
####################################
################ NFIX ##############
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "nfixtype", values = c("darkseagreen3", "darkgrey"))
colFill <- scale_fill_manual(name = "nfixtype", values = c("darkseagreen3", "darkgrey"))

# end/non dat
n.endnon <- readRDS("data/merge_nfix_data_2025.RDS") %>%
  select(c("entity_ID","flora","nfix","nonfix","age_Ma")) %>%
  filter(!flora == "native") %>%
  pivot_longer(cols = c(nfix,nonfix), names_to = "nfixtype", values_to = "spcount") %>%
  pivot_wider(names_from = flora, values_from = spcount) 

# exp rich dat
n.expsprich <- readRDS("data/deficit_nfix_latitude_data_2025.RDS") %>%
  select(c("entity_ID","nfix_exp","nonfix_exp")) %>%
  pivot_longer(cols = c(nfix_exp, nonfix_exp), names_to = "nfixtype", values_to = "exp_sp") %>%
  mutate(nfixtype = ifelse(nfixtype == "nfix_exp", "nfix", "nonfix"))

# endemisp data quality
n.check <- readRDS("data/nfix_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# vars are scaled.
n.debt <- readRDS("data/deficit_nfix_latitude_data_2025.RDS") %>%
  select(-c("nfix_exp", "nonfix_exp", "propnfix_exp", "nfix.diff", "nonfix.diff", "sprichdiff", "Tdebt","C_nfix.debt", "C_nonfix.debt", "nfix", "nonfix")) %>%
  filter(entity_class2 == "Oceanic") %>%
  left_join(n.endnon, by ="entity_ID") %>%
  left_join(n.check, by = "entity_ID") %>%
  left_join(n.expsprich, by = c("entity_ID", "nfixtype")) %>%
  mutate(pend = endemic/(endemic + nonendemic)) %>%
  mutate(logitend = log(endemic + 0.01/(endemic+ 0.01 + nonendemic+ 0.01)/(nonendemic+ 0.01/(endemic+ 0.01 + nonendemic+ 0.01)))) %>%
  mutate(end_exists = as.integer(ifelse(endemic>0, 1,0))) %>%
  mutate(
    age_Ma = as.vector(scale(log10((age_Ma)+.01))),
    sc_nonendemic = as.vector(scale(log10((nonendemic)+.01))),
    area = as.vector(scale(log10((area)+.01))), 
    dist = as.vector(scale(dist)),
    temp = as.vector(scale(log10((temp)+.01))),
    prec = as.vector(scale(log10((prec)+.01))),
    elev_range = as.vector(scale(log10((elev_range)+.01))),
    abslatitude = as.vector(scale(log10((abslatitude)+.01)))) %>%
  drop_na() 

#########################
##### PRED END OCC ######
#### ONLY NONENDEMIC ####
#########################

table(dat$end_exists, dat$entity_ID)
# Some entity_ID have all same values of end_exists

dat <- n.debt
mod <- glm(end_exists ~ nfixtype + sc_nonendemic + age_Ma + dist + area + elev_range + #(1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:nfixtype + 
             nfixtype:age_Ma + nfixtype:dist + nfixtype:area + nfixtype:elev_range,
             family = binomial, weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Nfix_YN_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Nfix_YN_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model)

#########################
#### PRED END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

dat <- n.debt %>% filter(end_exists == 1) 
mod <- lmer(logitend ~ nfixtype + sc_nonendemic + age_Ma + dist + area + elev_range + (1|entity_ID) +
             sc_nonendemic:age_Ma + sc_nonendemic:dist + sc_nonendemic:area + sc_nonendemic:elev_range+ sc_nonendemic:nfixtype + 
             nfixtype:age_Ma + nfixtype:dist + nfixtype:area + nfixtype:elev_range,
             weights = pnata, data = dat)
summary(mod)

dredged_models <- dredge(mod, rank = AICc)
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)
saveRDS(averaged_model_results, "models/Nfix_DIVERS_ENDNonEnd_averaged_model_results.RDS")
averaged_model_results <- readRDS("models/Nfix_DIVERS_ENDNonEnd_averaged_model_results.RDS")
summary(averaged_model_results) 

best_model <- get.models(averaged_model_results, subset = 1)[[1]]
summary(best_model) 

# find simple slope by muttype
emtrends(best_model, ~ nfixtype, var="sc_nonendemic")

#########################
######### PLOT ##########
#########################

#########################
##### PLOT END OCC ######
#### ONLY NONENDEMIC ####
#########################

#########################
#### PLOT END DIVERS ####
#### ONLY NONENDEMIC ####
#########################

########## NFIXTYPE ##########
ref <- emmeans(best_model, pairwise ~ nfixtype, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$emmeans)

plot <- 
  ggplot(ref.table, aes(reorder(nfixtype, emmean), emmean, color = nfixtype, fill = nfixtype)) +  
  geom_bar(stat="identity", alpha = 0.2) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(nfixtype, logitend), position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Endemic prevalence") + 
  xlab("Mutualist type") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  scale_x_discrete(labels = c("Abiotic", "Biotic")) +
  colFill + colScale +
  ylim(-0.5,3)

png("figures/Nfix_DIVERS_ENDNonEnd_nfixtype.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

########## NFIXTYPE*NONENDEMIC ##########
sc_nonendemic.range <- seq(max(dat$sc_nonendemic),min(dat$sc_nonendemic),length.out=1000)
nfixtype.levels <- unique(dat$nfixtype)
ex.grid <- expand.grid(sc_nonendemic = sc_nonendemic.range, nfixtype = nfixtype.levels)

predictor_vars <- c("temp", "prec", "area", "pnata", "age_Ma", "elev_range", "dist")
pred.dat <- cbind(ex.grid, t(colMeans(na.omit(dat[, predictor_vars]), na.rm = TRUE)))
pred <- predict(best_model, newdata = pred.dat, se.fit = TRUE, type = "response", re.form = NA)

pred.nonend <- cbind(pred.dat, fit = pred$fit, se.fit = pred$se.fit) %>%
  mutate(conf.low = fit - 1.96 * se.fit, conf.high = fit + 1.96 * se.fit)

x_breaks <- seq(min(pred.nonend$sc_nonendemic), max(pred.nonend$sc_nonendemic), length.out = 5)
x_labels <- 10^(x_breaks * sd(log10(dat$nonendemic + 0.01), na.rm = TRUE) + mean(log10(dat$nonendemic + 0.01), na.rm = TRUE)) - 0.01

plot <- ggplot() +
  geom_line(data = pred.nonend, mapping = aes(x = sc_nonendemic, y = fit, color = nfixtype, fill = nfixtype), size = 1) +
  geom_ribbon(data = pred.nonend, aes(x = sc_nonendemic, ymin = conf.low, ymax = conf.high, color = nfixtype, fill = nfixtype),  alpha = 0.5) +
  geom_point(data = dat,aes(x = sc_nonendemic, y = logitend, color = nfixtype, fill = nfixtype),size = 3,shape = 16, alpha = 0.5) +
  theme_classic(base_size = 40) +
  ylab("Endemic prevalence") + 
  xlab("Colonist richness") +
  scale_x_continuous(breaks = x_breaks, labels = round(x_labels, 0)) +
  theme(legend.position="none") +
  colFill + colScale +
  ylim(-2.5,6)

png("figures/Nfix_DIVERS_ENDNonEnd_nfixypebysc_nonendemic.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

