############################
###### LOAD PACKAGES #######
############################

library(mgcv); library(gridExtra); library(betareg); library(MASS); library(lme4); library(lmerTest); library(lsmeans); library(ggeffects); library(spdep); library(ggplot2); library(ncf); library(ape); library(sjPlot); library(gridExtra); library(MuMIn);library(maps); library(sf); library(car);library(viridis);library(tidyverse)

############################
###### LOAD FUNCTIONS ######
############################

source("Master_modelling.R")

############################
######## READ DATA #########
############################

dat <-  readRDS("data/endnonend_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","endemic", "nonendemic", "latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland")) %>%
  mutate(entity_class2 = ifelse(is.na(geology) & entity_class == "Island", "Island", entity_class2)) %>%
  select(-geology) %>%          
  filter(area > 6) %>% 
  mutate(abslatitude = as.vector(abs(latitude))) %>%
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>%
  drop_na()

dat2 <-  readRDS("data/endnonend_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","endemic", "nonendemic","latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range","dist")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland")) %>%
  mutate(entity_class2 = ifelse(is.na(geology) & entity_class == "Island", "Island", entity_class2)) %>%
  select(-geology) %>%          
  filter(area > 6) %>%
  mutate(abslatitude = as.vector(abs(latitude))) %>%
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>%
  drop_na()

############################
######## ML PREDICT ########
############################

############################
######### MAINLAND #########
############################

dat.ml <- dat %>% filter(entity_class2 == "Mainland") 

gam.mod_sprich <- gam(sprich ~ s(abslatitude) , family = nb(link = "log"), data = dat.ml) 
summary(gam.mod_sprich)

mod <- gam.mod_sprich
new.dat.sprich <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred_sprich <- predict.gam(mod,newdata = new.dat.sprich, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.sprich$abslatitude)

# check assumptions
gam.check(gam.mod_sprich)
disp_check(gam.mod_sprich,dat)

############################
###### PREDICT EXP IS ######
######## ALL ISLANDS #######
############################

dat.is.min <- dat %>% 
  filter(!entity_class2 == "Mainland") %>%
  select(c('entity_ID',"abslatitude"))

pred_sprich <- predict.gam(gam.mod_sprich, newdata = dat.is.min, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_sprich_df <- cbind(dat.is.min,pred_sprich) %>% rename(sprich_exp = fit)

pred_is.dat <- pred_sprich_df %>%
  left_join(dat2, by= c('entity_ID','abslatitude')) %>%
  #keep only where there is actual end/non data
  filter(!(endemic == 0 & nonendemic == 0)) %>%
  #make sure nonendemic
  mutate(sprichdiff = sprich_exp - nonendemic) %>%
  mutate(Tdebt = (sprichdiff/sprich_exp)) %>%
  mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt))

############################
######## WRITE DATA ########
############################

# debt
saveRDS(pred_is.dat,"data/deficit_latitude_data_2025.RDS")
