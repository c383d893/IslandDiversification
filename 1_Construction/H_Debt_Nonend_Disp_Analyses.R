# This script:
# 1. calculates island species deficit for dispersal syndrome in nonendemics

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

dat <-  readRDS("data/merge_disp_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area",  "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range",
           "disp.b", "disp.ab")) %>%
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

dat2 <-  readRDS("data/merge_disp_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "elev_range", "dist", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec", 
           "disp.b", "disp.ab")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland")) %>%
  mutate(entity_class2 = ifelse(is.na(geology) & entity_class == "Island", "Island", entity_class2)) %>%
  select(-geology) %>%          
  filter(area > 6) %>% # based on paper Patrick shared
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>% 
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>% 
  drop_na()

############################
######### MAINLAND #########
############################

dat.ml <- dat %>% filter(entity_class2=="Mainland")

gam.mod_sprich <- gam(sprich ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod_sprich)

mod <- gam.mod_sprich
new.dat.sprich <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred_sprich <- predict.gam(mod,newdata = new.dat.sprich, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.sprich$abslatitude) 

# check assumptions
gam.check(gam.mod_sprich)
disp_check(gam.mod_sprich,dat)

#####################################
####### MAINLAND BIOTIC DISP ########
#####################################

gam.mod.disp.b <- gam(disp.b ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.disp.b)

mod <- gam.mod.disp.b
new.dat.disp.b <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.disp.b <- predict.gam(mod,newdata = new.dat.disp.b, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.disp.b$abslatitude, disp = "disp.b") 

# check assumptions
gam.check(gam.mod.disp.b)
disp_check(gam.mod.disp.b,dat)

#####################################
####### MAINLAND ABIOTIC DISP #######
#####################################

gam.mod.disp.ab <- gam(disp.ab ~ s(abslatitude) , family=nb(link="log"), data = dat.ml) 
summary(gam.mod.disp.ab)

mod <- gam.mod.disp.ab
new.dat.disp.ab <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.disp.ab <- predict.gam(mod, newdata = new.dat.disp.ab, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.disp.ab$abslatitude, disp = "disp.ab") 

# check assumptions
gam.check(gam.mod.disp.ab)
disp_check(gam.mod.disp.ab,dat)

############################
########## PLOT ############
##### LAT BY DISP TYPE #####
############################

pred.mainland <- rbind(pred.disp.b,pred.disp.ab)

dat.ml.cond <- dat.ml %>%
  select(abslatitude, disp.b, disp.ab) %>%
  gather(key = "disp", value = "sprich", disp.b, disp.ab)

# create a custom color scale
colScale <- scale_colour_manual(values=c("darkgrey", "brown4"))
fillScale <- scale_fill_manual(values=c("darkgrey","brown4"))

pred.mainland$disp <- ordered(pred.mainland$disp, levels = c("disp.b", "disp.ab"))
dat.ml.cond$disp <- ordered(dat.ml.cond$disp, levels = c("disp.b", "disp.ab"))

# rename disp for legend:
pred.mainland <- pred.mainland %>% mutate(disp = case_when(disp=="disp.b" ~ "Biotic",
                                                           disp=="disp.ab" ~ "Abiotic"))

dat.ml.cond <- dat.ml.cond %>% mutate(disp = case_when(disp=="disp.b" ~ "Biotic",
                                                       disp=="disp.ab" ~ "Abiotic"))
lat.disptype <-
  ggplot(pred.mainland, aes(x = abslatitude, y = fit,color = disp, fill = disp))+
  geom_line(size = 1) +
  geom_point(data = dat.ml.cond, aes(x = abslatitude, y = sprich, color = factor(disp)), alpha = 0.3, size = 4)+ 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.5) + 
  xlab("Absolute latitude") +
  ylab("Species richness")+
  theme_classic(base_size = 25)+
  ylim(0,3000)+
  colScale+
  fillScale+
  theme(legend.position = 'none')+
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  guides(color = guide_legend(title="dispination \nSyndrome",override.aes=list(fill =NA,size =3, alpha =0.7,linetype = c(0, 0))))+
  theme(axis.text.x = element_text(size =20),axis.text.y = element_text(angle = 45,size=20))

# write out
png("figures/disp_byabslat.jpg", width=10, height= 10, units='in', res=300)
lat.disptype
dev.off()

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

pred_disp.b <- predict.gam(gam.mod.disp.b,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_disp.b_df <- cbind(dat.is.min,pred_disp.b) %>% rename(disp.b_exp = fit)

pred_disp.ab <- predict.gam(gam.mod.disp.ab,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_disp.ab_df <- cbind(dat.is.min,pred_disp.ab) %>% rename(disp.ab_exp = fit)

pred_is.dat <- pred_sprich_df %>%
  left_join(dat2, by= c('entity_ID','abslatitude')) %>%
  mutate(sprichdiff = sprich_exp - sprich) %>%
  mutate(Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt))

pred_is.dat <- pred_sprich_df %>%
  left_join(dat2, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_disp.b_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_disp.ab_df, by= c('entity_ID','abslatitude')) %>%
  mutate_at(c('disp.b_exp','disp.ab_exp','sprich_exp'), as.integer) %>%
  mutate(propb_exp = disp.b_exp/(disp.b_exp + disp.ab_exp))

saveRDS(pred_is.dat,"data/exp_disp_latitude_data_2025.RDS")

############################
######## CREATE DATA #######
############################

pred_is.dat.shrink <- pred_is.dat %>%
  mutate(sprich = disp.b + disp.ab, sprich_exp = disp.b_exp + disp.ab_exp) %>%
  mutate(disp.b.diff = disp.b_exp - disp.b, disp.ab.diff = disp.ab_exp - disp.ab, sprichdiff = sprich_exp - sprich) %>%
  mutate(disp.b.debt = (disp.b.diff/disp.b_exp), disp.ab.debt = (disp.ab.diff/disp.ab_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(disp.b.debt = ifelse(disp.b.debt < 0, 0, disp.b.debt)) %>% mutate(disp.ab.debt = ifelse(disp.ab.debt <0, 0, disp.ab.debt)) %>% mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(disp.b.debt = ifelse(disp.b.debt > 1, 1, disp.b.debt)) %>% mutate(disp.ab.debt = ifelse(disp.ab.debt >1, 1, disp.ab.debt)) %>% mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt)) %>%
  mutate(C_disp.b.debt = (disp.b.diff/sprichdiff), C_disp.ab.debt = (disp.ab.diff/sprichdiff)) %>% 
  mutate(C_disp.b.debt = ifelse(C_disp.b.debt < 0, 0, C_disp.b.debt)) %>% mutate(C_disp.ab.debt = ifelse(C_disp.ab.debt <0, 0, C_disp.ab.debt)) %>% 
  mutate(C_disp.b.debt = ifelse(C_disp.b.debt > 1, 1, C_disp.b.debt)) %>% mutate(C_disp.ab.debt = ifelse(C_disp.ab.debt >1, 1, C_disp.ab.debt)) 

############################
######## CDEBT DATA ########
############################

pred_is.dat.alld <- pred_is.dat.shrink %>% 
  select(c('entity_ID','disp.b.debt','disp.ab.debt')) %>%
  gather(key = "disp", value = "debt", disp.b.debt, disp.ab.debt) %>%
  mutate(disp = case_when(disp == "disp.b.debt" ~ "disp.b",                                   
                          disp == "disp.ab.debt" ~ "disp.ab")) 

pred_is.dat.allcd <- pred_is.dat.shrink %>% 
  select(c('entity_ID','C_disp.b.debt','C_disp.ab.debt')) %>%
  gather(key = "disp", value = "debt.c", C_disp.b.debt, C_disp.ab.debt)%>%
  mutate(disp = case_when(disp == "C_disp.b.debt" ~ "disp.b",                                   
                          disp == "C_disp.ab.debt" ~ "disp.ab"))

pred_is.dat.all.diff <- pred_is.dat.shrink %>% 
  select(c('entity_ID','entity_class2', 'sprich','latitude','longitude','abslatitude','disp.b.diff','disp.ab.diff','dist','area','elev_range', 'temp','prec')) %>%
  gather(key = "disp", value = "diff", disp.b.diff, disp.ab.diff) %>%
  mutate(disp = case_when(disp == "disp.b.diff" ~ "disp.b",                                   
                          disp == "disp.ab.diff" ~ "disp.ab")) 

pred_is.dat.all.exp <- pred_is.dat.shrink %>% 
  select(c('entity_ID','disp.b_exp','disp.ab_exp')) %>%
  gather(key = "disp", value = "exp", disp.b_exp, disp.ab_exp) %>%
  mutate(disp = case_when(disp == "disp.b_exp" ~ "disp.b",                                   
                          disp == "disp.ab_exp" ~ "disp.ab")) 

pred_is.dat.all.obs <- pred_is.dat.shrink %>% 
  select(c('entity_ID', 'disp.b', 'disp.ab')) %>%
  gather(key = "disp", value = "obs", disp.b, disp.ab) 

pred_is.dat.all.T <- pred_is.dat.shrink %>% select(entity_ID, sprichdiff)

pred_is.dat.all <- pred_is.dat.all.diff %>%
  left_join(pred_is.dat.all.exp, by = c("disp", "entity_ID")) %>%
  left_join(pred_is.dat.all.obs, by = c("disp", "entity_ID")) %>%
  left_join(pred_is.dat.alld, by = c("disp", "entity_ID")) %>%
  left_join(pred_is.dat.allcd, by = c("disp", "entity_ID")) %>%
  left_join(pred_is.dat.all.T, by = c("entity_ID")) %>%
  mutate(debt.weights = exp, debt.c.weights = abs(sprichdiff)) %>%
  mutate(disp = as.factor(disp)) 

pred_is.dat.all <- within(pred_is.dat.all, disp <- relevel(disp, ref = "disp.ab"))

############################
######## WRITE DATA ########
############################

#debt
saveRDS(pred_is.dat.shrink,"data/deficit_disp_latitude_data_2025.RDS")

#cdebt
saveRDS(pred_is.dat.all,"data/cdeficit_disp_latitude_data_2025.RDS")
