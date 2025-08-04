# This script:
# 1. calculates island species deficit for nfixing status in nonendemics

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

dat <-  readRDS("data/merge_nfix_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area",  "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range",
           "nfix", "nonfix")) %>%
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

dat2 <-  readRDS("data/merge_nfix_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "elev_range", "dist", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec", 
           "nfix", "nonfix")) %>%
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

######################################
############ MAINLAND NFIX ###########
######################################

gam.mod.nfix <- gam(nfix ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.nfix)

mod <- gam.mod.nfix
new.dat.nfix <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000)))
pred.nfix <- predict.gam(mod,newdata = new.dat.nfix, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.nfix$abslatitude, nfix = "nfix") 

# check assumptions
gam.check(gam.mod.nfix)
disp_check(gam.mod.nfix,dat)

###############################
####### MAINLAND NONFIX #######
###############################

gam.mod.nonfix <- gam(nonfix ~ s(abslatitude) , family=nb(link="log"), data = dat.ml) 
summary(gam.mod.nonfix)

mod <- gam.mod.nonfix
new.dat.nonfix <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000)))
pred.nonfix <- predict.gam(mod, newdata = new.dat.nonfix, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.nonfix$abslatitude, nfix ="nonfix") 

# check assumptions
gam.check(gam.mod.nonfix)
disp_check(gam.mod.nonfix,dat)

############################
########## PLOT ############
##### LAT BY NFIX TYPE #####
############################

pred.mainland <- rbind(pred.nfix,pred.nonfix)

dat.ml.cond <- dat.ml %>%
  select(abslatitude, nfix, nonfix) %>%
  gather(key = "nfix", value = "sprich", nfix, nonfix)

# create a custom color scale
colScale <- scale_colour_manual(values=c("darkseagreen3", "darkgrey"))
fillScale <- scale_fill_manual(values=c("darkseagreen3", "darkgrey"))

pred.mainland$nfix <- ordered(pred.mainland$nfix, levels = c("nfix", "nonfix"))
dat.ml.cond$nfix <- ordered(dat.ml.cond$nfix, levels = c("nfix", "nonfix"))

# rename nfix for legend:

pred.mainland <- pred.mainland %>% mutate(nfix = case_when(nfix=="nfix" ~ "N-fixing",
                                                           nfix=="nonfix" ~ "Non N-fixing"))

dat.ml.cond <- dat.ml.cond %>% mutate(nfix = case_when(nfix=="nfix" ~ "N-fixing",
                                                       nfix=="nonfix" ~ "Non N-fixing"))

lat.nfixtype <-
  ggplot(pred.mainland, aes(x = abslatitude, y = fit, color = nfix, fill = nfix))+
  geom_line(size = 1) +
  geom_point(data = dat.ml.cond, aes(x = abslatitude, y = sprich, color = factor(nfix)), alpha = 0.3, size = 4)+ 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.5) + 
  xlab("Absolute latitude") +
  ylab("Species richness")+
  theme_classic(base_size = 25)+
  ylim(0,700)+
  colScale+
  fillScale+
  theme(legend.position = 'none')+
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  guides(color = guide_legend(title="N-fixing \nType",override.aes=list(fill = NA,size = 3, alpha = 0.7,linetype = c(0, 0))))+
  theme(axis.text.x = element_text(size =20),axis.text.y = element_text(angle = 45,size=20))

# write out
png("figures/nfix_byabslat.jpg", width=10, height= 10, units='in', res=300)
lat.nfixtype
dev.off()

############################
###### PREDICT EXP IS ######
######## ALL ISLANDS #######
############################

dat.is.min <- dat %>% 
  filter(!entity_class2 == "Mainland") %>%
  select(c('entity_ID',"abslatitude"))

pred_sprich <- predict.gam(gam.mod_sprich,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_sprich_df <- cbind(dat.is.min,pred_sprich) %>% rename(sprich_exp = fit)

pred_nfix <- predict.gam(gam.mod.nfix,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_nfix_df <- cbind(dat.is.min,pred_nfix) %>% rename(nfix_exp = fit)

pred_nonfix <- predict.gam(gam.mod.nonfix,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_nonfix_df <- cbind(dat.is.min,pred_nonfix) %>% rename(nonfix_exp = fit)

pred_is.dat <- pred_is.dat <- pred_sprich_df %>%
  left_join(dat2, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_nfix_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_nonfix_df, by= c('entity_ID','abslatitude')) %>%
  mutate_at(c('nfix_exp','nonfix_exp','sprich_exp'), as.integer) %>%
  mutate(propnfix_exp = nfix_exp/(nfix_exp + nonfix_exp))

saveRDS(pred_is.dat,"data/exp_nfix_latitude_data_2025.RDS")

############################
######## CREATE DATA #######
############################

pred_is.dat.shrink <- pred_is.dat %>%
  mutate(sprich = nfix + nonfix, sprich_exp = nfix_exp + nonfix_exp) %>%
  mutate(nfix.diff = nfix_exp - nfix, nonfix.diff = nonfix_exp - nonfix, sprichdiff = sprich_exp - sprich) %>%
  mutate(nfix.debt = (nfix.diff/nfix_exp), nonfix.debt = (nonfix.diff/nonfix_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(nfix.debt = ifelse(nfix.debt < 0, 0, nfix.debt)) %>% mutate(nonfix.debt = ifelse(nonfix.debt <0, 0, nonfix.debt)) %>% mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(nfix.debt = ifelse(nfix.debt > 1, 1, nfix.debt)) %>% mutate(nonfix.debt = ifelse(nonfix.debt >1, 1, nonfix.debt)) %>% mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt)) %>%
  mutate(C_nfix.debt = (nfix.diff/sprichdiff), C_nonfix.debt = (nonfix.diff/sprichdiff)) %>% 
  mutate(C_nfix.debt = ifelse(C_nfix.debt < 0, 0, C_nfix.debt)) %>% mutate(C_nonfix.debt = ifelse(C_nonfix.debt <0, 0, C_nonfix.debt)) %>% 
  mutate(C_nfix.debt = ifelse(C_nfix.debt > 1, 1, C_nfix.debt)) %>% mutate(C_nonfix.debt = ifelse(C_nonfix.debt >1, 1, C_nonfix.debt)) 

############################
######## CDEBT DATA ########
############################

pred_is.dat.alld <- pred_is.dat.shrink %>% 
  select(c('entity_ID','nfix.debt','nonfix.debt')) %>%
  gather(key = "nfix", value = "debt", nfix.debt, nonfix.debt) %>%
  mutate(nfix = case_when(nfix == "nfix.debt" ~ "nfix",                                   
                          nfix == "nonfix.debt" ~ "nonfix")) 

pred_is.dat.allcd <- pred_is.dat.shrink %>% 
  select(c('entity_ID','C_nfix.debt','C_nonfix.debt')) %>%
  gather(key = "nfix", value = "debt.c", C_nfix.debt, C_nonfix.debt)%>%
  mutate(nfix = case_when(nfix == "C_nfix.debt" ~ "nfix",                                   
                          nfix == "C_nonfix.debt" ~ "nonfix"))

pred_is.dat.all.diff <- pred_is.dat.shrink %>% 
  select(c('entity_ID','entity_class2', 'sprich','latitude','longitude','abslatitude','nfix.diff','nonfix.diff','dist','area','elev_range', 'temp','prec')) %>%
  gather(key = "nfix", value = "diff", nfix.diff, nonfix.diff) %>%
  mutate(nfix = case_when(nfix == "nfix.diff" ~ "nfix",                                   
                          nfix == "nonfix.diff" ~ "nonfix")) 

pred_is.dat.all.exp <- pred_is.dat.shrink %>% 
  select(c('entity_ID','nfix_exp','nonfix_exp')) %>%
  gather(key = "nfix", value = "exp", nfix_exp, nonfix_exp) %>%
  mutate(nfix = case_when(nfix == "nfix_exp" ~ "nfix",                                   
                          nfix == "nonfix_exp" ~ "nonfix")) 

pred_is.dat.all.obs <- pred_is.dat.shrink %>% 
  select(c('entity_ID', 'nfix', 'nonfix')) %>%
  gather(key = "nfix", value = "obs", nfix, nonfix) 

pred_is.dat.all.T <- pred_is.dat.shrink %>% select(entity_ID, sprichdiff)

pred_is.dat.all <- pred_is.dat.all.diff %>%
  left_join(pred_is.dat.all.exp, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.all.obs, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.alld, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.allcd, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.all.T, by = c("entity_ID")) %>%
  mutate(debt.weights = exp, debt.c.weights = abs(sprichdiff)) %>%
  mutate(nfix = as.factor(nfix))

pred_is.dat.all <- within(pred_is.dat.all, nfix <- relevel(nfix, ref = "nonfix"))

############################
######## WRITE DATA ########
############################

#debt
saveRDS(pred_is.dat.shrink,"data/deficit_nfix_latitude_data_2025.RDS")

#cdebt
saveRDS(pred_is.dat.all,"data/cdeficit_nfix_latitude_data_2025.RDS")
