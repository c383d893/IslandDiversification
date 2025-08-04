# This script:
# 1. calculates island species deficit for myc type in nonendemics

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

dat <-  readRDS("data/merge_myc_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range",
           "AM","EM","ORC","NM")) %>%
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

dat2 <-  readRDS("data/merge_myc_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "elev_range", "dist", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec", 
           "AM","EM","ORC","NM")) %>%
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
####### MAINLAND AM ########
############################

gam.mod.AM <- gam(AM ~ s(abslatitude),family=nb(link="log"), data = dat.ml) 
summary(gam.mod.AM)

mod <- gam.mod.AM
new.dat.AM <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.AM <- predict.gam(mod,newdata = new.dat.AM, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.AM$abslatitude, myctype = "AM") 

# check assumptions
gam.check(gam.mod.AM)
disp_check(gam.mod.AM,dat)

############################
####### MAINLAND EM ########
############################

gam.mod.EM <- gam(EM ~ s(abslatitude),family=nb(link="log"), data = dat.ml) 
summary(gam.mod.EM)

mod <- gam.mod.EM
new.dat.EM <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.EM <- predict.gam(mod,newdata = new.dat.EM, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.EM$abslatitude, myctype = "EM") 

# check assumptions
gam.check(gam.mod.EM)
disp_check(gam.mod.EM,dat)

############################
####### MAINLAND ORC #######
############################

gam.mod.ORC <- gam(ORC ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.ORC)

mod <- gam.mod.ORC
new.dat.ORC <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.ORC <- predict.gam(mod,newdata = new.dat.ORC, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.ORC$abslatitude, myctype = "ORC") 

# check assumptions
gam.check(gam.mod.ORC)
disp_check(gam.mod.ORC,dat)

############################
####### MAINLAND NM ########
############################

gam.mod.NM <- gam(NM ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.NM)

mod <- gam.mod.NM
new.dat.NM <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.NM <- predict.gam(mod,newdata = new.dat.NM, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.NM$abslatitude, myctype = "NM") 

# check assumptions
gam.check(gam.mod.NM)
disp_check(gam.mod.NM,dat)

############################
########## PLOT ############
##### LAT BY MYC TYPE ######
############################

#pred.mainland <- rbind(pred.AM,pred.EM,pred.ORC,pred.NM)
pred.mainland <- rbind(pred.AM,pred.EM,pred.NM)

dat.ml.cond <- dat.ml %>%
  select(abslatitude, AM,EM,NM) %>%
  gather(key="myctype", value="sprich", AM,EM,NM) 

# create a custom color scale
colScale <- scale_colour_manual(values=c("royalblue4", "royalblue2","darkgrey"))
fillScale <- scale_fill_manual(values=c("royalblue4", "royalblue2","darkgrey"))

pred.mainland$myctype <- ordered(pred.mainland$myctype, levels = c("AM","EM","NM"))
dat.ml.cond$myctype <- ordered(dat.ml.cond$myctype, levels = c("AM","EM","NM"))

lat.myctype <-
ggplot(pred.mainland,aes(x =abslatitude,y=fit,color =myctype, fill =myctype))+
  geom_line(size=1) +
  geom_point(data=dat.ml.cond, aes(x = abslatitude,y=sprich, color =factor(myctype)), alpha=0.3,size=4)+ 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.5) + 
  xlab("Absolute latitude") +
  ylab("Species richness")+
  theme_classic(base_size = 25)+
  ylim(0,4500)+
  colScale+
  fillScale+
  theme(legend.position = 'none')+
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  guides(color = guide_legend(title="Mycorrhizal \nType",override.aes=list(fill =NA,size =3, alpha =0.7,linetype = c(0, 0, 0))))+
  theme(axis.text.x = element_text(size =20),axis.text.y = element_text(angle = 45,size=20))

# write out
png("figures/myc_byabslat.jpg", width=10, height= 10, units='in', res=300)
lat.myctype
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

pred_AM <- predict.gam(gam.mod.AM,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_AM_df <- cbind(dat.is.min,pred_AM) %>% rename(AM_exp = fit)

pred_EM <- predict.gam(gam.mod.EM,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_EM_df <- cbind(dat.is.min,pred_EM) %>% rename(EM_exp = fit)

pred_ORC <- predict.gam(gam.mod.ORC,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_ORC_df <- cbind(dat.is.min,pred_ORC) %>% rename(ORC_exp = fit)

pred_NM <- predict.gam(gam.mod.NM,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_NM_df <- cbind(dat.is.min,pred_NM) %>% rename(NM_exp = fit)

pred_is.dat <- pred_sprich_df %>%
  left_join(dat2, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_AM_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_EM_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_ORC_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_NM_df, by= c('entity_ID','abslatitude')) %>%
  mutate_at(c('AM_exp','EM_exp','ORC_exp', 'NM_exp','sprich_exp'), as.integer) %>%
  mutate(propAM_exp = AM_exp/(AM_exp + EM_exp + ORC_exp + NM_exp)) %>%
  mutate(propEM_exp = EM_exp/(AM_exp + EM_exp + ORC_exp + NM_exp)) %>%
  mutate(propORC_exp = ORC_exp/(AM_exp + EM_exp + ORC_exp + NM_exp))

# write out
saveRDS(pred_is.dat,"data/exp_myc_latitude_data_2025.RDS")

############################
######## CREATE DATA #######
############################

pred_is.dat.shrink <- pred_is.dat %>%
  mutate(sprich = AM + EM + ORC + NM, sprich_exp = AM_exp + EM_exp + ORC_exp + NM_exp) %>%
  mutate(AMdiff = AM_exp - AM, EMdiff = EM_exp - EM, ORCdiff = ORC_exp - ORC, NMdiff = NM_exp - NM, sprichdiff = sprich_exp - sprich) %>%
  mutate(AMdebt = (AMdiff/AM_exp), EMdebt = (EMdiff/EM_exp), ORCdebt = (ORCdiff/ORC_exp), NMdebt = (NMdiff/NM_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(AMdebt = ifelse(AMdebt < 0, 0, AMdebt)) %>% mutate(EMdebt = ifelse(EMdebt <0, 0, EMdebt)) %>% mutate(ORCdebt = ifelse(ORCdebt <0, 0, ORCdebt)) %>% mutate(NMdebt = ifelse(NMdebt <0, 0, NMdebt)) %>% mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(AMdebt = ifelse(AMdebt > 1, 1, AMdebt)) %>% mutate(EMdebt = ifelse(EMdebt >1, 1, EMdebt)) %>% mutate(ORCdebt = ifelse(ORCdebt >1, 1, ORCdebt)) %>% mutate(NMdebt = ifelse(NMdebt >1, 1, NMdebt)) %>% mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt)) %>%
  mutate(C_AMdebt = (AMdiff/sprichdiff), C_EMdebt = (EMdiff/sprichdiff), C_ORCdebt = (ORCdiff/sprichdiff), C_NMdebt = (NMdiff/sprichdiff)) %>% 
  mutate(C_AMdebt = ifelse(C_AMdebt < 0, 0, C_AMdebt)) %>% mutate(C_EMdebt = ifelse(C_EMdebt <0, 0, C_EMdebt)) %>% mutate(C_ORCdebt = ifelse(C_ORCdebt <0, 0, C_ORCdebt)) %>% mutate(C_NMdebt = ifelse(C_NMdebt <0, 0, C_NMdebt)) %>%
  mutate(C_AMdebt = ifelse(C_AMdebt > 1, 1, C_AMdebt)) %>% mutate(C_EMdebt = ifelse(C_EMdebt >1, 1, C_EMdebt)) %>% mutate(C_ORCdebt = ifelse(C_ORCdebt >1, 1, C_ORCdebt)) %>% mutate(C_NMdebt = ifelse(C_NMdebt >1, 1, C_NMdebt)) 

############################
######## CDEBT DATA ########
############################

pred_is.dat.alld <- pred_is.dat.shrink %>% 
  select(c('entity_ID','AMdebt','EMdebt','ORCdebt','NMdebt')) %>%
  gather(key = "myctype", value = "debt", AMdebt, EMdebt, ORCdebt, NMdebt) %>%
  mutate(myctype = case_when(myctype=="AMdebt" ~ "AM",                                   
                             myctype=="EMdebt" ~ "EM", 
                             myctype=="ORCdebt" ~ "ORC", 
                             myctype=="NMdebt" ~ "NM")) 

pred_is.dat.allcd <- pred_is.dat.shrink %>% 
  select(c('entity_ID','C_AMdebt','C_EMdebt','C_ORCdebt','C_NMdebt')) %>%
  gather(key = "myctype", value = "debt.c", C_AMdebt, C_EMdebt, C_ORCdebt, C_NMdebt) %>%
  mutate(myctype = case_when(myctype == "C_AMdebt" ~ "AM",                                   
                             myctype == "C_EMdebt" ~ "EM", 
                             myctype == "C_ORCdebt" ~ "ORC", 
                             myctype == "C_NMdebt" ~ "NM"))

pred_is.dat.all.diff <- pred_is.dat.shrink %>% 
  select(c('entity_ID','entity_class2', 'sprich','latitude','longitude','abslatitude','AMdiff','EMdiff','ORCdiff','NMdiff','dist','area','elev_range', 'temp','prec')) %>%
  gather(key = "myctype", value = "diff", AMdiff, EMdiff, ORCdiff, NMdiff) %>%
  mutate(myctype = case_when(myctype == "AMdiff" ~ "AM",                                   
                             myctype == "EMdiff" ~ "EM", 
                             myctype == "ORCdiff" ~ "ORC", 
                             myctype == "NMdiff" ~ "NM")) 

pred_is.dat.all.exp <- pred_is.dat.shrink %>% 
  select(c('entity_ID','AM_exp','EM_exp','ORC_exp','NM_exp')) %>%
  gather(key = "myctype", value = "exp", AM_exp, EM_exp, ORC_exp,NM_exp) %>%
  mutate(myctype = case_when(myctype == "AM_exp" ~ "AM",                                   
                             myctype == "EM_exp" ~ "EM", 
                             myctype == "ORC_exp" ~ "ORC", 
                             myctype == "NM_exp" ~ "NM")) 


pred_is.dat.all.obs <- pred_is.dat.shrink %>% 
  select(c('entity_ID', 'AM', 'EM', 'ORC', 'NM')) %>%
  gather(key = "myctype", value = "obs", AM, EM, ORC, NM) 

pred_is.dat.all.T <- pred_is.dat.shrink %>% select(entity_ID, sprichdiff)

pred_is.dat.all <- pred_is.dat.all.diff %>%
  left_join(pred_is.dat.all.exp, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.all.obs, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.alld, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.allcd, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.all.T, by = c("entity_ID")) %>%
  mutate(debt.weights = exp, debt.c.weights = abs(sprichdiff)) %>%
  mutate(myctype = as.factor(myctype))

pred_is.dat.all <- within(pred_is.dat.all, myctype <- relevel(myctype, ref = "NM"))

############################
######## WRITE DATA ########
############################

#debt
saveRDS(pred_is.dat.shrink,"data/deficit_myc_latitude_data_2025.RDS")

#cdebt
saveRDS(pred_is.dat.all,"data/cdeficit_myc_latitude_data_2025.RDS")
