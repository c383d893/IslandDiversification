# This script:
# 1. calculates island species deficit for pollination syndrome in nonendemics

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

dat <-  readRDS("data/merge_poll_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area",  "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range", "poll.b", "poll.ab")) %>%
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

dat2 <-  readRDS("data/merge_poll_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>%
  filter(flora == "nonendemic") %>%
  filter(sprich > 0) %>%
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "elev_range", "dist", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec", 
           "poll.b", "poll.ab")) %>%
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
####### MAINLAND BIOTIC POLL ########
#####################################

gam.mod.poll.b <- gam(poll.b ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.poll.b)

mod <- gam.mod.poll.b
new.dat.poll.b <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.poll.b <- predict.gam(mod,newdata = new.dat.poll.b, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.poll.b$abslatitude, poll = "poll.b") 

# check assumptions
gam.check(gam.mod.poll.b)
disp_check(gam.mod.poll.b,dat)

#####################################
####### MAINLAND ABIOTIC POLL #######
#####################################

gam.mod.poll.ab <- gam(poll.ab ~ s(abslatitude) , family=nb(link="log"), data = dat.ml) 
summary(gam.mod.poll.ab)

mod <- gam.mod.poll.ab
new.dat.poll.ab <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.poll.ab <- predict.gam(mod, newdata = new.dat.poll.ab, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.poll.ab$abslatitude, poll = "poll.ab") 

# check assumptions
gam.check(gam.mod.poll.ab)
disp_check(gam.mod.poll.ab,dat)

############################
########## PLOT ############
##### LAT BY POLL TYPE #####
############################

pred.mainland <- rbind(pred.poll.b,pred.poll.ab)

dat.ml.cond <- dat.ml %>%
  select(abslatitude, poll.b, poll.ab) %>%
  gather(key = "poll", value = "sprich", poll.b, poll.ab)

# create a custom color scale
colScale <- scale_colour_manual(values=c("darkgrey", "deeppink3"))
fillScale <- scale_fill_manual(values=c("darkgrey","deeppink3"))

pred.mainland$poll <- ordered(pred.mainland$poll, levels = c("poll.b", "poll.ab"))
dat.ml.cond$poll <- ordered(dat.ml.cond$poll, levels = c("poll.b", "poll.ab"))

# rename poll for legend:
pred.mainland <- pred.mainland %>% mutate(poll = case_when(poll=="poll.b" ~ "Biotic",
                                                           poll=="poll.ab" ~ "Abiotic"))

dat.ml.cond <- dat.ml.cond %>% mutate(poll = case_when(poll=="poll.b" ~ "Biotic",
                                                       poll=="poll.ab" ~ "Abiotic"))
lat.polltype <-
  ggplot(pred.mainland, aes(x = abslatitude, y = fit,color = poll, fill = poll))+
  geom_line(size = 1) +
  geom_point(data = dat.ml.cond, aes(x = abslatitude, y = sprich, color = factor(poll)), alpha = 0.3, size = 4)+ 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.5) + 
  xlab("Absolute latitude") +
  ylab("Species richness")+
  theme_classic(base_size = 25)+
  ylim(0,2500)+
  colScale+
  fillScale+
  theme(legend.position = 'none')+
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  guides(color = guide_legend(title="Pollination \nSyndrome",override.aes=list(fill =NA,size =3, alpha =0.7,linetype = c(0, 0))))+
  theme(axis.text.x = element_text(size =20),axis.text.y = element_text(angle = 45,size=20))

# write out
png("figures/poll_byabslat.jpg", width=10, height= 10, units='in', res=300)
lat.polltype
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

pred_poll.b <- predict.gam(gam.mod.poll.b,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_poll.b_df <- cbind(dat.is.min,pred_poll.b) %>% rename(poll.b_exp = fit)

pred_poll.ab <- predict.gam(gam.mod.poll.ab,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_poll.ab_df <- cbind(dat.is.min,pred_poll.ab) %>% rename(poll.ab_exp = fit)

pred_is.dat <- pred_sprich_df %>%
  left_join(dat2, by= c('entity_ID','abslatitude')) %>%
  mutate(sprichdiff = sprich_exp - sprich) %>%
  mutate(Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt))

pred_is.dat <- pred_sprich_df %>%
  left_join(dat2, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_poll.b_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_poll.ab_df, by= c('entity_ID','abslatitude')) %>%
  mutate_at(c('poll.b_exp','poll.ab_exp','sprich_exp'), as.integer) %>%
  mutate(propb_exp = poll.b_exp/(poll.b_exp + poll.ab_exp))

saveRDS(pred_is.dat,"data/exp_poll_latitude_data_2025.RDS")

############################
######## CREATE DATA #######
############################

pred_is.dat.shrink <- pred_is.dat %>%
  mutate(sprich = poll.b + poll.ab, sprich_exp = poll.b_exp + poll.ab_exp) %>%
  mutate(poll.b.diff = poll.b_exp - poll.b, poll.ab.diff = poll.ab_exp - poll.ab, sprichdiff = sprich_exp - sprich) %>%
  mutate(poll.b.debt = (poll.b.diff/poll.b_exp), poll.ab.debt = (poll.ab.diff/poll.ab_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(poll.b.debt = ifelse(poll.b.debt < 0, 0, poll.b.debt)) %>% mutate(poll.ab.debt = ifelse(poll.ab.debt <0, 0, poll.ab.debt)) %>% mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(poll.b.debt = ifelse(poll.b.debt > 1, 1, poll.b.debt)) %>% mutate(poll.ab.debt = ifelse(poll.ab.debt >1, 1, poll.ab.debt)) %>% mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt)) %>%
  mutate(C_poll.b.debt = (poll.b.diff/sprichdiff), C_poll.ab.debt = (poll.ab.diff/sprichdiff)) %>% 
  mutate(C_poll.b.debt = ifelse(C_poll.b.debt < 0, 0, C_poll.b.debt)) %>% mutate(C_poll.ab.debt = ifelse(C_poll.ab.debt <0, 0, C_poll.ab.debt)) %>% 
  mutate(C_poll.b.debt = ifelse(C_poll.b.debt > 1, 1, C_poll.b.debt)) %>% mutate(C_poll.ab.debt = ifelse(C_poll.ab.debt >1, 1, C_poll.ab.debt)) 

############################
######## CDEBT DATA ########
############################

pred_is.dat.alld <- pred_is.dat.shrink %>% 
  select(c('entity_ID','poll.b.debt','poll.ab.debt')) %>%
  gather(key = "poll", value = "debt", poll.b.debt, poll.ab.debt) %>%
  mutate(poll = case_when(poll == "poll.b.debt" ~ "poll.b",                                   
                          poll == "poll.ab.debt" ~ "poll.ab")) 

pred_is.dat.allcd <- pred_is.dat.shrink %>% 
  select(c('entity_ID','C_poll.b.debt','C_poll.ab.debt')) %>%
  gather(key = "poll", value = "debt.c", C_poll.b.debt, C_poll.ab.debt)%>%
  mutate(poll = case_when(poll == "C_poll.b.debt" ~ "poll.b",                                   
                          poll == "C_poll.ab.debt" ~ "poll.ab"))

pred_is.dat.all.diff <- pred_is.dat.shrink %>% 
  select(c('entity_ID','entity_class2', 'sprich','latitude','longitude','abslatitude','poll.b.diff','poll.ab.diff','dist','area','elev_range', 'temp','prec')) %>%
  gather(key = "poll", value = "diff", poll.b.diff, poll.ab.diff) %>%
  mutate(poll = case_when(poll == "poll.b.diff" ~ "poll.b",                                   
                          poll == "poll.ab.diff" ~ "poll.ab")) 

pred_is.dat.all.exp <- pred_is.dat.shrink %>% 
  select(c('entity_ID','poll.b_exp','poll.ab_exp')) %>%
  gather(key = "poll", value = "exp", poll.b_exp, poll.ab_exp) %>%
  mutate(poll = case_when(poll == "poll.b_exp" ~ "poll.b",                                   
                          poll == "poll.ab_exp" ~ "poll.ab")) 

pred_is.dat.all.obs <- pred_is.dat.shrink %>% 
  select(c('entity_ID', 'poll.b', 'poll.ab')) %>%
  gather(key = "poll", value = "obs", poll.b, poll.ab) 

pred_is.dat.all.T <- pred_is.dat.shrink %>% select(entity_ID, sprichdiff)

pred_is.dat.all <- pred_is.dat.all.diff %>%
  left_join(pred_is.dat.all.exp, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.all.obs, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.alld, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.allcd, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.all.T, by = c("entity_ID")) %>%
  mutate(debt.weights = exp, debt.c.weights = abs(sprichdiff)) %>%
  mutate(poll = as.factor(poll)) 

pred_is.dat.all <- within(pred_is.dat.all, poll <- relevel(poll, ref = "poll.ab"))

############################
######## WRITE DATA ########
############################

#debt
saveRDS(pred_is.dat.shrink,"data/deficit_poll_latitude_data_2025.RDS")

#cdebt
saveRDS(pred_is.dat.all,"data/cdeficit_poll_latitude_data_2025.RDS")
