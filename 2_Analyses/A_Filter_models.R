# This script models and plots results for:
# 1. native and non-endemic filters for each mutualism (4)

############################
###### LOAD PACKAGES #######
############################

library(mgcv); library(gridExtra); library(betareg); library(MASS); library(lme4); library(lmerTest); library(lsmeans); library(ggeffects); library(spdep); library(ggplot2); library(ncf); library(ape); library(sjPlot); library(gridExtra); library(MuMIn);library(maps); library(sf); library(car);library(viridis);library(tidyverse);library(glmmTMB);library(DHARMa);library(purrr)
options(na.action = "na.fail")

source("Master_modelling.R")

# set manual color
colScalef <- scale_colour_manual(name = "polltype", values = c("grey20","darkslategray4", "darkslategray3"))
colFillf <- scale_fill_manual(name = "polltype", values = c("grey20","darkslategray4", "darkslategray3"))

####################################
####################################
####################################
############ POLLFILTER ############
####################################
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "polltype", values = c("darkgrey","deeppink3"))
colFill <- scale_fill_manual(name = "polltype", values = c("darkgrey","deeppink3"))

# endemism data quality
p.check <- readRDS("data/poll_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# input data not scaled
dat.start <-  readRDS("data/merge_poll_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>% 
  dplyr::select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range","poll.b", "poll.ab","flora")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class == "Mainland" ~ "Mainland")) %>%
  dplyr::select(-geology) %>%          
  filter(area > 6) %>% 
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>% 
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range), 1, elev_range)) %>% 
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), 
         temp = as.vector(scale(log10((temp)+.01))),
         prec = as.vector(scale(log10((prec)+.01))),
         elev_range = as.vector(scale(log10((elev_range)+.01)))) %>% 
  drop_na() %>%
  left_join(p.check, by = "entity_ID") %>%
  mutate(pnata = ifelse(is.na(pnata), 0, pnata)) 

# non-endemic dat: 479
dat.nonend <- dat.start %>% filter(flora == "nonendemic") %>%
  filter(!(poll.b == 0 & poll.ab == 0)) %>%
  mutate(prop.b = poll.b/ (poll.b+ poll.ab))

############################
### NONEND FILTER MODELS ###
############################

dat <- dat.nonend
mod <- glmer(cbind(poll.b, poll.ab) ~ entity_class2 + abslatitude +  prec + temp + area + elev_range + (1|entity_class2:entity_ID), weights = pnata,
           data = dat, family = binomial)

# calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

# generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 5000, resamp = 100, latlon=TRUE)
plot(sp)

# rac model
rac <- Spat.cor.rep(mod,dat,2000)
mod.rac <- glmer(cbind(poll.b, poll.ab) ~ entity_class2 + abslatitude +  prec + temp + area + elev_range + rac + (1|entity_class2:entity_ID), weights = pnata,
               data = dat, family = binomial)
summary(mod.rac)

# full plot proportion 
ref <- emmeans(mod.rac, pairwise ~ entity_class2, type = "response")
ref.table <- as.data.frame(ref$emmeans)

plot <- 
  ggplot(ref.table, aes(entity_class2, prob, color = entity_class2, fill = entity_class2)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(entity_class2, prop.b), position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Proportion biotically pollinated") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  colFillf + colScalef +
  scale_x_discrete(labels = c("Mainland", "Non-oceanic", "Oceanic")) +
  ylim(0,1) 

png("figures/islandfilter_poll_nonend.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

####################################
####################################
########### DISP FILTER ############
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "disptype", values = c("darkgrey","brown4"))
colFill <- scale_fill_manual(name = "disptype", values = c("darkgrey","brown4"))

# endemism data quality
d.check <- readRDS("data/disp_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# input data not scaled
dat.start <-  readRDS("data/merge_disp_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>% 
  dplyr::select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range","disp.b", "disp.ab","flora")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class == "Mainland" ~ "Mainland")) %>%
  dplyr::select(-geology) %>%          
  filter(area > 6) %>% 
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>% 
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>% 
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), 
         temp = as.vector(scale(temp)), 
         prec = as.vector(scale(log10((prec)+.01))),
         elev_range = as.vector(scale(log10((elev_range)+.01))),
         abslatitude = as.vector(scale(abslatitude))) %>% 
  drop_na() %>%
  left_join(d.check, by = "entity_ID") %>%
  mutate(pnata = ifelse(is.na(pnata), 0, pnata)) 

# non-endemic dat: 480
dat.nonend <- dat.start %>% filter(flora == "nonendemic") %>%
  filter(!(disp.b == 0 & disp.ab == 0)) %>%
  mutate(prop.b = disp.b/ (disp.b+ disp.ab))

############################
### NONEND FILTER MODELS ###
############################

dat <- dat.nonend
mod <- glmer(cbind(disp.b, disp.ab) ~ entity_class2 + abslatitude +  prec + temp + area + elev_range + (1|entity_class2:entity_ID), weights = pnata,
           data = dat, family = binomial)

# calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

# generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 5000, resamp = 100, latlon=TRUE)
plot(sp)

# rac model
rac <- Spat.cor.rep(mod,dat,1500)
mod.rac <- glmer(cbind(disp.b, disp.ab) ~ entity_class2 + abslatitude +  prec + temp + area + elev_range + rac + (1|entity_class2:entity_ID), weights = pnata,
               data = dat, family = binomial)
summary(mod.rac)

# full plot proportion 
ref <- emmeans(mod.rac, pairwise ~ entity_class2, type = "response")
ref.table <- as.data.frame(ref$emmeans)

plot <- 
  ggplot(ref.table, aes(entity_class2, prob, color = entity_class2, fill = entity_class2)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(entity_class2, prop.b), position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Proportion biotically dispersed") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  colFillf + colScalef +
  scale_x_discrete(labels = c("Mainland", "Non-oceanic", "Oceanic")) +
  ylim(0,1) 

png("figures/islandfilter_disp_nonend.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

####################################
####################################
############ MYC FILTER ############
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "myctype", values = c("darkgrey", "royalblue4"))
colFill <- scale_fill_manual(name = "myctype", values = c("darkgrey", "royalblue4"))

# endemism data quality
m.check <- readRDS("data/myc_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# input data not scaled
dat.start <-  readRDS("data/merge_myc_data_2025.RDS") %>%
  filter(!entity_class == "undetermined") %>% 
  dplyr::select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range","AM","EM","ORC","NM","flora")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class == "Mainland" ~ "Mainland")) %>%
  dplyr::select(-geology) %>%          
  filter(area > 6) %>%
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>% 
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>% 
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), 
         temp = as.vector(scale(temp)), 
         prec = as.vector(scale(log10((prec)+.01))),
         elev_range = as.vector(scale(log10((elev_range)+.01))),
         abslatitude = as.vector(scale(abslatitude))) %>% 
  drop_na() %>%
  left_join(m.check, by = "entity_ID") %>%
  mutate(pnata = ifelse(is.na(pnata), 0, pnata)) 

# non-endemic dat: 480
dat.nonend <- dat.start %>% filter(flora == "nonendemic") %>% 
  filter(!(AM == 0 & NM == 0))  %>%
  mutate(prop.AM = AM/ (AM + EM + ORC + NM))

############################
### NONEND FILTER MODELS ###
############################

dat <- dat.nonend 
mod <- glmer(cbind(AM, NM) ~ entity_class2 + abslatitude + prec + temp + area + elev_range + (1|entity_class2:entity_ID), weights = pnata,
           data = dat, family = binomial)

# calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

# generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 5000, resamp = 100, latlon=TRUE)
plot(sp)

# rac model
rac <- Spat.cor.rep(mod,dat,1000)
mod.rac <- glmer(cbind(AM, NM) ~ entity_class2 + abslatitude +  prec + temp + area + elev_range + rac + (1|entity_class2:entity_ID), weights = pnata,
               data = dat, family = binomial)
summary(mod.rac)

# full plot proportion 
ref <- emmeans(mod.rac, pairwise ~ entity_class2, type = "response")
ref.table <- as.data.frame(ref$emmeans)

plot <- 
  ggplot(ref.table, aes(entity_class2, prob, color = entity_class2, fill = entity_class2)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(entity_class2, prop.AM), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Proportion AM plant species") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  colFillf + colScalef +
  scale_x_discrete(labels = c("Mainland", "Non-oceanic", "Oceanic")) +
  ylim(0,1) 

png("figures/islandfilter_myc_nonend.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

####################################
####################################
############ NFIX FILTER ###########
####################################
####################################

# set manual color
colScale <- scale_colour_manual(name = "nfixtype", values = c("darkseagreen3", "darkgrey"))
colFill <- scale_fill_manual(name = "nfixtype", values = c("darkseagreen3", "darkgrey"))

# endemism data quality
n.check <- readRDS("data/nfix_endnon_assignment_2025.RDS") %>% filter(!propnatassigned == 0) %>% rename(pnata = propnatassigned)

# input data not scaled
dat.start <-  readRDS("data/merge_nfix_data_2025.RDS") %>%
  filter(entity_ID %in% n.check$entity_ID) %>%
  filter(!entity_class == "undetermined") %>% 
  dplyr::select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range","nfix", "nonfix","flora")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class == "Mainland" ~ "Mainland")) %>%
  dplyr::select(-geology) %>%          
  filter(area > 6) %>%
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>% 
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>% 
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), 
         temp = as.vector(scale(temp)), 
         prec = as.vector(scale(log10((prec)+.01))),
         elev_range = as.vector(scale(log10((elev_range)+.01))),
         abslatitude = as.vector(scale(abslatitude))) %>% 
  drop_na() %>%
  left_join(n.check, by = "entity_ID") %>%
  mutate(pnata = ifelse(is.na(pnata), 0, pnata)) 

# non-endemic dat: 477
dat.nonend <- dat.start %>% filter(flora == "nonendemic") %>% 
  filter(!(nfix == 0 & nonfix == 0)) %>%
  mutate(prop.nfix = nfix/ (nfix + nonfix))

############################
### NONEND FILTER MODELS ###
############################

dat <- dat.nonend
mod <- glmer(cbind(nfix, nonfix) ~ entity_class2 + abslatitude +  prec + temp + area + elev_range + (1|entity_class2:entity_ID), weights = pnata,
           data = dat, family = binomial)

# calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

# generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 5000, resamp = 100, latlon=TRUE)
plot(sp)

# rac model
rac <- Spat.cor.rep(mod,dat,2000)
mod.rac <- glmer(cbind(nfix, nonfix) ~ entity_class2 + abslatitude +  prec + temp + area + elev_range + rac + (1|entity_class2:entity_ID), weights = pnata,
               data = dat, family = binomial)
summary(mod.rac)

# full plot proportion 
ref <- emmeans(mod.rac, pairwise ~ entity_class2, type = "response")
ref.table <- as.data.frame(ref$emmeans)

plot <- 
  ggplot(ref.table, aes(entity_class2, prob, color = entity_class2, fill = entity_class2)) +  
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size = 4) + 
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(entity_class2, prop.nfix), position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  theme_minimal(base_size = 40) +
  ylab("Proportion Nfixing plant species") +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="none") +
  colFillf + colScalef +
  scale_x_discrete(labels = c("Mainland", "Non-oceanic", "Oceanic")) +
  ylim(0,1) 

png("figures/islandfilter_nfix_nonend.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()
