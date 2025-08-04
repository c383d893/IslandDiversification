# This script:
# 1. generates proportions of endemics assigned 
# 2. generates proportions of endemics assigned per mutualism (4)

##########################
###### LOAD PACKAGES #####
##########################

library(tidyverse)

###########################
##### CHECK ENDNONEND #####
###########################

dat <- readRDS("data/endnonend_data_2025.RDS")

# check proportion assignment of end + non
dat.check <- dat %>%
  mutate(endnonend = endemic + nonendemic) %>%
  mutate(propnatassigned = endnonend/sprich) %>%
  select(entity_ID, propnatassigned)

saveRDS(dat.check, "data/endnonend_assignment_2025.RDS")

###########################
##### CHECK ENDNONEND #####
########### POLL ##########
###########################

p.nat <- readRDS("data/native_poll_data_2025.RDS") %>%
  mutate(poll.as = poll.b + poll.ab) %>%
  select(c('entity_ID', 'sprich','poll.as'))

p.end <- readRDS("data/end_poll_data_2025.RDS") %>%
  mutate(endemic = poll.b + poll.ab) %>%
  select(c('entity_ID', 'endemic'))

p.nend <- readRDS("data/nonend_poll_data_2025.RDS") %>%
  mutate(nonendemic = poll.b + poll.ab) %>%
  select(c('entity_ID', 'nonendemic'))

p.check <- p.nat %>% left_join(p.end) %>% left_join(p.nend) %>%
  mutate_at(vars(c('endemic','nonendemic')), ~replace(., is.na(.), 0)) %>%
  mutate(endnonend = endemic + nonendemic) %>%
  mutate(propnatassigned = endnonend/sprich) %>%
  select(entity_ID, propnatassigned)

saveRDS(p.check, "data/poll_endnon_assignment_2025.RDS")

###########################
##### CHECK ENDNONEND #####
########### DISP ##########
###########################

d.nat <- readRDS("data/native_disp_data_2025.RDS") %>%
  mutate(disp.as = disp.b + disp.ab) %>%
  select(c('entity_ID', 'sprich','disp.as'))

d.end <- readRDS("data/end_disp_data_2025.RDS") %>%
  mutate(endemic = disp.b + disp.ab) %>%
  select(c('entity_ID', 'endemic'))

d.nend <- readRDS("data/nonend_disp_data_2025.RDS") %>%
  mutate(nonendemic = disp.b + disp.ab) %>%
  select(c('entity_ID', 'nonendemic'))

d.check <- d.nat %>% left_join(d.end) %>% left_join(d.nend) %>%
  mutate_at(vars(c('endemic','nonendemic')), ~replace(., is.na(.), 0)) %>%
  mutate(endnonend = endemic + nonendemic) %>%
  mutate(propnatassigned = endnonend/sprich) %>%
  select(entity_ID, propnatassigned)

saveRDS(d.check, "data/disp_endnon_assignment_2025.RDS")

###########################
##### CHECK ENDNONEND #####
########### MYC ###########
###########################

m.nat <- readRDS("data/native_myc_data_2025.RDS")  %>%
  mutate(myc.as = AM + EM + ORC + NM) %>%
  select(c('entity_ID', 'sprich','myc.as'))

m.end <- readRDS("data/end_myc_data_2025.RDS") %>%
  mutate(endemic = AM + EM + ORC + NM) %>%
  select(c('entity_ID', 'endemic'))

m.nend <- readRDS("data/nonend_myc_data_2025.RDS") %>%
  mutate(nonendemic = AM + EM + ORC + NM) %>%
  select(c('entity_ID', 'nonendemic'))

m.check <- m.nat %>% left_join(m.end) %>% left_join(m.nend) %>%
  mutate_at(vars(c('endemic','nonendemic')), ~replace(., is.na(.), 0)) %>%
  mutate(endnonend = endemic + nonendemic) %>%
  mutate(propnatassigned = endnonend/sprich) %>%
  select(entity_ID, propnatassigned) 

saveRDS(m.check, "data/myc_endnon_assignment_2025.RDS")

###########################
##### CHECK ENDNONEND #####
########### NFIX ##########
###########################

n.nat <- readRDS("data/native_nfix_data_2025.RDS") %>%
  mutate(nfix.as = nfix + nonfix) %>%
  select(c('entity_ID', 'sprich','nfix.as'))

n.end <- readRDS("data/end_nfix_data_2025.RDS") %>%
  mutate(endemic = nfix + nonfix) %>%
  select(c('entity_ID', 'endemic'))

n.nend <- readRDS("data/nonend_nfix_data_2025.RDS") %>%
  mutate(nonendemic = nfix + nonfix) %>%
  select(c('entity_ID', 'nonendemic'))

n.check <- n.nat %>% left_join(n.end) %>% left_join(n.nend) %>%
  mutate_at(vars(c('endemic','nonendemic')), ~replace(., is.na(.), 0)) %>%
  mutate(endnonend = endemic + nonendemic) %>%
  mutate(propnatassigned = endnonend/sprich) %>%
  select(entity_ID, propnatassigned)

saveRDS(n.check, "data/nfix_endnon_assignment_2025.RDS")
