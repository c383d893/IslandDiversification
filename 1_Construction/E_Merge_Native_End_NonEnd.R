# This script 
# 1. joins native, endemic and nonend counts for each mutualism (4)

##########################
###### LOAD PACKAGES #####
##########################

library(tidyverse)

###################
#### EXTENDED #####
###################

###################
#### LOAD DATA ####
###################

# poll
poll.nat <- readRDS("data/native_poll_data_2025.RDS") %>% mutate(flora = "native") 
poll.nat.geo <- poll.nat %>% select(-c('sprich', 'poll.b', 'poll.ab', 'flora'))

poll.nat.end <- readRDS("data/end_poll_data_2025.RDS") %>% mutate(flora = "endemic") %>% 
  select(c('entity_ID', 'sprich', 'poll.b', 'poll.ab', 'flora'))
poll.nat.end <- poll.nat.geo %>% left_join(poll.nat.end, by = "entity_ID") %>% mutate(flora = "endemic")

poll.nat.nonend <- readRDS("data/nonend_poll_data_2025.RDS") %>% mutate(flora = "nonendemic") %>% 
  select(c('entity_ID', 'sprich', 'poll.b', 'poll.ab', 'flora'))
poll.nat.nonend <- poll.nat.geo %>% left_join(poll.nat.nonend, by = "entity_ID") %>% mutate(flora = "nonendemic")

poll <- rbind(poll.nat, poll.nat.end, poll.nat.nonend) %>%
  mutate(sprich = ifelse(is.na(sprich), 0, sprich)) %>%
  mutate(poll.b = ifelse(is.na(poll.b), 0, poll.b)) %>%
  mutate(poll.ab = ifelse(is.na(poll.ab), 0, poll.ab)) 

saveRDS(poll, "data/merge_poll_data_2025.RDS")

# disp
disp.nat <- readRDS("data/native_disp_data_2025.RDS") %>% mutate(flora = "native") 
disp.nat.geo <- disp.nat %>% select(-c('sprich', 'disp.b', 'disp.ab', 'flora'))

disp.nat.end <- readRDS("data/end_disp_data_2025.RDS") %>% mutate(flora = "endemic") %>% 
  select(c('entity_ID', 'sprich', 'disp.b', 'disp.ab', 'flora'))
disp.nat.end <- disp.nat.geo %>% left_join(disp.nat.end, by = "entity_ID") %>% mutate(flora = "endemic")

disp.nat.nonend <- readRDS("data/nonend_disp_data_2025.RDS") %>% mutate(flora = "nonendemic") %>% 
  select(c('entity_ID', 'sprich', 'disp.b', 'disp.ab', 'flora'))
disp.nat.nonend <- disp.nat.geo %>% left_join(disp.nat.nonend, by = "entity_ID") %>% mutate(flora = "nonendemic")

disp <- rbind(disp.nat, disp.nat.end, disp.nat.nonend) %>%
  mutate(sprich = ifelse(is.na(sprich), 0, sprich)) %>%
  mutate(disp.b = ifelse(is.na(disp.b), 0, disp.b)) %>%
  mutate(disp.ab = ifelse(is.na(disp.ab), 0, disp.ab)) 

saveRDS(disp, "data/merge_disp_data_2025.RDS")

# myc
myc.nat <- readRDS("data/native_myc_data_2025.RDS") %>% mutate(flora = "native") 
myc.nat.geo <- myc.nat %>% select(-c('sprich', 'AM', 'EM', 'ORC','NM', 'flora'))

myc.nat.end <- readRDS("data/end_myc_data_2025.RDS") %>% mutate(flora = "endemic") %>% 
  select(c('entity_ID', 'sprich', 'AM', 'EM', 'ORC','NM', 'flora'))
myc.nat.end <- myc.nat.geo %>% left_join(myc.nat.end, by = "entity_ID") %>% mutate(flora = "endemic")

myc.nat.nonend <- readRDS("data/nonend_myc_data_2025.RDS") %>% mutate(flora = "nonendemic") %>% 
  select(c('entity_ID', 'sprich', 'AM', 'EM', 'ORC','NM', 'flora'))
myc.nat.nonend <- myc.nat.geo %>% left_join(myc.nat.nonend, by = "entity_ID") %>% mutate(flora = "nonendemic")

myc <- rbind(myc.nat, myc.nat.end, myc.nat.nonend) %>%
  mutate(sprich = ifelse(is.na(sprich), 0, sprich)) %>%
  mutate(AM = ifelse(is.na(AM), 0, AM)) %>%
  mutate(EM = ifelse(is.na(EM), 0, EM)) %>%
  mutate(ORC = ifelse(is.na(ORC), 0, ORC)) %>%
  mutate(NM = ifelse(is.na(NM), 0, NM))

saveRDS(myc, "data/merge_myc_data_2025.RDS")

# nfix
nfix.nat <- readRDS("data/native_nfix_data_2025.RDS") %>% mutate(flora = "native") 
nfix.nat.geo <- nfix.nat %>% select(-c('sprich', 'nfix', 'nonfix', 'flora'))

nfix.nat.end <- readRDS("data/end_nfix_data_2025.RDS") %>% mutate(flora = "endemic") %>% 
  select(c('entity_ID', 'sprich', 'nfix', 'nonfix', 'flora'))
nfix.nat.end <- nfix.nat.geo %>% left_join(nfix.nat.end, by = "entity_ID") %>% mutate(flora = "endemic")

nfix.nat.nonend <- readRDS("data/nonend_nfix_data_2025.RDS") %>% mutate(flora = "nonendemic") %>% 
  select(c('entity_ID', 'sprich', 'nfix', 'nonfix', 'flora'))
nfix.nat.nonend <- nfix.nat.geo %>% left_join(nfix.nat.nonend, by = "entity_ID") %>% mutate(flora = "nonendemic")

nfix <- rbind(nfix.nat, nfix.nat.end, nfix.nat.nonend) %>%
  mutate(sprich = ifelse(is.na(sprich), 0, sprich)) %>%
  mutate(nfix = ifelse(is.na(nfix), 0, nfix)) %>%
  mutate(nonfix = ifelse(is.na(nonfix), 0, nonfix)) 

saveRDS(nfix, "data/merge_nfix_data_2025.RDS")
