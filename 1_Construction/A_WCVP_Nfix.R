# This script:
# 1. generates Nfix species harmonized to WCVP taxonomy

##########################
###### LOAD PACKAGES #####
##########################

library(tidyverse)
#library(devtools)
#devtools::install_github("matildabrown/rWCVP")
#remotes::install_github("matildabrown/rWCVPdata")
library(rWCVP)
library(rWCVPdata)

###################
#### LOAD DATA ####
###################

ndat <- read.csv("data/Werner_NFix.csv", header = TRUE) %>%               
  rename(ndatspecies = species)%>%
  select(ndatspecies, data_fixing)%>% distinct(ndatspecies, .keep_all=TRUE)

# match Nfix to GIFT (WCVP nomenclature)
set.seed(156)
matches <- wcvp_match_names(ndat,
                            name_col="ndatspecies",
                            fuzzy=TRUE,
                            progress_bar=FALSE)

# subset to accepted and synonyms and choose best match
matches_w_synonyms <- matches %>%
  filter(wcvp_status %in% c("Accepted", "Synonym", "Orthographic")) %>%
  arrange(ndatspecies, match_type, wcvp_status) %>%  
  group_by(ndatspecies) %>%     
  slice(1) %>%  
  ungroup() 

# add family information
wcvp_fam <- wcvp_names %>% 
  mutate(species = paste(genus,species,sep=" ")) %>%
  select(c("species","family")) %>%
  distinct(species, .keep_all = TRUE) 

# join together
matches_all <- left_join(matches_w_synonyms, wcvp_fam, by = c("ndatspecies"="species")) %>%
  rename(family_wcvp=family) 

# clean 
wcvp_final <- matches_all %>% select(c("ndatspecies","data_fixing","family_wcvp")) %>%
  rename(species = ndatspecies, family = family_wcvp) %>%
  drop_na()

########################
####### SAVE DATA ######
########################

write_csv(wcvp_final, "data/Werner_NFix_WVCP.csv")
