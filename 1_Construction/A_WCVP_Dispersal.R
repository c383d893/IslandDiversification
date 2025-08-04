# This script:
# 1. generates TRY derived species harmonized to WCVP taxonomy
# 2. generates family proportions for dispersal syndrome

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

# Initial TRY data: 2446500 lines
# TRY data was cleaned: spaces to underscores, empty cells to NA; removed 0,1,2,3 for trait value
# Syndromes synonymized using: https://nutnet.org/sites/default/files/Leda_dispersability%20traits.pdf
# should be 457418 traits

try.dat <- read.table("data/Try_Dispersal.txt", header = TRUE) %>%
  mutate(dispersal = case_when(Dispersal == "abiotic" ~ "abiotic",
                               Dispersal == "animal" ~ "biotic",   
                               Dispersal == "human" ~ "abiotic",
                               Dispersal == "self" ~ "abiotic",
                               Dispersal == "wind" ~ "abiotic",
                               Dispersal == "water" ~ "abiotic")) %>%
  rename(tryspecies = SpeciesName) %>%
  select(c("tryspecies","dispersal")) %>%
  # split species and keep only first two columns; remove problematic rows
  separate(tryspecies, c("genus", "species"), sep = "_", remove = TRUE) %>%
  filter(species != "") %>%
  filter(!str_detect(species, "\\d")) %>%
  filter(species != "(~bushii)") %>%
  filter(species != "&") %>%
  filter(species != "(Panicum)") %>%
  filter(species != "(L.)") %>%
  filter(species != "A.") %>%
  filter(species != "A.E.") %>%
  mutate(tryspecies = paste(genus,species,sep=" ")) %>%
  drop_na() %>%
  select(tryspecies,dispersal) %>%
  #find consensus dispersal based on majority assigned
  group_by(tryspecies) %>% 
  summarise(dispersal = names(which.max(table(dispersal))), .groups = "drop")

# match TRY (TPL nomenclature) to GIFT (WCVP nomenclature): https://matildabrown.github.io/rWCVP/articles/redlist-name-matching.html
set.seed(156)
matches <- wcvp_match_names(try.dat,
                            name_col="tryspecies",
                            fuzzy=TRUE,
                            progress_bar=FALSE)

# subset to accepted and synonyms and choose best match
matches_w_synonyms <- matches %>%
  filter(wcvp_status %in% c("Accepted", "Synonym", "Orthographic")) %>%
  arrange(tryspecies, match_type, wcvp_status) %>%  
  group_by(tryspecies) %>%     
  slice(1) %>%  
  ungroup() 

# add family information
wcvp_fam <- wcvp_names %>% 
  mutate(species = paste(genus,species,sep=" ")) %>%
  select(c("species","family")) %>%
  distinct(species, .keep_all = TRUE) 

# join together
matches_all <- left_join(matches_w_synonyms, wcvp_fam, by = c("tryspecies"="species")) %>%
  rename(family_wcvp=family) 

# clean
wcvp_final <- matches_all %>% select(c("tryspecies","dispersal","family_wcvp")) %>%
  rename(species = tryspecies, family = family_wcvp) %>%
  drop_na()

wcvp_final_sp <- wcvp_final %>% select(c("species","dispersal"))

# generate consensus family level dispersal syndrome
matches_fam <- wcvp_final %>%
  group_by(family) %>% 
  summarise(dispersal = names(which.max(table(dispersal))), .groups = "drop")

########################
####### SAVE DATA ######
########################

saveRDS(wcvp_final_sp, "data/dispersal_species_WVCP.RDS")
saveRDS(matches_fam, "data/dispersal_family_WVCP.RDS")
