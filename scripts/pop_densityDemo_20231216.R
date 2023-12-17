# Jeremy Summers
# April 2023

# This script calculates the number of adult females and males
# within a well-censused area of the study tract for all years. 

# this version of the script works for Demo tract only

# set working directory
setwd("C:/Users/jtsum/OneDrive/2019-0-spring/elasticity_project/rotation3/final_fileset/")

# libraries required
library(tidyverse)

FullLOI <- read_tsv("data/FullLOI.txt") %>% data.frame()
TerrsToKeep <- read_tsv("data/TerrsToKeep.txt") %>% data.frame()
terr.map <- read_tsv("data/TerrMap.txt") %>% data.frame()

### calculate area of outline of all included territories
StudyAreas <- filter(terr.map, 
                     TerrYr %in% TerrsToKeep$TerrYr) %>%
  group_by(Year) %>%
  dplyr::summarize(area = sum(ha))

# filter LOI for territories within study extent and for adults
LOI.clean <- filter(FullLOI, 
                    social.class != "juveniles",
                    TerrYr %in% TerrsToKeep$TerrYr)

# calculate the number of adult females and females breeders per year
densityCalc <- group_by(LOI.clean, Year) %>% 
  dplyr::summarize(adults = n(), 
                   breeders = sum(pair.status != "none"),
                   males = sum(Sex == "M", na.rm = TRUE),
                   females = sum(Sex == "F", na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(area = StudyAreas$area,
         adults.den = adults/area,
         breeders.den = breeders/area,
         male.den = males/area,
         female.den = females/area) %>%
  data.frame()

save(densityCalc,TerrsToKeep,
     file = "data/generated_data/densityCalcDemo.rdata")
