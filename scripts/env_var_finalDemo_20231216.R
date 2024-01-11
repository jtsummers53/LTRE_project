# Jeremy Summers
# This script generates the environmental variables used for my LTRE
# March 2023

library(tidyverse)

# Acorn Data
acorns <- read_tsv("data/acorns_update.txt") %>% data.frame()

# Density Data
load("data/generated_data/densityCalcDemo.rdata")

# Burn Area Data
## load burn data
FireArea <- read_tsv("data/TerrYrBurnArea.txt") %>% data.frame()
burnArea <- filter(FireArea, TerrYr %in% TerrsToKeep$TerrYr) %>%
  group_by(Year) %>%
  dplyr::summarize(tot_area = sum(area_ha))

burnArea$TSF5_area <- filter(FireArea, TSF_yr >= 2, TSF_yr <= 5,
                             TerrYr %in% TerrsToKeep$TerrYr) %>%
  group_by(Year) %>% {dplyr::summarize(., area = sum(area_ha))$area}

burnArea$TSF9_area <- filter(FireArea, TSF_yr >= 2, TSF_yr <= 9,
                             TerrYr %in% TerrsToKeep$TerrYr) %>%
  group_by(Year) %>% {dplyr::summarize(., area = sum(area_ha))$area}

burnArea <- mutate(burnArea, area5 = TSF5_area/tot_area, area9 = TSF9_area/tot_area) %>%
  data.frame()

# EQSOI Data
EQSOI <- read.table(file = "data/reqsoi_update.txt") %>% data.frame()
colnames(EQSOI) <- c("Year", 1:12)
EQSOI.long <- pivot_longer(EQSOI, !Year, names_to = "month",
                           values_to = "EQSOI") %>%
  mutate(CensusYear = ifelse(month < 6, Year - 1, Year)) %>%
  # filter for dry season
  filter(month %in% c(1:3, 10:12)) %>%
  group_by(CensusYear) %>%
  filter(CensusYear <= 2021) %>%
  dplyr::summarize(EQSOIdry = mean(EQSOI))

# calculate EQSOI for full year
EQSOI.long <- pivot_longer(EQSOI, !Year, names_to = "month",
                           values_to = "EQSOI") %>%
  mutate(CensusYear = ifelse(month < 6, Year - 1, Year)) %>%
  group_by(CensusYear) %>%
  dplyr::summarize(EQSOI = mean(EQSOI)) %>%
  filter(CensusYear <= 2021) %>%
                     full_join(EQSOI.long)
EQSOI.long <- data.frame(EQSOI.long)

# load vital rates to test for correlations with lagged environmental factors
load("data/generated_data/vr_clean_F_4stageDemo.rdata")

vrAnnualDemo.F <- apply(vr.mat, c(2, 3), mean)

load("data/generated_data/vr_clean_M_4stageDemo.rdata")

vrAnnualDemo.M <- apply(vr.mat, c(2, 3), mean)


## acorns correlation
cor(acorns[1:33, 4], t(vrAnnualDemo.F))
cor(acorns[1:33, 4], t(vrAnnualDemo.M))

cor(acorns[1:32, 4], t(vrAnnualDemo.F[, 2:33]))
cor(acorns[1:32, 4], t(vrAnnualDemo.M[, 2:33]))

### test acorn autocorrelation
cor.test(acorns[1:32, 4], acorns[2:33, 4])
#### cor = -0.01091039, p-value = 0.9527

### test for correlation between past acorn and breeder survival
cor.test(acorns[1:32, 4], vrAnnualDemo.F["Pb", 2:33])
#### cor = -0.4447469, p-value = 0.01076

### test for correlation between past acorn and breeder survival
cor.test(acorns[1:32, 4], vrAnnualDemo.M["Pb", 2:33])
#### cor = -0.3691391, p-value = 0.0376


## density correlation
cor(densityCalc[1:33, 7], t(vrAnnualDemo.F))
cor(densityCalc[1:33, 7], t(vrAnnualDemo.M))

cor(densityCalc[1:32, 7], t(vrAnnualDemo.F[, 2:33]))
cor(densityCalc[1:32, 7], t(vrAnnualDemo.M[, 2:33]))

### test density autocorrelation
cor.test(densityCalc[1:32, 7], densityCalc[2:33, 7])
#### cor = -0.06010695, p-value = 0.7438

### test for correlation between past density and juvenile survival
cor.test(densityCalc[1:32, 7], vrAnnualDemo.F["Pj", 2:33])
#### cor = -0.3239353, p-value = 0.0705

### test for correlation between past density and juvenile survival
cor.test(densityCalc[1:32, 7], vrAnnualDemo.M["Pj", 2:33])
#### cor = -0.39815, p-value = 0.02401

## burn correlation
cor(burnArea[1:33, 6], t(vrAnnualDemo.F))
cor(burnArea[1:33, 6], t(vrAnnualDemo.M))

cor(burnArea[1:32, 6], t(vrAnnualDemo.F[, 2:33]))
cor(burnArea[1:32, 6], t(vrAnnualDemo.M[, 2:33]))

### no major correlations between past burn and vital rates

## EQSOI correlation
cor(EQSOI.long[41:73, 2], t(vrAnnualDemo.F))
cor(EQSOI.long[41:73, 2], t(vrAnnualDemo.M))

cor(EQSOI.long[40:72, 2], t(vrAnnualDemo.F))
cor(EQSOI.long[40:72, 2], t(vrAnnualDemo.M))

### test density autocorrelation
cor.test(EQSOI.long[1:73, 2], EQSOI.long[2:74, 2])
#### cor = 0.3927946, p-value = 0.0005874

### test for correlation between past EQSOI and fecundity
cor.test(EQSOI.long[40:72, 2], vrAnnualDemo.F["Fn",])
#### cor = -0.3254959, p-value = 0.06454

### test for correlation between past density and juvenile survival
cor.test(EQSOI.long[40:72, 2], vrAnnualDemo.F["Fo",])
#### cor = -0.2196853, p-value = 0.2193

## There are no significant time-lagged correlations for EQSOI

# combine data
env.var.final <- data.frame(Year = 1988:2021, 
                            acorns = acorns$acorns.mean[1:34],
                            burn = burnArea$area9,
                            EQSOI = EQSOI.long$EQSOI[41:74],
                            density = densityCalc$adults.den)

write_tsv(env.var.final, 
          file = "data/generated_data/env_var_updateDemo.txt")
