# Jeremy Summers
# This script generates model estimates for juvenile and helper survival rates
# using the same sex assignment as those used in the vital rates estimation.

# this version of the script uses only data from the Demo tract
# this version also generates both female only and male only models
# this version is designed for the simplified 4-stage model
library(glmmTMB)
library(tidyverse)
library(DHARMa)
library(performance)
setwd("C:/Users/jtsum/OneDrive/2019-0-spring/elasticity_project/rotation3/final_fileset/")

# load required files
## List of individuals
FullLOI <- read_tsv("data/FullLOI.txt") %>%
  data.frame()
## Territories to keep and density data
load("data/generated_data/densityCalcDemo.rdata")

# filter FullLOI for the study period
StudySiteLOI <- filter(FullLOI, TerrYr %in% TerrsToKeep$TerrYr,
                       Year %in% 1988:2021)

# convert identify immigrant individuals as individuals new to the study site
ImmFirstYear <- filter(StudySiteLOI,
                       !USFWS %in% filter(StudySiteLOI, 
                                          social.class == "juvenile")$USFWS) %>%
  group_by(USFWS) %>%
  dplyr::summarize(FirstYear = min(Year)) %>%
  # filter out 1988 as a possible first year, use already decided
  # immigrant status
  filter(FirstYear != 1988)

StudySiteLOI[paste(StudySiteLOI$USFWS, StudySiteLOI$Year) %in%
               paste(ImmFirstYear$USFWS, ImmFirstYear$FirstYear),
             "category"] <- "immigrant"

# check for individuals who left the study tract but later came back
IndvLastYear <- StudySiteLOI %>%
  group_by(USFWS) %>%
  dplyr::summarize(FirstYear = min(Year), LastYear = max(Year))

MissingYears <- filter(FullLOI, !TerrYr %in% TerrsToKeep$TerrYr,
                       Year %in% 1988:2021,
                       USFWS %in% StudySiteLOI$USFWS) %>%
  mutate(FirstYear = IndvLastYear[match(USFWS, IndvLastYear$USFWS),]$FirstYear,
         LastYear = IndvLastYear[match(USFWS, IndvLastYear$USFWS),]$LastYear) %>%
  filter(Year >= FirstYear, Year <= LastYear)

# there are 163 records for individuals who eventually came back to the study tract
# add these records back into StudySiteLOI
StudySiteLOI.clean <- bind_rows(StudySiteLOI, 
                                dplyr::select(MissingYears, -c("FirstYear",
                                                                      "LastYear")))

# convert survival metrics for individuals who leave the study area
StudySiteLOI.clean[paste(StudySiteLOI.clean$USFWS,
                         StudySiteLOI.clean$Year) %in%
                     paste(IndvLastYear$USFWS,
                           IndvLastYear$LastYear),
                   "survival"] <- 0

# correct vital rates dependent on survival
StudySiteLOI.clean[StudySiteLOI.clean$survival == 0, 
                   "transition"] <- NA

StudySiteLOI.clean[StudySiteLOI.clean$survival == 0, 
                   "divorce"] <- NA

StudySiteLOI.clean[StudySiteLOI.clean$survival == 0, 
                   "back_tansition"] <- NA

# load inbreeding values
ped.update <- read_tsv("data/pedInbr.txt") %>% data.frame()

# load environmental data
env.fac <- read_tsv("data/generated_data/env_var_updateDemo.txt") %>%
  data.frame() %>%
  filter(Year <= 2020) %>%
  # scale environmental factors
  scale() %>% data.frame() %>%
  mutate(Year.scale = Year, Year = 1988:2020)

# calculate pairwise inbreeding
library(ribd)
FSJped <- ped(ped.update$USFWS, ped.update$Dad,
              ped.update$Mom, ped.update$Sex)

Relatedness <- kinship(FSJped)

Relatedness.df <- reshape2::melt(Relatedness)
Relatedness.df$Var2 <- paste(Relatedness.df$Var2)
Relatedness.df$Var1 <- paste(Relatedness.df$Var1)

StudySiteLOI.clean$USFWS <- paste(StudySiteLOI.clean$USFWS)

LOI.fec <- left_join(StudySiteLOI.clean,  
                     Relatedness.df, 
                     by = c("USFWS" = "Var1",
                            "mate_USFWS" = "Var2")) %>%
  mutate(Year = Year - 1)

save(Relatedness, LOI.fec, 
     file = "data/generated_data/kinship_coef_Demo.rdata")

StudySiteLOI.clean$pedF <- ped.update[match(StudySiteLOI.clean$USFWS, 
                                            ped.update$USFWS), "pedF"]

remove(Relatedness.df, Relatedness, FSJped)

ped.update$USFWS <- paste(ped.update$USFWS)

# attach environmental data to census data
LOI.env <- left_join(dplyr::select(StudySiteLOI.clean, -pedF),
                     env.fac, by = "Year") %>%
  left_join(dplyr::select(ped.update, USFWS, pedF), by = "USFWS") %>%
  # scale inbreeding coef
  mutate(pedF = (pedF - 
           attr(scale(filter(ped.update, USFWS %in% StudySiteLOI.clean$USFWS)$pedF),
                "scaled:center"))/
             attr(scale(filter(ped.update, USFWS %in% StudySiteLOI.clean$USFWS)$pedF),
                  "scaled:scale"))

LOI.imm.env <- left_join(dplyr::select(mutate(StudySiteLOI.clean, Year = Year - 1),
                                       -pedF), env.fac, by = "Year") %>%
  filter(category == "immigrant", Year >= 1988)

# get pedF scale values
pedF.values <- filter(ped.update, USFWS %in% StudySiteLOI.clean$USFWS) %>%
  dplyr::summarize(mean = mean(pedF), sd = sd(pedF))

# get values to scale pair.Fped
pair.Fped.values <- filter(LOI.fec, Sex == "F", pair.status != "none",
                         Year >= 1988) %>%
  dplyr::select(USFWS, mate_USFWS, value) %>% distinct() %>%
  dplyr::summarize(mean = mean(value), sd = sd(value))

# create offset LOI
LOI.env.fec <- left_join(LOI.fec, env.fac, 
                         by = "Year") %>%
  left_join(ped.update, by = "USFWS") %>%
  filter(pair.status != "none")

# convert kinship matrix to a data frame
colnames(LOI.env.fec)[colnames(LOI.env.fec) == "value"] <- "pair.Fped"
colnames(LOI.env.fec)[colnames(LOI.env.fec) == "Sex.x"] <- "Sex"

# scale pair.Fped
LOI.env.fec <- mutate(LOI.env.fec,
                      pair.Fped = (pair.Fped - pair.Fped.values$mean)/
                        pair.Fped.values$sd)

##### SURVIVAL

# Check for associations between environmental variables and survival

LOI.env %>% dplyr::select(USFWS, Year, Sex, social.class,
                          survival,
                          acorns, burn, EQSOI, density,
                          pedF) %>%
  filter(Year < 2021) %>%
  pivot_longer(-c("USFWS", "Sex", "social.class",
                  "survival")) %>%
  ggplot(aes(x = value, y = survival, col = social.class)) +
  geom_smooth() + 
  facet_wrap(.~name, scales = "free_x")

# generate survival models

### Juvenile survival
Pj.df <- filter(LOI.env, social.class == "juvenile", Year < 2021)

Pj.model <- glmmTMB(survival ~ acorns + burn + Year.scale + EQSOI +
                    density + pedF +
                     (1|Year) + (1|Terr),
                  data = Pj.df, family = binomial(), REML = FALSE)

### Helper survival
Ph.df <- filter(LOI.env, social.class == "helper",
                Year < 2021)

Ph.model <- glmmTMB(survival ~ acorns + burn + Year.scale + EQSOI +
                    density + pedF + Sex + (1|Year) +
                    (1|USFWS) + (1|Terr), data = Ph.df, family = binomial(),
                    REML = FALSE)

### Breeder survival
Pb.df <- filter(LOI.env, pair.status != "none",
                Year < 2021)

Pb.model <- glmmTMB(survival ~ acorns + burn + Year.scale + EQSOI +
                    density + pedF + Sex + (1|Year) +
                    (1|USFWS) + (1|Terr), data = Pb.df, family = binomial(),
                    REML = FALSE)

# Pb.model without 1997

Pb.model.no1997 <- glmmTMB(survival ~ acorns + burn + Year.scale + EQSOI +
                             density + pedF + Sex + (1|Year) +
                             (1|USFWS) + (1|Terr), data = filter(Pb.df, Year != 1997), 
                           family = binomial(),
                           REML = FALSE)

#### TRANSITION TO BREEDER

LOI.env %>% dplyr::select(USFWS, Year, Sex, social.class,
                          transition, survival,
                          acorns, burn, EQSOI, density,
                          pedF) %>%
  filter(social.class %in% c("juvenile", "helper", Year < 2021), survival == 1) %>%
  pivot_longer(-c("USFWS", "Sex", "social.class",
                  "transition", "survival")) %>%
  ggplot(aes(x = value, y = transition, col = social.class)) +
  geom_smooth() + 
  facet_wrap(.~name, scales = "free_x")

# resident transition into breeder
B.model <- glmmTMB(transition ~ acorns + burn + Year.scale + EQSOI +
                     density + pedF + Sex + social.class + 
                     (1|Year) + (1|Terr) + (1|USFWS), 
                   data = filter(LOI.env, survival == 1, social.class != "breeder"),
                   family = binomial(), REML = FALSE)

#### IMMIGRATION

LOI.env %>% dplyr::select(USFWS, Year, Sex, social.class,
                          category,
                          acorns, burn, EQSOI, density) %>%
  filter(social.class %in% c("helper", "novice_breeder")) %>%
  pivot_longer(-c("USFWS", "Sex", "social.class",
                  "category")) %>%
  ggplot(aes(x = value, y = as.numeric(category == "immigrant"), 
             col = social.class)) +
  geom_smooth() + 
  facet_wrap(.~name, scales = "free_x")


# estimate the number of immigrants arriving per year
ImmRate <- filter(StudySiteLOI, category == "immigrant", Sex == "F",
                  social.class != "juvenile") %>%
  group_by(Year) %>%
  dplyr::summarize(I = n()) %>% 
  mutate(Year = Year - 1) %>%
  full_join(env.fac, by = "Year") %>%
  filter(Year <= 2020, Year >= 1988)

ImmRate[is.na(ImmRate)] <- 0

ImmRate.scale <- filter(ImmRate, !is.na(acorns)) %>%
  dplyr::select(Year, acorns, burn, EQSOI, density) %>%
  scale() %>% data.frame() %>%
  mutate(Year.scale = Year, Year = 1988:2020) %>%
  left_join(ImmRate, by = "Year", suffix = c("", ".y"))

IF.model <- glmmTMB(I ~ acorns + burn + Year.scale + EQSOI + density,
             data = ImmRate.scale, family = poisson(), REML = FALSE)

## using female only density metric improve model fit and finds extremely significant
## the effect of density. However, unless density is de-trended, the impact of year is
## lost. This is likely due to there being more unsexed adults in the earlier portion
## of the study period than the later. De-trending density recovers the impact of year.

# male immigration rate
ImmRateM <- filter(StudySiteLOI, category == "immigrant", Sex == "M",
                  social.class != "juvenile") %>%
  group_by(Year) %>%
  dplyr::summarize(I = n()) %>%
  mutate(Year = Year - 1) %>%
  full_join(env.fac, by = "Year") %>%
  filter(Year <= 2020, Year >= 1988)

ImmRateM[is.na(ImmRateM)] <- 0

ImmRateM.scale <- filter(ImmRateM, !is.na(acorns)) %>%
  dplyr::select(Year, acorns, burn, EQSOI, density) %>%
  scale() %>% data.frame() %>%
  mutate(Year.scale = Year, Year = 1988:2020) %>%
  left_join(ImmRateM, by = "Year", suffix = c("", ".y"))

IM.model <- glmmTMB(I ~ acorns + burn + Year.scale + EQSOI + density,
               data = ImmRateM.scale, family = poisson(), REML = FALSE)

#### Fecundity

F.df <- filter(LOI.env.fec, Year >= 1988, Sex == "F") %>%
  mutate(NestSuccess = as.numeric(fecundity > 0))

# create plot to test associations between env factors and fecundity
F.plot <- F.df %>% dplyr::select(fecundity, NestSuccess,
                                 pair.status, pair.Fped, acorns,
                                 burn, EQSOI, density, Year.scale) %>%
  pivot_longer(-c("fecundity", "NestSuccess", "pair.status"))

# colors for plot
colors <- RColorBrewer::brewer.pal(12, "Paired")[c(7, 8, 3, 4, 5, 
                                                   6, 1, 2, 9 , 10)]

fec.interact <- ggplot(F.plot, aes(x = value, y = fecundity, 
                                   col = pair.status)) +
  geom_jitter(alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_color_manual(values = colors[c(6, 8)],
                               labels = c("New", "Established")) +
  facet_wrap(.~name, scales = "free_x",
             labeller = labeller(name = c(acorns = "Acorns",
                                          burn = "Burned Area",
                                          density = "Pop. Density",
                                          EQSOI = "EQSOI",
                                          pair.Fped = "Relatedness",
                                          Year.scale = "Year"))) +
  labs(x = "Scaled Factor", y = "# of Offspring", col = "Pair Status") +
  theme_light() +
  theme(text = element_text(size = 20))

pdf(file = "figures/fec_interaction_plot_20231105.pdf",
    width = 14,
    height = 8.5)

fec.interact

dev.off()

# need to plot for pair.Fped differently
pair.Fped.fec <- filter(F.plot, NestSuccess == 1,
                        name == "pair.Fped") %>%
  mutate(Fped.cat = ifelse(value <= 0,
                           0,
                           ifelse(value <= 0.5,
                                  1,
                           ifelse(value <= 1,
                                  2,
                                  ifelse(value <= 2,
                                         3,
                                         ifelse(value <= 3,
                                                4,
                                                ifelse(value <= 4,
                                                       5,
                                                       6))))))) %>%
  group_by(Fped.cat, pair.status) %>%
  dplyr::summarize(fec = mean(fecundity),
                   sd_fec = sd(fecundity))

ggplot(pair.Fped.fec, aes(x = Fped.cat, y = fec, col = pair.status)) +
  geom_point(position = position_dodge(1)) +
  geom_errorbar(aes(ymin = fec - sd_fec, ymax = fec + sd_fec),
                position = position_dodge(1)) +
  labs(x = "Pair Relatedness Interval", y = "Mean # of Offspring")

# only major difference between old and new is 
F.model <- glmmTMB(fecundity ~ acorns + burn + Year.scale + EQSOI +
                      density + pair.Fped + pair.status + 
                      (1|Year) + (1|Terr), 
                   ziformula = ~ acorns + burn + Year.scale + EQSOI +
                     density + pair.Fped + pair.status + (1|Year) +
                     (1|Terr),
                    family = compois(),
                    REML = FALSE, data = F.df)

F.model.full <- glmmTMB(fecundity ~ acorns + burn + Year.scale + EQSOI +
                     density + pair.Fped + pair.status + 
                     (1|Year) + (1|Terr), 
                   ziformula = ~ acorns + burn + Year.scale + EQSOI +
                     density + pair.Fped*pair.status + (1|Year) +
                     (1|Terr),
                   family = compois(),
                   REML = FALSE, data = F.df)


# compare competing models
F.model.comparison <- compare_performance(F.model, F.model.full)
# delta AIC = 13.1, selecting full model

# for additional environmental factor, calculate average pairwise relatedness per year
AnnualRelatedness <- filter(LOI.env.fec, Sex == "F", pair.status != "none") %>%
  group_by(Year) %>%
  dplyr::summarize(pair.Fped = mean(pair.Fped)) %>%
  merge(LOI.env %>%
              group_by(Year, social.class) %>%
              dplyr::summarize(Fped = mean(pedF)) %>%
              pivot_wider(id_cols = Year, names_from = social.class,
                          values_from = Fped),
            by = "Year")


### Model Validation with DHARMa
set.seed(47)

# create residuals
Pj.sim <- simulateResiduals(Pj.model)
Ph.sim <- simulateResiduals(Ph.model)
Pb.sim <- simulateResiduals(Pb.model)
B.sim <- simulateResiduals(B.model)
F.sim <- simulateResiduals(F.model.full)
IM.sim <- simulateResiduals(IM.model)
IF.sim <- simulateResiduals(IF.model)

# test multicollinearity
multicollinearity.test <- multicollinearity(Pj.model) %>% data.frame() %>%
  mutate(model = "Pj") %>%
  bind_rows(mutate(data.frame(multicollinearity(Ph.model)), model = "Ph")) %>%
  bind_rows(mutate(data.frame(multicollinearity(Pb.model)), model = "Pb")) %>%
  bind_rows(mutate(data.frame(multicollinearity(B.model)), model = "B")) %>%
  bind_rows(mutate(data.frame(multicollinearity(F.model.full)), model = "F")) %>%
  bind_rows(mutate(data.frame(multicollinearity(IM.model)), model = "IM")) %>%
  bind_rows(mutate(data.frame(multicollinearity(IF.model)), model = "IF"))
# all VIF < 1.7, low correlation

# check residuals
Pj.residual.test <- testResiduals(Pj.sim)
# non-significant results for all test
Ph.residual.test <- testResiduals(Ph.sim)
# non-significant results for all test
Pb.residual.test <- testResiduals(Pb.sim)
# non-significant results for all test
B.residual.test <- testResiduals(B.sim)
# non-significant results for all test
F.residual.test <- testResiduals(F.sim)
# significant results: outliers: outliers at margins, p-value = 0.02549
## test influence of outliers
F.outliers <- outliers(F.sim)
F.model.nooutlier <- glmmTMB(fecundity ~ acorns + burn + Year.scale + EQSOI +
                               density + pair.Fped + pair.status + 
                               (1|Year) + (1|Terr), 
                             ziformula = ~ acorns + burn + Year.scale + EQSOI +
                               density + pair.Fped*pair.status + (1|Year) +
                               (1|Terr),
                             family = compois(),
                             REML = FALSE, data = F.df[-F.outliers,])
F.nooutlier.sim <- simulateResiduals(F.model.nooutlier)
F.nooutlier.residual.test <- testResiduals(F.nooutlier.sim)
## removing outliers causes minor change to fixed effect of density on
## conditional model and KS test becomes significant
## can stick with model without outliers removed
IM.residual.test <- testResiduals(IM.sim)
# non-significant results for all test
IF.residual.test <- testResiduals(IF.sim)
# non-significant results for all test

# check model performance
check_predictions(Pj.model)
check_predictions(Ph.model)
check_predictions(Pb.model)
check_predictions(B.model)
check_predictions(F.model.full)
check_predictions(IM.model)
check_predictions(IF.model)

# save models
save(LOI.env, LOI.env.fec, pedF.values, pair.Fped.values, AnnualRelatedness,
     Pj.model, Ph.model, Pb.model, Pb.model.no1997, B.model,
     IF.model, IM.model, F.model.nooutlier, F.model, F.model.full,
     file = "data/generated_data/vr_modelsDemo.rdata")

# save model validation results
save(multicollinearity.test, Pj.residual.test, Ph.residual.test,
     Pb.residual.test, B.residual.test, F.residual.test, F.nooutlier.residual.test,
     IM.residual.test, IF.residual.test, F.model.comparison,
     file = "data/generated_data/vr_modelsDemo_validation.rdata")
