# Jeremy Summers
# This script generates the supplementary plots needed for the LTRE paper
# October 2023

setwd("C:/Users/jtsum/OneDrive/2019-0-spring/elasticity_project/rotation3/final_fileset/")
library(tidyverse)
library(kableExtra)
library(glmmTMB)
library(Kendall)
library(performance)
library(sjPlot)

## Base vital rate table
# load required data
load("data/generated_data/densityCalcDemo.rdata")
FullLOI <- read_tsv("data/FullLOI.txt") %>%
  data.frame()
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
         LastYear = IndvLastYear[match(USFWS, IndvLastYear$USFWS),]$LastYear,
         social.class = "helper", pair.status = "none") %>%
  filter(Year >= FirstYear, Year <= LastYear)

dim(MissingYears)
unique(MissingYears$USFWS) %>% length()
# there are 121 records for individuals who eventually came back to the study tract
# add these records back into StudySiteLOI
# these records represent 67 unique individuals
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

# calculate base vital rates
vrObsSex <- StudySiteLOI.clean %>% group_by(Year) %>%
  dplyr::summarize(Pj = sum(social.class == "juvenile" & survival == 1)/
                    sum(social.class == "juvenile"),
                  Ph.M = sum(social.class == "helper" & survival == 1 & Sex == "M",
                            na.rm = TRUE)/
                    sum(social.class == "helper" & Sex == "M", na.rm = TRUE),
                  Ph.F = sum(social.class == "helper" & survival == 1 & Sex == "F",
                            na.rm = TRUE)/
                    sum(social.class == "helper" & Sex == "F", na.rm = TRUE),
                  Pb.M = sum(social.class == "breeder" & survival == 1 & Sex == "M",
                            na.rm = TRUE)/
                    sum(social.class == "breeder" & Sex == "M", na.rm = TRUE),
                  Pb.F = sum(social.class == "breeder" & survival == 1 & Sex == "F",
                            na.rm = TRUE)/
                    sum(social.class == "breeder" & Sex == "F",
                        na.rm = TRUE),
                  Bj = sum(social.class == "juvenile" & transition == 1, na.rm = TRUE)/
                    sum(social.class == "juvenile" & survival == 1),
                  Bh.F = sum(social.class == "helper" & transition == 1 & Sex == "F", 
                            na.rm = TRUE)/
                    sum(social.class == "helper" & survival == 1 & Sex == "F",
                        na.rm = TRUE),
                  Bh.M = sum(social.class == "helper" & transition == 1 & Sex == "M", 
                            na.rm = TRUE)/
                    sum(social.class == "helper" & survival == 1 & Sex == "M",
                        na.rm = TRUE),
                  Bb = sum(social.class == "breeder" & back_transition == 1,
                           na.rm = TRUE)/
                    sum(social.class == "breeder" & !is.na(back_transition), 
                        na.rm = TRUE),
                  Bi.M = sum(pair.status == "new" & category == "immigrant" & 
                              Sex == "M", na.rm = TRUE)/
                    sum(category == "immigrant" & Sex == "M", na.rm = TRUE),
                  Bi.F = sum(pair.status == "new" & category == "immigrant" &
                              Sex == "F", na.rm = TRUE)/
                    sum(category == "immigrant" & Sex == "F", na.rm = TRUE),
                  Dn = sum(pair.status == "new" & divorce == 1, na.rm = TRUE)/
                    sum(pair.status == "new"),
                  Do = sum(pair.status == "old" & divorce == 1, na.rm = TRUE)/
                    sum(pair.status == "old"),
                  Fn = sum(fecundity*(pair.status == "new"))/
                    sum(pair.status == "new"),
                  Fo = sum(fecundity*(pair.status == "old"))/
                    sum(pair.status == "old"),
                  I.M = sum(category == "immigrant" & Sex == "M", na.rm = TRUE)/
                    sum(Sex == "M", na.rm = TRUE),
                  I.F = sum(category == "immigrant" & Sex == "F", na.rm = TRUE)/
                    sum(Sex == "F", na.rm = TRUE))

vrObs <- StudySiteLOI.clean %>% group_by(Year) %>%
  dplyr::summarize(Pj = sum(social.class == "juvenile" & survival == 1)/
                     sum(social.class == "juvenile"),
                   Ph = sum(social.class == "helper" & survival == 1,
                              na.rm = TRUE)/
                     sum(social.class == "helper", na.rm = TRUE),
                   Pb = sum(social.class == "breeder" & survival == 1,
                              na.rm = TRUE)/
                     sum(social.class == "breeder", na.rm = TRUE),
                   Bj = sum(social.class == "juvenile" & transition == 1, na.rm = TRUE)/
                     sum(social.class == "juvenile" & survival == 1),
                   Bh = sum(social.class == "helper" & transition == 1, 
                              na.rm = TRUE)/
                     sum(social.class == "helper" & survival == 1,
                         na.rm = TRUE),
                   Bb = sum(social.class == "breeder" & back_transition == 1,
                            na.rm = TRUE)/
                     sum(social.class == "breeder" & !is.na(back_transition), 
                         na.rm = TRUE),
                   Bi = sum(pair.status == "new" & category == "immigrant",
                            na.rm = TRUE)/
                     sum(category == "immigrant", na.rm = TRUE),
                   Dn = sum(pair.status == "new" & divorce == 1, na.rm = TRUE)/
                     sum(pair.status == "new"),
                   Do = sum(pair.status == "old" & divorce == 1, na.rm = TRUE)/
                     sum(pair.status == "old"),
                   Fn = sum(fecundity*(pair.status == "new"))/
                     sum(pair.status == "new"),
                   Fo = sum(fecundity*(pair.status == "old"))/
                     sum(pair.status == "old"),
                   I = sum(category == "immigrant")/n())

# edit to line up vital rates in time
vrObs.clean <- vrObs %>% filter(Year < 2021) %>%
  mutate(Fn = vrObs$Fn[2:length(vrObs$Year)], Fo = vrObs$Fo[2:length(vrObs$Year)],
         I = vrObs$I[2:length(vrObs$Year)])

# create table with vital rate values
vrObsTable <- vrObs.clean %>%
  # round values
  mutate_at(colnames(vrObs.clean)[2:13], round, 3) %>%
  kbl(col.names = c("Year",
                    "Juvenile",
                    "Helper",
                    "Breeder",
                    "Juvenile",
                    "Helper", "Breeder", 
                    "Immigrant", "New Pair", "Established Pair", 
                    "New Pair", "Established Pair", " ")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T) %>%
  add_header_above(c(" " = 1, "Survival" = 3, "Transition" = 4, 
                     "Divorce" = 2,
                     "Fecundity" = 2, "Immigration" = 1))

save_kable(vrObsTable, file = "figures/vrObsTable_20231030.png",
           zoom = 5)
## correlations across vital rates

vrObs.clean.cor <- vrObs.clean %>% dplyr::select(-Ireal) %>%
  cor() %>% round(3) %>% data.frame()
vrObs.clean.cor[upper.tri(vrObs.clean.cor, diag = TRUE)] <- ""
rownames(vrObs.clean.cor) <- c()
vrObs.clean.cor <- data.frame(rowname = c("Year",
                                                       "Juvenile",
                                                       "Helper",
                                                       "Breeder",
                                                       "Juvenile",
                                                       "Helper", "Breeder", 
                                                       "Immigrant", "New Pair", 
                                                       "Established Pair", 
                                                       "New Pair", "Established Pair", 
                                                       "Immigration")) %>%
  bind_cols(vrObs.clean.cor)
vrObs.clean.cor <- vrObs.clean.cor[2:13,1:13]
rownames(vrObs.clean.cor) <- c()

vrObs.corTable <- vrObs.clean.cor %>%
  kbl(col.names = c("","Year",
                    "Juvenile",
                    "Helper",
                    "Breeder",
                    "Juvenile",
                    "Helper", "Breeder", 
                    "Immigrant", "New Pair", "Established Pair", 
                    "New Pair", "Established Pair")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  pack_rows("Survival", 1, 3) %>%
  pack_rows("Transition", 4, 7) %>%
  pack_rows("Divorce", 8, 9) %>%
  pack_rows("Fecundity", 10, 11) %>%
  column_spec(1, bold = T, border_right = T) %>%
  add_header_above(c(" " = 2, "Survival" = 3, "Transition" = 4, 
                     "Divorce" = 2,
                     "Fecundity" = 2))

save_kable(vrObs.corTable, file = "figures/final_paper/vrObscorTable_20231106.png",
           zoom = 5)

## Unsexed-filled vital rate tables
# load vr results
load("data/generated_data/vr_clean_F_4stageDemo.rdata")
vr.mat.F <- vr.mat
pop.vec.F <- pop.vec

load("data/generated_data/vr_clean_M_4stageDemo.rdata")
vr.mat.M <- vr.mat
pop.vec.M <- pop.vec

# generate means across replications
vr.mean.F <- apply(vr.mat.F, c(2, 3), mean)
vr.sd.F <- apply(vr.mat.F, c(2, 3), sd)
pop.mean.F <- apply(pop.vec.F, c(2, 3), mean)
pop.sd.F <- apply(pop.vec.F, c(2, 3), sd)
colnames(pop.mean.F) <- 1988:2021
colnames(pop.sd.F) <- 1988:2021

vr.mean.M <- apply(vr.mat.M, c(2, 3), mean)
vr.sd.M <- apply(vr.mat.M, c(2, 3), sd)
pop.mean.M <- apply(pop.vec.M, c(2, 3), mean)
pop.sd.M <- apply(pop.vec.M, c(2, 3), sd)
colnames(pop.mean.M) <- 1988:2021
colnames(pop.sd.M) <- 1988:2021

# convert vital rate results into printable table

vrTable.M <- matrix(paste(round(vr.mean.M, 3), round(vr.sd.M, 3),
                          sep = "±"),
                    nrow = nrow(vr.mean.M), ncol = ncol(vr.mean.M),
                    dimnames = list(rownames(vr.mean.M),colnames(vr.mean.M))) %>%
  t()


vrTable.M.print <- vrTable.M %>% data.frame() %>%
  tibble::rownames_to_column("Year") %>%
  dplyr::select(Year, Pj, Ph, Pb, Bj, Bh, Bb, Bi, Dn, Do, Fn, Fo, I) %>%
  # round values
  kbl(col.names = c("Year", "Juvenile",
                    "Helper",
                    "Breeder",
                    "Juvenile",
                    "Helper", "Breeder", 
                    "Immigrant", "New Pair", "Established Pair", 
                    "New Pair", "Established Pair", " ")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T) %>%
  add_header_above(c(" " = 1, "Survival" = 3, "Transition" = 4, 
                     "Divorce" = 2,
                     "Fecundity" = 2, "Immigration" = 1))

vrTable.F <- matrix(paste(round(vr.mean.F, 3), round(vr.sd.F, 3),
                          sep = "±"),
                    nrow = nrow(vr.mean.F), ncol = ncol(vr.mean.F),
                    dimnames = list(rownames(vr.mean.F),colnames(vr.mean.F))) %>%
  t()


vrTable.F.print <- vrTable.F %>% data.frame() %>%
  tibble::rownames_to_column("Year") %>%
  dplyr::select(Year, Pj, Ph, Pb, Bj, Bh, Bb, Bi, Dn, Do, Fn, Fo, I) %>%
  # round values
  kbl(col.names = c("Year", "Juvenile",
                    "Helper",
                    "Breeder",
                    "Juvenile",
                    "Helper", "Breeder", 
                    "Immigrant", "New Pair", "Established Pair", 
                    "New Pair", "Established Pair", " ")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T) %>%
  add_header_above(c(" " = 1, "Survival" = 3, "Transition" = 4, 
                     "Divorce" = 2,
                     "Fecundity" = 2, "Immigration" = 1))

save_kable(vrTable.M.print, file = "figures/vrMTable_20231102.png")
save_kable(vrTable.F.print, file = "figures/vrFTable_20231102.png")

# correlations between male and female vital rates
### calculate correlations between sexed vital rates
vr.sex.cor <- cor(t(vr.mean.F[c("Pj", "Ph", "Pb", "Bj", "Bh", "Bb",
                                "Bi", "Dn", "Do", "Fn", "Fo", "I"),]), 
                  t(vr.mean.M[c("Pj", "Ph", "Pb", "Bj", "Bh", "Bb",
                                "Bi", "Dn", "Do", "Fn", "Fo", "I"),])) %>% round(3)
# create rownames as a column
rownames(vr.sex.cor) <- c()
vr.sex.cor <- data.frame(rownames = c(
                                      "Juvenile",
                                      "Helper",
                                      "Breeder",
                                      "Juvenile",
                                      "Helper", "Breeder", 
                                      "Immigrant", "New Pair", 
                                      "Established Pair", 
                                      "New Pair", "Established Pair", 
                                      "Immigration")) %>%
  bind_cols(vr.sex.cor)

vr.sex.corTable <- data.frame(vr.sex.cor) %>%
  kbl(col.names = c("",
                    "Juvenile",
                    "Helper",
                    "Breeder",
                    "Juvenile",
                    "Helper", "Breeder", 
                    "Immigrant", "New Pair", "Established Pair", 
                    "New Pair", "Established Pair", "Immigration")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  pack_rows("Female-Only Model", 1, 12) %>%
  pack_rows("Survival", 1, 3) %>%
  pack_rows("Transition", 4, 7) %>%
  pack_rows("Divorce", 8, 9) %>%
  pack_rows("Fecundity", 10, 11) %>%
  column_spec(1, bold = T, border_right = T) %>%
  add_header_above(c(" " = 1, "Survival" = 3, "Transition" = 4, 
                     "Divorce" = 2,
                     "Fecundity" = 2, " " = 1)) %>%
  add_header_above(c(" " = 1, "Male-Only Model" = 12))

save_kable(vr.sex.corTable, file = "figures/vrSexcorTable_20231106.png",
           zoom = 5)

## population distribution
Obs.pop.vec <- StudySiteLOI.clean %>% group_by(Year) %>%
  dplyr::summarize(j = sum(social.class == "juvenile")/n(),
                   h = sum(social.class == "helper")/n(),
                   n = sum(pair.status == "new")/n(),
                   o = sum(pair.status == "old")/n())

Obs.pop.vecTable <- round(Obs.pop.vec, 3) %>%
  kbl(col.names = c("Year", "Juvenile", "Helper", "New Pair", "Established Pair")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T) %>%
  add_header_above(c(" " = 1, "% Population in Stage" = 4))

save_kable(Obs.pop.vecTable, file = "figures/popDist.png",
           zoom = 5)

# test for temporal trend in stage distribution
pop.vec.MannKendall <- apply(Obs.pop.vec[, c(2:5)], 2, MannKendall)
pop.vec.MannKendallTable <- data.frame(sc = c("Juvenile", "Helper", 
                                              "New Pair", "Established Pair"),
                                       tau = c(pop.vec.MannKendall$j$tau[[1]],
                                               pop.vec.MannKendall$h$tau[[1]],
                                               pop.vec.MannKendall$n$tau[[1]],
                                               pop.vec.MannKendall$o$tau[[1]]),
                                       p = c(pop.vec.MannKendall$j$sl[[1]],
                                             pop.vec.MannKendall$h$sl[[1]],
                                             pop.vec.MannKendall$n$sl[[1]],
                                             pop.vec.MannKendall$o$sl[[1]])) %>%
  mutate_at(c("tau", "p"), round, 3) %>%
  kbl(col.names = c("Social Class", "Tau", "2-sided p-value")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T)

save_kable(pop.vec.MannKendallTable, file = "figures/popDistMannKendall.png",
           zoom = 5)

## population growth rate table
load("data/generated_data/Demo_LTRE_results_20231106.rdata")

pop.summaryTable <- mutate_at(pop.summary, c("geoMean", "sd", 
                                             "popmean", "popsd"), round, 3) %>%
  mutate(geoMean = paste(geoMean, sd, sep = "±"),
         popMean = paste(popmean, popsd, sep = "±")) %>%
  dplyr::select(Model = Sex, geoMean, popMean) %>%
  kbl(col.names = c("Model", "Geometric Mean\nGrowth Rate",
                                      "Mean Annual\nPopulation Size")) %>%
    kable_classic(full_width = F, html_font = "Arial") %>%
    column_spec(1, bold = T, border_right = T)

# test for temporal trend in growth rate
lambda.obs %>% MannKendall()

save_kable(pop.summaryTable, file = "figures/popSummary.png",
           zoom = 5)

## Environmental correlation table
env.fac <- read_tsv("data/generated_data/env_var_updateDemo.txt") %>%
  data.frame()

load("data/generated_data/vr_modelsDemo.rdata")

env.fac.full <- right_join(env.fac, AnnualRelatedness, by = "Year")

env.cor <- cor(env.fac.full)
rownames(env.cor) <- c("Year", "Acorn Abundance", "Proportion Burned Area",
                       "EQSOI", "Population Density", "Pair Relatedness",
                       "Breeder Inbreeding", "Helper Inbreeding",
                       "Juvenile Inbreeding")

env.cor <- round(env.cor, 3)
env.cor[lower.tri(env.cor, diag = TRUE)] <- ""

env.cor.Table <- kbl(env.cor[1:8,2:9], col.names = rownames(env.cor[2:9,])) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T) %>%
  pack_rows("Mean Inbreeding/\nRelatendess Coef.", 6, 8) %>%
  add_header_above(c(" " = 5, "Mean Inbreeding/Relatendess Coef." = 4))

save_kable(env.cor.Table, file = "figures/env_corTable.png",
           zoom = 5)

## include VIF test

## Contribution of vital rates to population growth
# load final LTRE results
load("data/generated_data/Demo_LTRE_results.rdata")

# update labeling of vital rates for changing old to established
vrNameChange <- c("Fec", "Feq", "Dec", "Deq", "Ne")
names(vrNameChange) <- c("Foc", "Foq", "Doc", "Doq", "No")

# convert factor back to character
contrib.means$vr <- as.character(contrib.means$vr)
contrib.means$co_vr <- as.character(contrib.means$co_vr)

# update labels
contrib.means[contrib.means$vr %in% names(vrNameChange), 
              "vr"] <- contrib.means[contrib.means$vr %in% names(vrNameChange),]$vr %>%
  {vrNameChange[.]}

contrib.means[contrib.means$co_vr %in% names(vrNameChange), 
              "co_vr"] <- contrib.means[contrib.means$co_vr %in% 
                                          names(vrNameChange),]$co_vr %>%
  {vrNameChange[.]}

# convert back to factor
contrib.means$vr <- factor(contrib.means$vr,
                              levels = c("Pjc", "Phc", "Pbc", "Bjc", "Bhc", "Bbc",
                                         "Bic", "Dnc", "Dec", "Ic", "Fnc", "Fec",
                                         "Pjq", "Phq", "Pbq", "Bjq", "Bhq", "Bbq",
                                         "Biq", "Dnq", "Deq", "Iq", "Fnq", "Feq",
                                         "Nj","Nh","Nn","Ne"))

contrib.means$co_vr <- factor(contrib.means$co_vr,
                              levels = c("Pjc", "Phc", "Pbc", "Bjc", "Bhc", "Bbc",
                                         "Bic", "Dnc", "Dec", "Ic", "Fnc", "Fec",
                                         "Pjq", "Phq", "Pbq", "Bjq", "Bhq", "Bbq",
                                         "Biq", "Dnq", "Deq", "Iq", "Fnq", "Feq",
                                         "Nj","Nh","Nn","Ne"))

# produce a table that display the % contribution made by each vital rate interaction
limit <- max(abs(contrib.means$contrib)*100) * c(-1, 1)

raw_contrib_plot <- ggplot(contrib.means, aes(x = vr, y = co_vr)) +
  geom_tile(aes(fill = contrib*100), col = "gray40") +
  scale_fill_distiller(palette = "PRGn", limit = limit) +
  scale_x_discrete(labels = gsub("c|q", "", levels(contrib.means$vr))) +
  scale_y_discrete(labels = gsub("c|q", "", levels(contrib.means$vr))) +
  labs(fill = "% Contribution\nto Var. in Pop. Growth") +
  coord_cartesian(ylim = c(1, 28), xlim = c(1, 28), clip = "off") +
  facet_wrap(.~sex, labeller = labeller(sex = c(F = "Female-Only Model",
                                                M = "Male-Only Model"))) +
  labs(x = "", y = "") +
  geom_text(aes(label = ifelse(abs(round(contrib*100, 1)) > 0,
                               round(contrib*100, 1), "")),
            size = 2) +
  geom_text(x = 6, y = -0.9, label = "Current", col = "royalblue3") +
  geom_text(x = 18, y = -0.9, label = "Past", col = "tomato3") +
  geom_text(data = filter(contrib.means, sex == "F"),
            x = -1.5, y = 6, label = "Current", col = "royalblue3", angle = 90) +
  geom_text(data = filter(contrib.means, sex == "F"),
            x = -1.5, y = 18, label = "Past", col = "tomato3", angle = 90) +
  geom_segment(aes(x = 0.5, xend = 12.5, y = -0.5, yend = -0.5),
               col = "royalblue3", linewidth = 1) +
  geom_segment(aes(x = 12.5, xend = 24.5, y = -0.5, yend = -0.5),
               col = "tomato3", linewidth = 1) +
  geom_segment(data = filter(contrib.means, sex == "F"),
    aes(x = -1, xend = -1, y = 0.5, yend = 12.5),
               col = "royalblue3", linewidth = 1) +
  geom_segment(data = filter(contrib.means, sex == "F"),
    aes(x = -1, xend = -1, y = 12.5, yend = 24.5),
               col = "tomato3", linewidth = 1) +
  theme_light() +
  theme(strip.text = element_text(size = 20))

pdf(file = "figures/raw_contrib_plot_20231205.pdf",
    width = 14,
    height = 8.5)

raw_contrib_plot

dev.off()

## Tables of linear model results

# create functions to extract effects from models
get_est <- function(model){
  est <- coef(summary(model))$cond %>% data.frame() %>%
    mutate(ub = Estimate + 1.96*Std..Error,
           lb = Estimate - 1.96*Std..Error) %>%
    rownames_to_column("factor") %>%
    dplyr::select(factor, est = Estimate, se = Std..Error, p_value = Pr...z..,
                  ub, lb)
  return(est)
}

get_est_zi <- function(model){
  est <- coef(summary(model))$zi %>% data.frame() %>%
    mutate(ub = Estimate + 1.96*Std..Error,
           lb = Estimate - 1.96*Std..Error) %>%
    rownames_to_column("factor") %>%
    dplyr::select(factor, est = Estimate, se = Std..Error, p_value = Pr...z..,
                  ub, lb)
  return(est)
}

# create key to adjust factor
envFactorKey <- c("Intercept", "Acorn Abund.", "Burned Area", "Year", "EQSOI",
                  "Pop. Density", "Inbreeding Coef.", "Sex",
                  "Social Class", "Pair Relatedness",
                  "Pair Status", "Pair Relatedness:Pair Status")
names(envFactorKey) <- c("(Intercept)", "acorns", "burn", "Year.scale",
                         "EQSOI", "density", "pedF",
                         "SexM", "social.classjuvenile",
                         "pair.Fped", "pair.statusold", "pair.Fped:pair.statusold")

# gather estimate values
model_est <- get_est(Pj.model) %>% mutate(model = "Pj") %>%
  bind_rows(mutate(get_est(Ph.model), model = "Ph")) %>%
  bind_rows(mutate(get_est(Pb.model), model = "Pb")) %>%
  bind_rows(mutate(get_est(B.model), model = "B")) %>%
  bind_rows(mutate(get_est(IM.model), model = "IM")) %>%
  bind_rows(mutate(get_est(IF.model), model = "IF")) %>%
  bind_rows(mutate(get_est(F.model.full), model = "F2")) %>%
  bind_rows(mutate(get_est_zi(F.model.full), model = "F1")) %>%
  # adjust factor names
  mutate(fac = envFactorKey[factor],
         Sex = ifelse(grepl("Sex", factor), "M",
                      ifelse(factor == "density" & model == "I",
                             "F", "both")),
         pair = ifelse(!grepl("F", model), "none",
                       ifelse(factor == "pair.Fped", "new",
                              ifelse(factor == "pair.Fped:pair.statusold", "old",
                                     "both"))))

# create table for survival models

surv.modelTable <- filter(model_est, grepl("P", model)) %>%
  dplyr::select(fac, est, se, p_value) %>%
  mutate_at(c("est", "se", "p_value"), round, 3) %>%
  mutate(p_value = ifelse(p_value != 0, paste(p_value),
                          "< 0.001")) %>%
  kbl(col.names = c("Factor", "Estimate", "Standard Error",
                    "P-Value")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  pack_rows("Juvenile Survival", 1, 7) %>%
  pack_rows("Helper Survival", 8, 15) %>%
  pack_rows("Breeder Survival", 16, 23) %>%
  column_spec(1, bold = T, border_right = T)

save_kable(surv.modelTable, file = "figures/SurvModelTable_20231210.png",
           zoom = 5)

# create table for transition models
B.modelTable <- filter(model_est, grepl("B", model)) %>%
  dplyr::select(fac, est, se, p_value) %>%
  mutate_at(c("est", "se", "p_value"), round, 3) %>%
  mutate(p_value = ifelse(p_value != 0, paste(p_value),
                          "< 0.001")) %>%
  kbl(col.names = c("Factor", "Estimate", "Standard Error",
                    "P-Value")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T)

save_kable(B.modelTable, file = "figures/BModelTable_20231210.png",
           zoom = 5)

# create table for fecundity models
fec.modelTable <- filter(model_est, model %in% c("F1", "F2")) %>%
  dplyr::select(fac, est, se, p_value) %>%
  mutate_at(c("est", "se", "p_value"), round, 3) %>%
  mutate(p_value = ifelse(p_value != 0, paste(p_value),
                          "< 0.001")) %>%
  kbl(col.names = c("Factor", "Estimate", "Standard Error",
                    "P-Value")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  pack_rows("Conditional Model", 1, 8) %>%
  pack_rows("Zero-Inflated Model", 9, 17) %>%
  column_spec(1, bold = T, border_right = T)

save_kable(fec.modelTable, file = "figures/FecModelTable_20231210.png",
           zoom = 5)

# create table for immigration models
im.modelTable <- filter(model_est, grepl("I", model)) %>%
  dplyr::select(fac, est, se, p_value) %>%
  mutate_at(c("est", "se", "p_value"), round, 3) %>%
  mutate(p_value = ifelse(p_value != 0, paste(p_value),
                          "< 0.001")) %>%
  kbl(col.names = c("Factor", "Estimate", "Standard Error",
                    "P-Value")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  pack_rows("Male Immigration", 1, 6) %>%
  pack_rows("Female Immigration", 7, 12) %>%
  column_spec(1, bold = T, border_right = T)

save_kable(im.modelTable, file = "figures/ImModelTable_20231210.png",
           zoom = 5)

# create table of model AIC
fec.model.compare <- data.frame(model = c("Base Fecundity Model",
                                          "Fecundity Model with Interaction Term"),
                                AIC = c(summary(F.model)$AICtab[1],
                                        summary(F.model.full)$AICtab[1])) %>%
  mutate(deltaAIC = AIC - min(AIC))

fec.model.compareTable <- mutate_at(fec.model.compare, 
                                    c("AIC", "deltaAIC"), round, 3) %>%
  kbl(col.names = c("Model", "AIC", "deltaAIC")) %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  column_spec(1, bold = T, border_right = T)

save_kable(fec.model.compareTable, 
           file = "figures/fecModelCompareTable_20231210.png",
           zoom = 5)

# survival model table
tab_model(Pj.model, Ph.model, Pb.model,
          dv.labels = c("Juvenile Survival",
                        "Helper Survival",
                        "Breeder Survival"),
          CSS = list(css.table = '+font-family: Arial;')) %>%
  save_kable(file = "figures/SurvModelTable_20231210.png")

# transition model table
tab_model(B.model,
          dv.labels = c("Prob. to Pair"),
          CSS = list(css.table = '+font-family: Arial;')) %>%
  save_kable(file = "figures/BModelTable_20231210.png")

# fecundity model table
F.Table <- tab_model(F.model.full,
          dv.labels = c("Fecundity"),
          CSS = list(css.table = '+font-family: Arial;')) %>%
  save_kable(file = "figures/fecModelTable_20231210.png")

# immigration model table
Im.Table <- tab_model(IM.model, IF.model,
                      dv.labels = c("Male Immigration", "Female Immigrataion"),
                      CSS = list(css.table = '+font-family: Arial;')) %>%
  save_kable(file = "figures/ImModelTable_20231210.png")
