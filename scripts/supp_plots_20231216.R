# Jeremy Summers
# This script generates the supplementary plots needed for the LTRE paper
# October 2023

library(tidyverse)
library(kableExtra)
library(glmmTMB)
library(Kendall)
library(performance)
library(Hmisc)
library(modelsummary)
library(flextable)

options(knitr.kable.NA = '')

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

## population distribution
Obs.pop.vec <- StudySiteLOI.clean %>% group_by(Year) %>%
  dplyr::summarize(Nj = sum(social.class == "juvenile")/n(),
                   Nh = sum(social.class == "helper")/n(),
                   Nn = sum(pair.status == "new")/n(),
                   Ne = sum(pair.status == "old")/n())

# create plots of vital rates

# re-format vital rate data frame
vrObs.clean.long <- pivot_longer(vrObs.clean, -Year) %>%
  # add identifiers of vital rate and social class
  mutate(vr = substr(name, 1, 1) %>% factor(levels = c("F", "P",
                                                           "I", "B", "D", "N")),
         class = ifelse(grepl("j", name), "j",
                        ifelse(grepl("h", name), "h",
                               ifelse(grepl("b", name), "b",
                                      ifelse(grepl("i|I", name), "I",
                                             ifelse(grepl("n", name), "n",
                                                    ifelse(grepl("o", name), "e",
                                                           NA)))))) %>%
           factor(levels = c("j", "h", "b", "I", "n", "e")))

# re-format the population distribution data
Obs.pop.vec.long <- pivot_longer(Obs.pop.vec, -Year) %>%
  mutate(vr = factor("N", levels = c("F", "P",
                                     "I", "B", "D", "N")),
         class = substr(name, 2, 2) %>%
           factor(levels = c("j", "h", "b", "I", "n", "e")))

vrObs.clean.long <- bind_rows(vrObs.clean.long, Obs.pop.vec.long)

vr_plot <- ggplot(vrObs.clean.long, aes(x = Year, y = value, col = class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(labels = c("Juvenile",
                                "Helper",
                                "Breeder",
                                "Immigrant",
                                "New Pair",
                                "Established Pair"),
                     values = c("#FF7F00",
                                "#33A02C",
                                "#CAB2D6",
                                "#7570B3",
                                "#E31A1C",
                                "#1F78B4")) +
  facet_grid(vr~., scales = "free_y",
             labeller = labeller(vr = c(F = "Fecundity",
                                        P = "Survival",
                                        I = "Immigration",
                                        B = "Prob. to Pair",
                                        D = "Divorce",
                                        N = "Stage Dist."))) +
  theme_light() +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.title = element_blank())

ggsave("figures/vr_plot.png",
       vr_plot,
       width = 6,
       height = 7,
       units = "in",
       dpi = 700)

ggsave("figures/vr_plot.pdf",
       vr_plot,
       width = 6,
       height = 7,
       units = "in")

## correlations across vital rates

vrObs.clean.cor <- vrObs.clean %>%
  cor() %>% round(3) %>% data.frame()

vrObs.clean.cor <- rcorr(as.matrix(vrObs.clean), type = "spearman")
vr.key <- c("Year", "Juvenile", "Helper",
            "Breeder", "Juvenile",
            "Helper", "Breeder",
            "Immigrant", "New Pair",
            "Est. Pair", "New Pair", "Est. Pair",
            "Immigration")
names(vr.key) <- colnames(vrObs.clean)

# convert correlation matrix into a data frame for plotting
vrObs.clean.cor.df <- reshape2::melt(vrObs.clean.cor$r) %>%
  bind_cols(reshape2::melt(vrObs.clean.cor$P) %>% select(p = value)) %>%
  # add stars to mark significance
  mutate(stars = ifelse(is.na(p), NA,
                        ifelse(p > 0.05, "",
                               ifelse(p > 0.01, "*",
                                      ifelse( p > 0.001, "**", "***")))))

# split dataset into lower and upper triangle, and diagonal
vrObs.clean.cor.low <- filter(vrObs.clean.cor.df[lower.tri(vrObs.clean.cor$r),], 
                              Var1 != Var2)
vrObs.clean.cor.diag <- filter(vrObs.clean.cor.df, Var1 == Var2)

vrObs.clean.cor.labels <- as.character(unique(vrObs.clean.cor.df$Var2))

vrObs.clean.cor.plot <- ggplot(vrObs.clean.cor.low, 
                               aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  # add text with significance stars
  geom_text(aes(label = paste0(round(value, 2), stars)),
            size = 3) +
  scale_fill_distiller(palette = "PRGn", limit = c(-1, 1)) +
  # re-order axes
  scale_y_discrete(limits = vrObs.clean.cor.labels[
    1:(length(vrObs.clean.cor.labels) - 1)],
                   labels = vr.key) +
  scale_x_discrete(limits = vrObs.clean.cor.labels[
    length(vrObs.clean.cor.labels):2],
                   labels = vr.key) +
  coord_cartesian(ylim = c(1, 12), xlim = c(1, 12), clip = "off") +
   # add vital rate labels for x axis
  geom_segment(aes(x = 1.5, xend = 3.5, y = -2.3, yend = -2.3),
               col = "#D95F02", linewidth = 2) +
   geom_text(x = 2.5, y = -1.95, label = "Fecundity", col = "#D95F02",
             size = 5) +
  geom_segment(aes(x = 3.5, xend = 5.5, y = -2.3, yend = -2.3),
               col = "#A6761D", linewidth = 2) +
   geom_text(x = 4.5, y = -1.95, label = "Divorce", col = "#A6761D",
             size = 5) +
  geom_segment(aes(x = 5.5, xend = 9.5, y = -2.3, yend = -2.3),
               col = "#1B9E77", linewidth = 2) +
   geom_text(x = 7.5, y = -1.95, label = "Prob. to Pair", col = "#1B9E77",
             size = 5) +
  geom_segment(aes(x = 9.5, xend = 12.5, y = -2.3, yend = -2.3),
               col = "#E7298A", linewidth = 2) +
   geom_text(x = 11, y = -1.95, label = "Survival", col = "#E7298A",
             size = 5) +
   # add vital rate labels for y axis
   geom_segment(aes(x = -1.75, xend = -1.75, y = 10.5, yend = 12.5),
                col = "#D95F02", linewidth = 2) +
   geom_text(x = -2.1, y = 11.5, label = "Fecundity", col = "#D95F02",
             size = 5, angle = 90) +
   geom_segment(aes(x = -1.75, xend = -1.75, y = 8.5, yend = 10.5),
                col = "#A6761D", linewidth = 2) +
   geom_text(x = -2.1, y = 9.5, label = "Divorce", col = "#A6761D",
             size = 5, angle = 90) +
   geom_segment(aes(x = -1.75, xend = -1.75, y = 4.5, yend = 8.5),
                col = "#1B9E77", linewidth = 2) +
   geom_text(x = -2.1, y = 6.5, label = "Prob. to Pair", col = "#1B9E77",
             size = 5, angle = 90) +
   geom_segment(aes(x = -1.75, xend = -1.75, y = 1.5, yend = 4.5),
                col = "#E7298A", linewidth = 2) +
   geom_text(x = -2.1, y = 3, label = "Survival", col = "#E7298A",
             size = 5, angle = 90) +
   labs(y = "", x = "") +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("figures/vr_cor_plot.png",
       vrObs.clean.cor.plot,
       width = 6.5,
       height = 6.5,
       units = "in",
       dpi = 700)

ggsave("figures/vr_cor_plot.pdf",
       vrObs.clean.cor.plot,
       width = 6.5,
       height = 6.5,
       units = "in")

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

# correlations between male and female vital rates
### calculate correlations between sexed vital rates
vr.sex.cor <- rcorr(t(vr.mean.F[c("Pj", "Ph", "Pb", "Bj", "Bh", "Bb",
                                "Bi", "Dn", "Do", "Fn", "Fo", "I"),]), 
                  t(vr.mean.M[c("Pj", "Ph", "Pb", "Bj", "Bh", "Bb",
                                "Bi", "Dn", "Do", "Fn", "Fo", "I"),]),
                  type = "spearman")

# convert correlations into plot-able data frame
# select on quadrant of correlations that captures relationship between
# sex-based models
vr.sex.cor.df <- reshape2::melt(vr.sex.cor$r[1:12, 13:24]) %>%
  bind_cols(reshape2::melt(vr.sex.cor$P[1:12, 13:24]) %>% select(p = value)) %>%
  # add stars to mark significance
  mutate(stars = ifelse(is.na(p), NA,
                        ifelse(p > 0.05, "",
                               ifelse(p > 0.01, "*",
                                      ifelse( p > 0.001, "**", "***")))))


vr.sex.cor.plot <- ggplot(vr.sex.cor.df, 
                               aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  # add text with significance stars
  geom_text(aes(label = paste0(round(value, 2), stars)),
            size = 3) +
  scale_fill_distiller(palette = "PRGn", limit = c(-1, 1.2)) +
  # re-order axes
  scale_y_discrete(limits = unique(vr.sex.cor.df$Var1)[c(12, 1:11)],
    labels = vr.key) +
  scale_x_discrete(limits = unique(vr.sex.cor.df$Var1)[12:1],
    labels = vr.key) +
  coord_cartesian(ylim = c(1, 12), xlim = c(1, 12), clip = "off") +
  # add vital rate labels for x axis
  geom_segment(aes(x = 1.5, xend = 3.5, y = -2.5, yend = -2.5),
               col = "#D95F02", linewidth = 2) +
  geom_text(x = 2.5, y = -2.15, label = "Fecundity", col = "#D95F02",
            size = 5) +
  geom_segment(aes(x = 3.5, xend = 5.5, y = -2.5, yend = -2.5),
               col = "#A6761D", linewidth = 2) +
  geom_text(x = 4.5, y = -2.15, label = "Divorce", col = "#A6761D",
            size = 5) +
  geom_segment(aes(x = 5.5, xend = 9.5, y = -2.5, yend = -2.5),
               col = "#1B9E77", linewidth = 2) +
  geom_text(x = 7.5, y = -2.15, label = "Prob. to Pair", col = "#1B9E77",
            size = 5) +
  geom_segment(aes(x = 9.5, xend = 12.5, y = -2.5, yend = -2.5),
               col = "#E7298A", linewidth = 2) +
  geom_text(x = 11, y = -2.15, label = "Survival", col = "#E7298A",
            size = 5) +
  # add vital rate labels for y axis
  geom_segment(aes(x = -2.15, xend = -2.15, y = 10.5, yend = 12.5),
               col = "#D95F02", linewidth = 2) +
  geom_text(x = -2.5, y = 11.5, label = "Fecundity", col = "#D95F02",
            size = 5, angle = 90) +
  geom_segment(aes(x = -2.15, xend = -2.15, y = 8.5, yend = 10.5),
               col = "#A6761D", linewidth = 2) +
  geom_text(x = -2.5, y = 9.5, label = "Divorce", col = "#A6761D",
            size = 5, angle = 90) +
  geom_segment(aes(x = -2.15, xend = -2.15, y = 4.5, yend = 8.5),
               col = "#1B9E77", linewidth = 2) +
  geom_text(x = -2.5, y = 6.5, label = "Prob. to Pair", col = "#1B9E77",
            size = 5, angle = 90) +
  geom_segment(aes(x = -2.15, xend = -2.15, y = 1.5, yend = 4.5),
               col = "#E7298A", linewidth = 2) +
  geom_text(x = -2.5, y = 3, label = "Survival", col = "#E7298A",
            size = 5, angle = 90) +
  labs(y = "Male-Only Model", x = "Female-Only Model") +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


ggsave("figures/vr_sex_cor_plot.png",
       vr.sex.cor.plot,
       width = 6.5,
       height = 6.5,
       units = "in",
       dpi = 700)

ggsave("figures/vr_sex_cor_plot.pdf",
       vr.sex.cor.plot,
       width = 6.5,
       height = 6.5,
       units = "in")

## population growth rate table
load("data/generated_data/Demo_LTRE_results_20231106.rdata")

pop.summary[pop.summary$Sex == "Raw Observations", 
            "Sex"] <- "Both-sex Observations" 

pop.summary <- arrange(pop.summary, Sex)

# create table in docx format
pop.summary %>% mutate(geoMean = paste0(round(geoMean, 3), " (", round(sd, 3), ")"),
                       popmean = paste0(round(popmean, 1), 
                                        " (", round(popsd, 1), ")")) %>%
  dplyr::select(Model = Sex, "Geometric Mean Growth Rate" = geoMean, 
                "Mean Population Size" = popmean) %>%
  flextable() %>%
  width(j = 1, 1.75) %>%
  width(j = 2, 2.3) %>%
  width(j = 3, 1.95) %>%
  theme_vanilla() %>%
  save_as_docx(path = "figures/popSummary.docx")

# test for temporal trend in growth rate
lambda.obs %>% MannKendall()

## Environmental correlation table
env.fac <- read_tsv("data/generated_data/env_var_updateDemo.txt") %>%
  data.frame()

load("data/generated_data/vr_modelsDemo.rdata")

# de-scale annual relatedness values
AnnualRelatedness[, c(3:5)] <- (AnnualRelatedness[, c(3:5)]*pedF.values$sd) +
  pedF.values$mean

AnnualRelatedness[, 2] <- (AnnualRelatedness[, 2]*pair.Fped.values$sd) +
  pair.Fped.values$mean

env.fac.full <- right_join(env.fac, AnnualRelatedness, by = "Year")

# calculate correlation between environmental values
env.cor <- rcorr(as.matrix(env.fac.full), type = "spearman")
env.key <- c("Year", "Acorn Abundance", "Proportion Burned Area",
                       "EQSOI", "Population Density", "Pair Relatedness",
                       "Breeder Inbreeding", "Helper Inbreeding",
                       "Juvenile Inbreeding")

names(env.key) <- rownames(env.cor$r)

# convert correlation matrix into a data frame for plotting
env.cor.df <- reshape2::melt(env.cor$r) %>%
  bind_cols(reshape2::melt(env.cor$P) %>% select(p = value)) %>%
  # add stars to mark significance
  mutate(stars = ifelse(is.na(p), NA,
                        ifelse(p > 0.05, "",
                               ifelse(p > 0.01, "*",
                                      ifelse( p > 0.001, "**", "***")))))

# split dataset into lower and upper triangle, and diagonal
env.cor.up <- filter(env.cor.df[upper.tri(env.cor$r),], Var1 != Var2)
env.cor.low <- filter(env.cor.df[lower.tri(env.cor$r),], Var1 != Var2)
env.cor.diag <- filter(env.cor.df, Var1 == Var2)

env.cor.labels <- as.character(unique(env.cor.df$Var2))

env.cor.plot <- ggplot(env.cor.low, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  # add text with significance stars
  geom_text(aes(label = paste0(round(value, 2), stars)),
            size = 3) +
  scale_fill_distiller(palette = "PRGn", limit = c(-1, 1)) +
  # re-order axes
  scale_y_discrete(limits = env.cor.labels[1:(length(env.cor.labels) - 1)],
                   labels = env.key) +
  scale_x_discrete(limits = env.cor.labels[length(env.cor.labels):2],
                   labels = env.key) +
  xlab(NULL) + 
  ylab(NULL) +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("figures/env_cor_plot.png",
       env.cor.plot,
       width = 6.5,
       height = 6,
       units = "in",
       dpi = 700)

## create plot of environmental factors over time
env.fac.long <- pivot_longer(env.fac.full, -Year) %>% 
  mutate(factor = ifelse(name %in% c("breeder",
                                     "juvenile",
                                     "helper",
                                     "pair.Fped"),
                         "relatedness",
                         name),
         class = ifelse(factor != "relatedness",
                        "none",
                        name) %>% 
           factor(levels = c("none", "juvenile", "helper", "breeder", "pair.Fped")))

env_fac_plot <- ggplot(env.fac.long, aes(x = Year, y = value)) +
  geom_line(data = filter(env.fac.long, factor != "relatedness"), size = 1) +
  geom_line(data = filter(env.fac.long, factor == "relatedness"),
            aes(col = class),
            size = 1) +
  geom_point(data = filter(env.fac.long, factor != "relatedness"), size = 2) +
  geom_point(data = filter(env.fac.long, factor == "relatedness"),
            aes(col = class),
            size = 2) +
  scale_color_manual(labels = c("Juveniles",
                                "Helpers",
                                "Breeders",
                                "Breeding\nPairs"),
                     values = c("#FF7F00",
                                "#33A02C",
                                "#CAB2D6",
                                "#1F78B4")) +
  facet_grid(factor~., scales = "free_y",
             labeller = labeller(factor = c(acorns = "Acorn Abund.",
                                            burn = "Burned Area",
                                            density = "Pop. Density",
                                            EQSOI = "EQSOI",
                                            relatedness = "Inbreeding"))) +
  theme_light() +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.title = element_blank())

ggsave("figures/env_fac.png",
       env_fac_plot,
       width = 6,
       height = 7,
       units = "in",
       dpi = 700)

## Contribution of vital rates to population growth
# load final LTRE results
load("data/generated_data/Demo_LTRE_results_20231106.rdata")

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
limit <- max(abs(contrib.means$contrib)*100) * c(-1.2, 1.2)

 
raw_contrib_plot <- ggplot(contrib.means, aes(x = vr, y = co_vr)) +
  geom_tile(aes(fill = contrib*100), col = "gray40") +
  scale_fill_distiller(palette = "PRGn", limit = limit) +
  scale_x_discrete(labels = gsub("c|q", "", levels(contrib.means$vr))) +
  scale_y_discrete(labels = gsub("c|q", "", levels(contrib.means$vr))) +
  labs(fill = expression("% Contribution to Var. in" ~ lambda)) +
  coord_cartesian(ylim = c(1, 28), xlim = c(1, 28), clip = "off") +
  facet_wrap(.~sex, labeller = labeller(sex = c(F = "Female-only Model",
                                                M = "Male-only Model")),
             nrow = 2, ncol = 1) +
  labs(x = NULL, y = "") +
  geom_text(aes(label = ifelse(abs(round(contrib*100, 1)) > 0,
                               round(contrib*100, 1), "")),
            size = 2.5) +
  # add text labels for past of present vital rates
  geom_text(x = 6, y = -2.1, label = "Current", col = "royalblue3") +
  geom_text(x = 18, y = -2.1, label = "Past", col = "tomato3") +
  geom_text(x = -1.35, y = 6, label = "Current", col = "royalblue3", angle = 90) +
  geom_text(x = -1.35, y = 18, label = "Past", col = "tomato3", angle = 90) +
  # add line segments to label vital rates
  geom_segment(aes(x = 0.5, xend = 12.5, y = -1.5, yend = -1.5),
               col = "royalblue3", linewidth = 1) +
  geom_segment(aes(x = 12.5, xend = 24.5, y = -1.5, yend = -1.5),
               col = "tomato3", linewidth = 1) +
  geom_segment(aes(x = -1, xend = -1, y = 0.5, yend = 12.5),
               col = "royalblue3", linewidth = 1) +
  geom_segment(aes(x = -1, xend = -1, y = 12.5, yend = 24.5),
               col = "tomato3", linewidth = 1) +
  theme_light() +
  theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 10),
        legend.title  = element_text(size = 16),
        legend.position = "bottom")

ggsave(file = "figures/raw_contrib_plot.png",
       raw_contrib_plot,
       dpi = 700,
       width = 6.5,
       height = 8,
       units = "in")

ggsave(file = "figures/raw_contrib_plot.pdf",
       raw_contrib_plot,
       dpi = 700,
       width = 6.5,
       height = 8)

## Tables of linear model results

# create key to adjust factor
envFactorKey <- c("Intercept", "Acorn Abundance", 
                  "Burned Area", "Time", 
                  "EQSOI",
                  "Population Density", "Inbreeding Coefficient", "Sex",
                  "Social Class", "Pair Relatedness",
                  "Pair Status", "Pair Relatedness:Pair Status",
                  "Year", "Territory ID", "Individual ID")
names(envFactorKey) <- c("(Intercept)", "acorns", "burn", "Year.scale",
                         "EQSOI", "density", "pedF",
                         "SexM", "social.classjuvenile",
                         "pair.Fped", "pair.statusold", "pair.Fped:pair.statusold",
                         "SD (Intercept Year)", "SD (Intercept Terr)", 
                         "SD (Intercept USFWS)")

# create table for survival models
modelsummary(list("Juvenile Survival" = Pj.model,
                                     "Helper Survival" = Ph.model,
                                     "Breeder Survival" = Pb.model),
                                stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
                                coef_rename = envFactorKey,
                                gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(16, 19)) %>%
  autofit() %>%
  save_as_docx(path = "figures/PModelTable.docx")

# create table for transition models
modelsummary(list("Probability to Pair" = B.model),
                                stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
                                coef_rename = envFactorKey,
                                gof_omit = "ICC",
                             output = "flextable") %>%
  hline(c(18, 21)) %>%
  autofit() %>%
  save_as_docx(path = "figures/BModelTable.docx")

# create table for fecundity models
modelsummary(list("Full Fecundity Model" = F.model.full,
                  "No Interaction Terms" = F.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_rename = envFactorKey,
             gof_omit = "ICC",
             shape = component + term + statistic ~ model,
             group_map = c("conditional" = "Conditional",
                           "zero_inflated" = "Zero-Inflated",
                           "dispersion" = "Dispersion"),
             output = "flextable") %>%
  hline(c(16, 18, 34, 36, 38, 39)) %>%
  autofit() %>%
  save_as_docx(path = "figures/FecModelTable.docx")

# create table comparing interaction term tests
modelsummary(list("Acorn Abundance" = F.lm.acorns, 
                  "Burned Area" = F.lm.burn, 
                  "Population Density" = F.lm.density, 
                  "EQSOI" = F.lm.EQSOI, 
                  "Pair Relatedness" = F.lm.Fped,
                  "Time" = F.lm.Year),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = c("acorns:pair.statusold" = "Acorns",
                          "burn:pair.statusold" = "Burned Area",
                          "density:pair.statusold" = "Population Density",
                          "EQSOI:pair.statusold" = "EQSOI",
                          "pair.Fped:pair.statusold" = "Pair Relatedness",
                          "Year.scale:pair.statusold" = "Time"),
             gof_omit = "ICC",
             shape = component + term + statistic ~ model,
             group_map = c("conditional" = "Conditional",
                           "zero_inflated" = "Zero-Inflated"),
             output = "flextable") %>%
  hline(c(12, 24)) %>%
  autofit() %>%
  save_as_docx(path = "figures/FecInteractModelTable.docx")

# create table for immigration models
modelsummary(list("Female Immigration" = IF.model,
                  "Male Immigration" = IM.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_rename = envFactorKey,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(12) %>%
  autofit() %>%
  save_as_docx(path = "figures/IModelTable.docx")

## Plot of population growth rate over time

# create data frame of observed lambda
lambda.df <- data.frame(lambda = lambda.obs,
                        Year = 1989:2020)

lambda.plot <- ggplot(lambda.df, aes(x = Year, y = lambda)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linewidth = 1) +
  labs(y = expression(lambda)) + 
  geom_hline(yintercept = pop.summary$geoMean[3], linewidth = 1, col = "chartreuse3") +
  theme_light() +
  theme(text = element_text(size = 12))

ggsave("figures/pop_growth.png",
       lambda.plot,
       width = 5,
       height = 5,
       units = "in",
       dpi = 700)

ggsave("figures/pop_growth.pdf",
       lambda.plot,
       width = 5,
       height = 5,
       units = "in")

## test for temporal trend in stage distribution, vital rates,
## and environmental factors
fac.MannKendall <- Obs.pop.vec[1:33, 2:5] %>%
  bind_cols(env.fac.full[, 2:9]) %>%
  bind_cols(vrObs.clean[, 2:13]) %>%
  apply(2, MannKendall)


MannKendallData <- data.frame(fac = c("Juvenile",
                                       "Helper",
                                       "New Pair",
                                       "Established Pair",
                                       "Acorn Abundance",
                                       "Burned Area",
                                       "EQSOI",
                                       "Population Density",
                                       "Breeding Pairs",
                                       "Breeders",
                                       "Helpers",
                                       "Juveniles",
                                       "Juvenile",
                                       "Helper",
                                       "Breeder",
                                       "Juvenile",
                                       "Helper",
                                       "Breeder",
                                       "Immigrant",
                                       "New Pair",
                                       "Established Pair",
                                       "New Pair",
                                       "Established Pair",
                                       "Immigration"),
                               tau = plyr::laply(1:24, function(x)
                                 fac.MannKendall[[x]]$tau[[1]]),
                               p = plyr::laply(1:24, function(x)
                                 fac.MannKendall[[x]]$sl[[1]]),
                               cat = c(rep("dist", 4),
                                       rep("env", 8),
                                       rep("vr", 12))) %>%
  {.[c(5, 6, 8, 7, 12, 11, 10, 9, 13:24, 1:4), 1:3]} %>%
  mutate(cat = c(rep("Environmental Factors", 4),
                 rep("Inbreeding", 4),
                 rep("Survival", 3),
                 rep("Probability to Pair", 4),
                 rep("Divorce", 2),
                 rep("Fecundity", 2),
                 "Immigration",
                 rep("Stage Distribution", 4))) %>%
  dplyr::select(Category = cat, Factor = fac,
                Tau = tau, "2-sided p-value" = p) %>%
  remove_rownames()

flextable(MannKendallData) %>%
  merge_v(j = "Category") %>%
  colformat_double(j = 3:4, digits = 3) %>%
  bold(j = 3:4, i = MannKendallData$`2-sided p-value` < 0.05, bold = TRUE) %>%
  theme_vanilla() %>%
  width(j = 1:2, 1.5) %>%
  width(j = 3, 0.75) %>%
  width(j = 4, 1.5) %>%
  save_as_docx(path = "figures/MannKendallTable.docx")

## Create table for VIF

VIFtable <-multicollinearity(Pj.model)[,c(1:2)] %>% data.frame(model = "Pj") %>%
  bind_rows(data.frame(multicollinearity(Ph.model)[,c(1:2)], model = "Ph")) %>%
  bind_rows(data.frame(multicollinearity(Pb.model)[,c(1:2)], model = "Pb")) %>%
  bind_rows(data.frame(multicollinearity(B.model)[,c(1:2)], model = "B")) %>%
  bind_rows(data.frame(multicollinearity(F.model.full)[,c(1:2)], 
                       model = c(rep("F1", 7), rep("F2", 6)))) %>%
  bind_rows(data.frame(multicollinearity(IM.model)[,c(1:2)], model = "IM")) %>%
  bind_rows(data.frame(multicollinearity(IF.model)[,c(1:2)], model = "IF"))

# filter out social class and pair status, combine pair.Fped with pedF
VIFtable <- VIFtable %>%
  filter(!Term %in% c("Sex", "social.class", "pair.status")) %>%
  mutate(Term = ifelse(Term == "Year.scale", "Year", Term),
         VIF = round(VIF, 2)) %>%
  mutate(Term = env.key[Term]) %>%
  mutate(Term = ifelse(is.na(Term) | grepl("Relatedness", Term),
                       "Inbreeding Coef./Relatedness", Term)) %>%
  # convert to wide format for printing table
  pivot_wider(id_cols = model, values_from = VIF, names_from = Term)

colnames(VIFtable)[1] <- "Model"

VIFtable$Model <-  c("Juvenile Survival",
                        "Helper Survival", 
                        "Breeder Survival",
                        "Probability to Pair",
                        "Fecundity (Conditional)",
                        "Fecundity (Zero-Inflated)",
                        "Male Immigration",
                        "Female Immigration")

flextable(VIFtable) %>%
  autofit() %>%
  save_as_docx(path = "figures/VIFscore.docx")
