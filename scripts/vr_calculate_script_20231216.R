# This script calculates vital rates based on the new 4 stage model
# Jeremy Summers
# May 2023

library(tidyverse)
setwd("C:/Users/jtsum/OneDrive/2019-0-spring/elasticity_project/rotation3/final_fileset/")

# Inputs are the list of individuals filtered for the study tract,
# the sex the rates are calculated for, 
# the start and end year of the data, 
# the number of iterations, and the random seed for assigning sex to
# NA individuals
vr_estimate <- function(LOI, sex, startYr, endYr, iter, seed){
  #set seed
  set.seed(seed)
  #create array to store results
  #rows are the number of iterations, columns are the number of vital rates, 
  # and the third dimension is
  #the number of years
  results <- array(dim = c(iter, 12, endYr - startYr + 1),
                   dimnames = list(c(1:iter),
                                   c("Pj", "Ph", "Pb", "Bj", "Bh", "Bb", "Bi",
                                     "Dn", "Do", "I", "Fn", "Fo"), 
                                   c(startYr:endYr)))
  #also save the population vectors generated for each re-sampling
  pop.vec <- array(dim = c(iter, 4, endYr - startYr + 2),
                   dimnames = list(c(1:iter),
                                   c("j", "h", "n", "o")), 
                   c(startYr:(endYr + 1)))
  
  ## Sex assignment
  # create matrix of all unsexed individuals with randomly assigned "F" or "M" for sex
  unsexed <- filter(LOI, is.na(Sex))$USFWS %>% unique() %>%
    # generate random sexes equal to number of individuals per iteration
    {matrix(sample(c("F", "M"), length(.)*iter, replace = TRUE, prob = c(0.5, 0.5)),
            nrow = length(.), ncol = iter,
            dimnames = list(., c(1:iter)))}
  
  # run script for each iteration of sex assignments
  for(i in 1:iter){
    # assign unsexed individuals
    LOIassigned <- LOI
    LOIassigned[is.na(LOIassigned$Sex), 
                "Sex"] <- unsexed[paste(LOIassigned[is.na(LOIassigned$Sex), "USFWS"])
                                  , i]
    
    # filter for sex in analysis
    LOIassigned <- LOIassigned[LOIassigned$Sex == sex,]
    
    ### create population vector
    pop.vec[i,,] <- table(interaction(LOIassigned$social.class, LOIassigned$pair.status),
                          LOIassigned$Year)[c("juvenile.none", "helper.none", 
                                              "breeder.new", "breeder.old"),
                                            paste(startYr:(endYr + 1))]
    
    ### survival rates
    # find the number of individual who survive per social class per year
    results[i, c("Pb", "Ph", "Pj"),] <- (table(LOIassigned$social.class, 
                                              LOIassigned$survival, 
                                              LOIassigned$Year)[,2,]/
      # divide by the total number of individuals in that social class per year
      apply(table(LOIassigned$social.class, 
                  LOIassigned$survival, 
                  LOIassigned$Year), c(1, 3), sum)) %>%
      {.[,paste(startYr:endYr)]}
    
    ### transition rates
    # first calculate for helper and juveniles
    results[i, c("Bh", "Bj"),] <- (table(LOIassigned$social.class, 
                                               LOIassigned$transition, 
                                               LOIassigned$Year)[,2,]/
      # divide by the total number of individuals in that social class per year
      apply(table(LOIassigned$social.class, 
                  LOIassigned$survival, 
                  LOIassigned$Year), c(1, 3), sum)) %>%
      {.[c("helper", "juvenile"), paste(startYr:endYr)]}
    
    # calculate pair prob for breeders
    results[i, "Bb",] <- (table(LOIassigned$back_transition, LOIassigned$Year)[2,]/
      colSums(table(LOIassigned$back_transition, LOIassigned$Year))) %>%
      {.[paste(startYr:endYr)]}
    
    ### divorce rates
    results[i, c("Dn", "Do"),] <- (table(LOIassigned$pair.status,
                                         LOIassigned$divorce,
                                         LOIassigned$Year)[,2,]/
      apply(table(LOIassigned$pair.status,
                  LOIassigned$divorce,
                  LOIassigned$Year), c(1, 3), sum)) %>%
      {.[c("new", "old"), paste(startYr:endYr)]}
    
    ### immigration rate
    results[i, "I",] <- table(LOIassigned$category == "immigrant",
                              # immigrants arriving in one year are the result of
                              # events in the previous year
                              LOIassigned$Year - 1)[2, paste(startYr:endYr)]/
      # divide by the total number of individuals in the population the previous year
      colSums(pop.vec[i, ,1:dim(results)[3]])
    
    ### immigrant pair rate
    results[i, "Bi",] <- (table(LOIassigned$category,
                               LOIassigned$pair.status != "none",
                               LOIassigned$Year - 1)[1, 2,]/
      colSums(table(LOIassigned$category,
                    LOIassigned$pair.status != "none",
                    LOIassigned$Year - 1)[1,,])) %>%
      {.[paste(startYr:endYr)]}
    
    ### fecundity
    results[i, c("Fn", "Fo"),] <- filter(LOIassigned, pair.status != "none") %>%
      group_by(Year = (Year - 1), pair.status) %>%
      dplyr::summarize(Fec = mean(fecundity)) %>% pivot_wider(id_cols = pair.status,
                                                       names_from = Year,
                                                       values_from = Fec) %>%
      {as.matrix(.[,-1])} %>% {.[, paste(startYr:endYr)]}
  }
  return(list(results, pop.vec))
}

# calculate vital rates

# filter the LOI for Terrs within the study site
FullLOI <- read_tsv("data/FullLOI.txt") %>% data.frame()
load("data/generated_data/densityCalcDemo.rdata")
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

# there are 259 records for individuals who eventually came back to the study tract
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

# now try running my vr script on this new data frame
vrSampleResults <- vr_estimate(data.frame(StudySiteLOI.clean), "F", 
                               1988, 2020, 100, 5657)
vr.mat <- vrSampleResults[[1]]
pop.vec <- vrSampleResults[[2]]
save(vr.mat,pop.vec,vrSampleResults, 
     file = "data/generated_data/vr_clean_F_4stageDemo.rdata")

# repeat for males
vrSampleResults <- vr_estimate(data.frame(StudySiteLOI.clean), "M", 
                               1988, 2020, 100, 5657)
vr.mat <- vrSampleResults[[1]]
pop.vec <- vrSampleResults[[2]]
save(vr.mat,pop.vec,vrSampleResults, 
     file = "data/generated_data/vr_clean_M_4stageDemo.rdata")
