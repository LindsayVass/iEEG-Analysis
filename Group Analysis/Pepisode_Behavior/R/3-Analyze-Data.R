# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         19 June 2015
# Purpose:      This script will the clean data output by 2-Clean-Data.R 
#               Specifically, it will test whether pepisode during the
#               teleporter on trial N predicts the latency to enter the correct
#               arm on trial N+1.

library(dplyr)
library(broom)
library(permute)

load('Rda/allCleanData.Rda')


# Functions ---------------------------------------------------------------

cleanGlm <- function(inputData, electrode, frequencyBand) {
  glm(inputData$Latency ~ inputData$PrevTrialPepisode + inputData$TrialSpaceType + inputData$TrialTimeType,
      family = "quasipoisson") %>%
    tidy() %>%
    mutate(term = sub("inputData\\$", "", term),
           ElectrodeID = electrode,
           FrequencyBand = frequencyBand)
}

permAnalysis <- function(inputData, trueGlm, numPerm) {

  allPerms <- data.frame()
  pB <- txtProgressBar(min = 1, max = numPerm, initial = 1, style = 3)
  
  for (thisPerm in 1:numPerm) {
    
    setTxtProgressBar(pB, thisPerm)
    
    tempData <- inputData
    tempData$Latency <- tempData$Latency[shuffle(nrow(tempData))]
    
    permGlm <- cleanGlm(tempData, NA, NA)
    allPerms <- rbind(allPerms, permGlm)
    
  }
  
  allPerms <- rbind(allPerms, trueGlm)
  
  allPerms <- allPerms %>%
    select(term, statistic) %>%
    group_by(term)
  allPerms$term <- factor(allPerms$term)
  
  trueGlm$CorrectedP <- NA
  
  for (thisContrast in 1:nlevels(allPerms$term)) {
    
    # we expect PrevTrialPepisode to have negative relationship
    if (thisContrast == 2) {
      tempData <- allPerms %>%
        filter(term == levels(term)[thisContrast]) %>%
        arrange(statistic)
    } else {
      tempData <- allPerms %>%
        filter(term == levels(term)[thisContrast]) %>%
        arrange(desc(statistic))
    }
    
    trueStatistic <- trueGlm %>%
      filter(term == levels(allPerms$term)[thisContrast]) %>%
      select(statistic) %>%
      unlist()
    correctedP <- which(tempData$statistic == trueStatistic) / nrow(tempData)
    
    trueGlm$CorrectedP[thisContrast] = correctedP
    
  }
  
  return(trueGlm)
  
}

# Extract valid data ------------------------------------------------------

withinSessionData <- cleanData %>%
  group_by(TrialSpaceType, TrialTimeType, ElectrodeID, TrialNumber, FrequencyBand, TimeBin) %>%
  summarise(MeanPepisode = mean(Pepisode), Latency = mean(Latency))

# initialize output
validData <- data.frame()

for (thisElectrode in 1:nlevels(withinSessionData$ElectrodeID)) {
  for (thisFreqBand in 1:nlevels(withinSessionData$FrequencyBand)) {
      
    thisData <- withinSessionData %>%
      filter(ElectrodeID == levels(ElectrodeID)[thisElectrode] &
               FrequencyBand == levels(FrequencyBand)[thisFreqBand] &
               TimeBin == "Tele" &
               is.na(Latency) == FALSE) %>%
      ungroup() %>%
      arrange(TrialNumber)
    
    
    # keep only trials for which we have the previous trial's pepisode data
    thisTrialNumber <- thisData$TrialNumber
    prevTrialNumber <- thisTrialNumber - 1
    goodTrials      <- intersect(thisTrialNumber, prevTrialNumber) + 1
    
    thisValidData <- thisData %>%
      filter(TrialNumber %in% goodTrials) %>%
      select(-MeanPepisode) %>%
      mutate(PrevTrialPepisode = NA)
    
    # add in the previous trial's mean pepisode
    for (thisTrial in 1:nrow(thisValidData)) {
      
      thisValidData$PrevTrialPepisode[thisTrial] <- thisData %>%
        filter(TrialNumber == thisValidData$TrialNumber[thisTrial] - 1) %>%
        select(MeanPepisode) %>%
        unlist()
      
    }
    
    # append to output array
    validData <- rbind(validData, thisValidData)  
    
  }

}


# Perform GLM for each electrode ------------------------------------------

numPerm <- 1000
glmOutput <- data.frame()

for (thisElectrode in 1:nlevels(validData$ElectrodeID)) {
  
  cat(paste0('\n\n-----------------------------------------\n',
             levels(validData$ElectrodeID)[thisElectrode],
             '\n'))
  
  for (thisFreqBand in 1:nlevels(validData$FrequencyBand)) {
    
    cat(paste0(levels(validData$FrequencyBand)[thisFreqBand],
               '\n'))
    
    thisData <- validData %>%
      filter(ElectrodeID == levels(validData$ElectrodeID)[thisElectrode] &
               FrequencyBand == levels(validData$FrequencyBand)[thisFreqBand])
    
    # don't try to fit a model if fewer than 10 observations
    if (nrow(thisData) < 10)
      next
    
    trueGlm <- cleanGlm(thisData, levels(validData$ElectrodeID)[thisElectrode], levels(validData$FrequencyBand)[thisFreqBand])
    
    # run a permutation analysis
    correctedGlm <- permAnalysis(thisData, trueGlm, numPerm)
    glmOutput <- rbind(glmOutput, correctedGlm)
  }
}

glmOutput$ElectrodeID <- factor(glmOutput$ElectrodeID)

# Visualize the data for each electrode
p <- validData %>%
  ggplot(aes(x = PrevTrialPepisode, y = Latency)) +
  geom_point() +
  facet_wrap(~ElectrodeID)
dir.create('Figures')
ggsave('Figures/Latency_By_Teleporter_Pepisode_Each_Electrode.png', width = 16, height = 12)

# Save data
save(file = 'Rda/allAnalyzedData.Rda', list = c('validData', 'glmOutput', 'numPerm'))