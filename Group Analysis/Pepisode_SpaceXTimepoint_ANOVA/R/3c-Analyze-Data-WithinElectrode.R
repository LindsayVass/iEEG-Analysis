# Script Name:  3c-Analyze-Data-WithinElectrode.R
# Author:       Lindsay Vass
# Date:         29 July 2015
# Purpose:      This script will analyze the data output by "2-Clean-Data.R". It
#               will perform 2-way repeated measure ANOVAs with permutation
#               analyses to identify whether pepisode varies as a function of
#               Space (NS/FS) and Timepoint (Pre/Tele/Post). It collapses across 
#               Delta/Theta and performs all analyses within-electrode.

library(dplyr)
library(permute)

load('Rda/allCleanData.Rda')


# Functions ---------------------------------------------------------------

# return data from a specific frequency band, electrode, and trial time type
filterData <- function(episodeData, timeType, frequencyBand, electrode) {
  filteredData <- episodeData %>%
    filter(TrialTimeType == timeType,
           FrequencyBand == frequencyBand,
           ElectrodeID == electrode)
}

# sample equal number of trials from each spatial condition
sampleData <- function(thisData, numTrials) {
  thisData <- thisData %>%
    ungroup() %>%
    group_by(TrialSpaceType) %>%
    select(TrialNumber, TrialSpaceType) %>%
    unique() %>%
    sample_n(min(numTrials$Count)) %>%
    inner_join(thisData, by = c("TrialNumber", "TrialSpaceType"))
  
}

# clean up anova results
cleanAnova <- function(aovResult) {
  spaceF <- summary(aovResult)$"Error: TrialNumber"[[1]][["F value"]][1]
  timeF  <- summary(aovResult)$"Error: TrialNumber:TimeBin"[[1]][["F value"]][1]
  intxF  <- summary(aovResult)$"Error: TrialNumber:TimeBin"[[1]][["F value"]][2]
  output <- data.frame(Contrast = c("Space", "Timepoint", "Interaction"),
                       FVal = c(spaceF, timeF, intxF))
  return(output)
}

# shuffle data
shuffleData <- function(thisData) {
  plots  <- Plots(thisData$TrialNumber, type = "free")
  ctrl   <- how(plots = plots, within = Within(type = "free"))
  pOrder <- shuffle(nrow(thisData), ctrl)
  thisData$Pepisode <- thisData$Pepisode[pOrder]
  return(thisData)
}

# compare true F to distribution to get corrected P value
getCorrP <- function(trueAnova, permAnova, contrast) {
  trueData <- trueAnova %>%
    filter(Contrast == contrast)
  permData <- permAnova %>%
    filter(Contrast == contrast) %>%
    select(-Iteration) %>%
    rbind(trueData) %>%
    arrange(desc(FVal))
  corrP <- min(which(permData$FVal == trueData$FVal)) / nrow(permData)
  output <- trueData %>%
    mutate(CorrP = corrP)
  return(output)
}

# get p values associated with each # of significant electrodes
sortMaxElec <- function(maxElec) {
  rowPos <- data.frame(Row = seq(1, nrow(maxElec)))
  maxElec <- maxElec %>%
    select(-Iteration) %>%
    arrange(desc(NSigElectrodes))
  nVals <- data.frame(NSigElectrodes = unique(maxElec$NSigElectrodes))
  pVals <- data.frame()
  for (thisN in 1:nrow(nVals)) {
    nInd  <- min(which(maxElec == nVals$NSigElectrodes[thisN]))
    pVal  <- data.frame(P = nInd / nrow(maxElec))
    pVals <- rbind(pVals, pVal)
  }
  npVals <- cbind(nVals, pVals)
  output <- inner_join(maxElec, npVals, by = "NSigElectrodes")
}

# Get mean for each trial within frequency band ----------------------------

cleanData <- cleanData %>%
  group_by(ElectrodeID, TrialNumber, FrequencyBand, TrialSpaceType, TrialTimeType, TimeBin) %>%
  summarise(Pepisode = mean(Pepisode))
cleanData$TrialNumber <- factor(cleanData$TrialNumber)

# Run two-way ANOVAs ------------------------------------------------------

# analysis parameters
numPerm <- 1000 # for permutation testing
numCrossVal <- 100 # folds for true data ANOVA in case of uneven number of trials across space types
minTrials <- 5 # must have at least this many trials to run ANOVA

allTrueAnovaResults <- data.frame()
allPermAnovaResults <- data.frame()

for (thisFreqBand in 1:nlevels(cleanData$FrequencyBand)) {
  
  cat(paste0('\n\n', 
             levels(cleanData$FrequencyBand)[thisFreqBand],
             ' --------------------------------------------'))
  
  freqBand <- levels(cleanData$FrequencyBand)[thisFreqBand]
  
  for (thisTimeType in 1:nlevels(cleanData$TrialTimeType)) {
    
    cat(paste0('\n',
               levels(cleanData$TrialTimeType)[thisTimeType],
               '\n'))
    
    timeType <- levels(cleanData$TrialTimeType)[thisTimeType]
    
    pb <- txtProgressBar(min = 1, max = nlevels(cleanData$ElectrodeID), style = 3)
    for (thisElectrode in 1:nlevels(cleanData$ElectrodeID)) {
      setTxtProgressBar(pb, thisElectrode)
      
      electrode <- levels(cleanData$ElectrodeID)[thisElectrode]
      
      thisData <- filterData(cleanData, timeType, freqBand, electrode)
      
      # check that we have enough trials
      numTrials <- thisData %>%
        select(TrialNumber, TrialSpaceType) %>%
        unique() %>%
        ungroup() %>%
        group_by(TrialSpaceType) %>%
        summarise(Count = n())
      if (min(numTrials$Count) < minTrials) {
        thisResult <- data.frame(ElectrodeID = electrode,
                                 FrequencyBand = freqBand,
                                 TrialTimeType = timeType,
                                 SpaceF = NA,
                                 TimePointF = NA,
                                 InteractionF = NA,
                                 CorrP = NA)
        next
      }
      
      ### Run ANOVA on true data
      trueAnovaResults <- data.frame()
      # if equal number of NS and FS trials
      if (length(unique(numTrials$Count)) == 1) {
        aovResult  <- aov(Pepisode ~ TrialSpaceType * TimeBin + Error(TrialNumber / TimeBin), data = thisData) %>%
          cleanAnova() %>%
          mutate(ElectrodeID = electrode,
                 FrequencyBand = freqBand,
                 TrialTimeType = timeType)
        trueAnovaResults <- rbind(trueAnovaResults, aovResult)
      } else {
        for (i in 1:numCrossVal) {
          thisSample <- sampleData(thisData, numTrials)
          aovResult  <- aov(Pepisode ~ TrialSpaceType * TimeBin + Error(TrialNumber / TimeBin), data = thisSample) %>%
            cleanAnova() %>%
            mutate(ElectrodeID = electrode,
                   FrequencyBand = freqBand,
                   TrialTimeType = timeType)
          trueAnovaResults <- rbind(trueAnovaResults, aovResult)
        }
      }
      trueAnovaResults <- trueAnovaResults %>%
        group_by(Contrast) %>%
        summarise(FVal = mean(FVal)) %>%
        mutate(ElectrodeID = electrode,
               FrequencyBand = freqBand,
               TrialTimeType = timeType)
      
      ### Run ANOVA on shuffled data
      permAnovaResults <- data.frame()
      for (i in 1:numPerm) {
        setTxtProgressBar(pb, thisElectrode)
        thisPermData <- shuffleData(thisData) %>%
          sampleData(numTrials)
        permAov <- aov(Pepisode ~ TrialSpaceType * TimeBin + Error(TrialNumber / TimeBin), data = thisPermData) %>%
          cleanAnova() %>%
          mutate(ElectrodeID = electrode,
                 FrequencyBand = freqBand,
                 TrialTimeType = timeType,
                 Iteration = i)
        permAnovaResults <- rbind(permAnovaResults, permAov)
      }
      allPermAnovaResults <- rbind(allPermAnovaResults, permAnovaResults)
      
      ### Get corrected P values for true data
      trueCorrPResult <- data.frame()
      for (thisContrast in 1:nlevels(aovResult$Contrast)) {
        contrast <- levels(aovResult$Contrast)[thisContrast]
        thisResult <- getCorrP(trueAnovaResults, permAnovaResults, contrast)
        trueCorrPResult <- rbind(trueCorrPResult, thisResult)
      }
      allTrueAnovaResults <- rbind(allTrueAnovaResults, trueCorrPResult) 
    }
  }
}

# Calculate how many electrodes are significant by chance -----------------
allPermAnovaResults <- allPermAnovaResults %>%
  group_by(ElectrodeID, FrequencyBand, TrialTimeType, Contrast)
allPermAnovaResults$FrequencyBand <- factor(allPermAnovaResults$FrequencyBand)
allPermAnovaResults$TrialTimeType <- factor(allPermAnovaResults$TrialTimeType)
allPermAnovaResults$ElectrodeID   <- factor(allPermAnovaResults$ElectrodeID)
allPermAnovaResults$Contrast      <- factor(allPermAnovaResults$Contrast)
maxElec <- data.frame()
cat(paste0('\n\n', 
           'Determining number of electrodes significant by chance...\n',
           ' --------------------------------------------'))
pb <- txtProgressBar(min = 1, max = numPerm, style = 3)
for (i in 1:numPerm) {
  setTxtProgressBar(pb, i)
  thisSample <- allPermAnovaResults %>%
    sample_n(1)
  permCorrP <- data.frame()
  for (thisFreqBand in 1:nlevels(allPermAnovaResults$FrequencyBand)) {
    freqBand <- levels(allPermAnovaResults$FrequencyBand)[thisFreqBand]
    for (thisTimeType in 1:nlevels(allPermAnovaResults$TrialTimeType)) {
      timeType <- levels(allPermAnovaResults$TrialTimeType)[thisTimeType]
      for (thisElectrode in 1:nlevels(allPermAnovaResults$ElectrodeID)) {
        electrode <- levels(allPermAnovaResults$ElectrodeID)[thisElectrode]
        thisElecSample <- thisSample %>%
          filter(ElectrodeID == electrode,
                 FrequencyBand == freqBand,
                 TrialTimeType == timeType) %>%
          select(-Iteration) %>%
          ungroup()
        thisData  <- filterData(allPermAnovaResults, timeType, freqBand, electrode) %>%
          ungroup()
        for (thisContrast in 1:nlevels(allPermAnovaResults$Contrast)) {
          contrast <- levels(allPermAnovaResults$Contrast)[thisContrast]
          thisResult <- getCorrP(thisElecSample, thisData, contrast)
          permCorrP <- rbind(permCorrP, thisResult)
        }
      }
    }
  }
  permCorrP <- permCorrP %>%
    group_by(Contrast, FrequencyBand, TrialTimeType) %>%
    filter(CorrP < 0.05) %>%
    summarise(NSigElectrodes = n())
  thisMaxElec <- data.frame(Iteration = i,
                            NSigElectrodes = max(permCorrP$NSigElectrodes))
  maxElec <- rbind(maxElec, thisMaxElec)
}

### Calculate how many electrodes are significant and the associated significance
maxElecSig <- sortMaxElec(maxElec) %>%
  unique()
allSigElectrodes <- allTrueAnovaResults %>%
  filter(CorrP < 0.05) %>%
  group_by(Contrast, FrequencyBand, TrialTimeType) %>%
  summarise(NSigElectrodes = n()) %>%
  inner_join(maxElecSig)


# Get condition means for significant electrodes --------------------------

sigAnovas <- allSigElectrodes %>%
  filter(P < 0.05) %>%
  rename(NSigElectrodesP = P)

# Identify electrodes with significant effects
sigElectrodeList <- allTrueAnovaResults %>%
  filter(CorrP < 0.05) %>%
  inner_join(sigAnovas)

allConditionMeans <- data.frame()
for (thisAov in 1:nrow(sigElectrodeList)) {
  timeType  <- sigElectrodeList$TrialTimeType[thisAov]
  freqBand  <- sigElectrodeList$FrequencyBand[thisAov]
  electrode <- sigElectrodeList$ElectrodeID[thisAov]
  
  thisData <- filterData(cleanData, timeType, freqBand, electrode) %>%
    ungroup() %>%
    group_by(ElectrodeID, FrequencyBand, TrialSpaceType, TrialTimeType, TimeBin) %>%
    summarise(SE = sd(Pepisode) / sqrt(n()),
              Mean = mean(Pepisode))
  
  allConditionMeans <- rbind(allConditionMeans, thisData)  
}

save(file = 'Rda/allAnalyzedDataWithinElectrode.Rda', list = c("allSigElectrodes", "maxElecSig", "allTrueAnovaResults", "allPermAnovaResults", "sigElectrodeList", "sigAnovas", "allConditionMeans"))




