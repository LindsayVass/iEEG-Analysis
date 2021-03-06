# Script Name:  4c-Plot-Best-Electrode-Data.R
# Author:       Lindsay Vass
# Date:         28 August 2015
# Purpose:      This script will focus on producing more detailed plots of the 
#               best electrode's results. This includes a scatterplot of the
#               pepisode values at each frequency for short/long distance, and a
#               way of visualizing the classification for one specific iteration.

library(dplyr)
library(ggplot2)
library(ggthemes)
library(data.table)
library(grid)
library(R.matlab)

load('Rda/allClassificationResults.Rda')
load('Rda/allCleanData.Rda')

#dir.create('Figures/Single Electrode/')

# Functions ---------------------------------------------------------------

filterSample <- function(thisData, spaceType, trainSize, testSize) {
  trialList <- thisData %>%
    filter(TrialSpaceType == spaceType) %>%
    select(TrialNumber) %>%
    unique()
  trainList <- trialList %>%
    sample_n(trainSize)
  testList <- suppressMessages(anti_join(trialList, trainList)) %>%
    ungroup() %>%
    sample_n(testSize)
  return(list(train = trainList, test = testList))
  
}

getSampledData <- function(thisData, trialList) {
  
  output <- thisData %>%
    filter(TrialNumber %in% trialList$TrialNumber)
  
}

castWide <- function(inputData) {
  
  output <- inputData %>%
    dcast(TrialNumber + TrialSpaceType ~ Frequency, value.var = "Pepisode") %>%
    select(-(TrialNumber))
  
}

sampleData <- function(thisData, trainingSize, testingSize) {
  
  # Extract training and testing lists of trials
  nsLists <- filterSample(thisData, "NS", trainingSize, testingSize)
  fsLists <- filterSample(thisData, "FS", trainingSize, testingSize)
  
  nsTrain <- getSampledData(thisData, nsLists$train)
  nsTest  <- getSampledData(thisData, nsLists$test)
  fsTrain <- getSampledData(thisData, fsLists$train)
  fsTest  <- getSampledData(thisData, fsLists$test)
  
  # cast to wide
  nsTrainWide <- castWide(nsTrain)
  nsTestWide  <- castWide(nsTest)
  fsTrainWide <- castWide(fsTrain)
  fsTestWide  <- castWide(fsTest)
  
  # concatenate
  trainingData <- rbind(nsTrainWide, fsTrainWide)
  testingData  <- rbind(nsTestWide, fsTestWide)
  
  # get testing trial numbers
  nsTestTrials <- nsLists$test %>%
    mutate(TrialSpaceType = "NS")
  fsTestTrials <- fsLists$test %>%
    mutate(TrialSpaceType = "FS")
  allTestTrials <- rbind(nsTestTrials, fsTestTrials)
  
  return(list(train = trainingData, test = testingData, testTrials = allTestTrials))
}

# Prepare data for plotting -----------------------------------------------

# select the top classification results
numResults <- 2

topClass <- allClassificationResults %>%
  ungroup() %>%
  arrange(desc(Accuracy))
topClass <- topClass[1:numResults,] %>%
  select(ElectrodeID, Model) %>%
  rename(TrialTimeType = Model)

# get raw data for this classification
minFreq <- 1000 / (1830 / 3)
maxFreq <- 8
classData <- cleanData %>%
  inner_join(topClass, by = c("ElectrodeID", "TrialTimeType")) %>%
  filter(TimeBin == 'Tele',
         Frequency > minFreq,
         Frequency <= maxFreq) %>%
  group_by(ElectrodeID)
classData$ElectrodeID <- factor(classData$ElectrodeID)

# Make scatterplot of pepisode at each frequency --------------------------
topClass <- topClass %>%
  mutate(Frequency = NA,
         MeanNS = NA,
         MeanFS = NA)
for (thisResult in 1:nrow(topClass)) {
  meanFreqPepisode <- classData %>%
    filter(ElectrodeID == as.character(topClass$ElectrodeID[thisResult])) %>%
    group_by(TrialSpaceType, Frequency) %>%
    summarise(MeanPepisode = mean(Pepisode),
              SEM = sd(Pepisode) / sqrt(n()))
  scatterP <- meanFreqPepisode %>%
    ggplot(aes(x = Frequency,
               y = MeanPepisode,
               ymin = MeanPepisode - SEM,
               ymax = MeanPepisode + SEM,
               color = factor(TrialSpaceType, labels = c('Short', 'Long')))) +
    geom_point(size = 5) +
    geom_pointrange() +
    scale_x_log10(breaks = unique(meanFreqPepisode$Frequency), labels = round(unique(meanFreqPepisode$Frequency), digits = 2)) + 
    theme_few() +
    theme(text = element_text(size = 18),
          axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = -0.25),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.position = c(0.15, 0.75)) +
    labs(x = 'Frequency (Hz)',
         y = expression('Mean P'['Episode']),
         colour = 'Distance Teleported')
  scatterP
  ggsave(paste0('Figures/Single Electrode/Scatter_Pepisode_Each_Distance_Frequency_', as.character(topClass$ElectrodeID[thisResult]), '_', as.character(topClass$TrialTimeType[thisResult]), '.pdf'), useDingbats = FALSE)
  
  # find frequency showing biggest difference
  bigDiff <- meanFreqPepisode %>%
    dcast(Frequency ~ TrialSpaceType, value.var = 'MeanPepisode') %>%
    mutate(PepisodeDifference = abs(NS - FS)) %>%
    arrange(desc(PepisodeDifference))
  topClass$Frequency[thisResult] <- bigDiff$Frequency[1]
  topClass$MeanNS[thisResult]    <- bigDiff$NS[1]
  topClass$MeanFS[thisResult]    <- bigDiff$FS[1]
}

# Make plot of classification of individual trials ------------------------
numPerm          <- 1000
trainingPercent  <- 0.75 
for (thisResult in 1:nlevels(classData$ElectrodeID)) {
  
  # set classification parameters
  numObservations <- classData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisResult]) %>%
    group_by(ElectrodeID, TrialSpaceType) %>%
    summarise(Count = n() / length(unique(Frequency)))
  trainingSize     <- floor(trainingPercent * min(numObservations$Count))
  testingSize      <- min(numObservations$Count) - trainingSize

  
  oneElecResults <- vector(mode = "list", length = numPerm)
  thisData <- classData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisResult])
  for (i in 1:numPerm) {
    thisSample <- sampleData(thisData, trainingSize, testingSize)
    fit <- ksvm(TrialSpaceType ~ ., data = thisSample$train)
    predictions <- predict(fit, thisSample$test, type = "response") %>%
      cbind(thisSample$testTrials)
    names(predictions) <- c('Classification', 'ElectrodeID', 'TrialNumber', 'TrialSpaceType')
    oneElecResults[[i]] <- predictions
  }
  oneElecResults <- rbindlist(oneElecResults)
  
  meanTrialResults <- oneElecResults %>% 
    select(-ElectrodeID) %>%
    group_by(TrialNumber, TrialSpaceType, Classification) %>%
    summarise(Count = n()) %>%
    as.data.frame() %>%
    arrange(TrialNumber) %>%
    ungroup() %>%
    group_by(TrialNumber) %>%
    mutate(TotalCount = sum(Count)) %>%
    filter(Classification == "NS") %>%
    mutate(PropNS = Count / TotalCount - 0.5)
  meanTrialResults$TrialSpaceType <- factor(meanTrialResults$TrialSpaceType, levels = c('NS', 'FS'), labels = c('Short', 'Long'))
  
  trialAccuracy <- meanTrialResults %>%
    select(TrialNumber, TrialSpaceType, PropNS) %>%
    unique() %>%
    mutate(FinalClass = ifelse(PropNS > 0, "Short", "Long"),
           Accurate = ifelse(TrialSpaceType == FinalClass, TRUE, FALSE)) %>%
    group_by(Accurate) %>%
    summarise(Count = n())
  accuracyLabel <- paste0(trialAccuracy$Count[which(trialAccuracy$Accurate == TRUE)], '/', sum(trialAccuracy$Count), ' trials correctly classified')
  
  trialClassP <- meanTrialResults %>%
    transform(TrialNumber = reorder(TrialNumber, PropNS)) %>%
    ggplot(aes(x = factor(TrialNumber),
               y = PropNS,
               fill = TrialSpaceType)) +
    geom_bar(stat = 'identity', position = 'identity') +
    scale_y_continuous(limits = c(-0.5, 0.5),
                       breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                       labels = c(0, 25, 50, 75, 100)) +
    annotation_custom(grob = textGrob(accuracyLabel, gp = gpar(fontsize = 24)), xmin = -15, xmax = Inf, ymin = 0, ymax = Inf) +
    coord_flip() +
    theme_few() +
    theme(text = element_text(size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = c(0.15, 0.88),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(vjust = -0.5)) +
    labs(y = '% Iterations Classified as "Short"',
         x = 'Trial',
         fill = 'Distance Teleported') 
  trialClassP
  ggsave(paste0('Figures/Single Electrode/Bar_Trial_Classification_', as.character(levels(classData$ElectrodeID)[thisResult]), '_', as.character(topClass$Model[which(topClass$ElectrodeID == levels(classData$ElectrodeID)[thisResult])]), '.pdf'))

}


# Find best trials for raw data plots -------------------------------------
bestData <- cleanData %>%
  inner_join(topClass) %>%
  filter(TimeBin == "Tele") %>%
  mutate(MeanDiff = ifelse(TrialSpaceType == "NS", abs(MeanNS - Pepisode), abs(MeanFS - Pepisode))) %>%
  group_by(ElectrodeID, TrialSpaceType) %>%
  arrange(MeanDiff) %>%
  top_n(1, -MeanDiff) %>%
  group_by(ElectrodeID, TrialSpaceType) %>%
  sample_n(1) %>% # if multiple trials with same value
  select(ElectrodeID, TrialNumber, TrialSpaceType, TrialTimeType, Frequency, Pepisode)

#dir.create('mat')
writeMat('mat/bestPepisodeTrials.mat', bestData = bestData)
save(file = 'Rda/bestPepisodeTrials.Rda', list = 'bestData')


