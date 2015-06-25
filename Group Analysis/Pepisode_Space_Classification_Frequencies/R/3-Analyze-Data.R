# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         24 June 2015
# Purpose:      This script will analyze the data output by 2-Clean-Data.R. It
#               will test whether pepisode during the teleporter can be used to
#               classify the teleportation event based on its space type (NS/FS).
#               In contrast to the scripts in "Pepisode_Space_Classification",
#               this analysis will retain the frequency-specific information and
#               use all frequencies <= 8 Hz to create the pattern.

library(dplyr)
library(kernlab)
library(reshape2)

load('Rda/allCleanData.Rda')


# Functions ---------------------------------------------------------------

runClassifier <- function(trainingData, testingData) {
  
  tryCatch(fit <- ksvm(TrialSpaceType ~ ., data = trainingData),
           error = function(e) {
             fit <- ksvm(TrialSpaceType ~ ., data = trainingData)
           })
  predictions <- predict(fit, testingData[, 2:ncol(testingData)], type = "response")
  confusionMatrix <- table(predictions, testingData$TrialSpaceType)
  accuracy <- sum(diag(confusionMatrix)) / sum(confusionMatrix)
}

filterSample <- function(thisData, spaceType, trainSize, testSize) {
 
  trialList <- thisData %>%
    filter(TrialSpaceType == spaceType) %>%
    select(TrialID) %>%
    unique()
  trainList <- trialList %>%
    sample_n(trainSize)
  testList <- suppressMessages(anti_join(trialList, trainList)) %>%
    sample_n(testSize)
  return(list(train = trainList, test = testList))
  
}

getSampledData <- function(thisData, trialList) {
  
  output <- thisData %>%
    filter(TrialID %in% trialList$TrialID)
  
}

castWide <- function(inputData) {
  
  output <- inputData %>%
    dcast(TrialNumber + TrialSpaceType ~ Frequency, value.var = "Pepisode") %>%
    select(-(TrialNumber))
  
}

classifyOnSpace <- function(thisData, trainingSize, testingSize) {
  
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
  
  accuracy <- runClassifier(trainingData, testingData)
}


# Perform classification --------------------------------------------------

minFreq <- 1000 / (1830 / 3)

cleanData <- cleanData %>%
  filter(Frequency >= minFreq & Frequency <= 8) %>%
  mutate(TrialID = paste(ElectrodeID, TrialNumber, sep = "_"))
cleanData$Frequency <- factor(cleanData$Frequency)

classificationResults <- data.frame()

## skipping electrodes because not enough observations
badElectrodes <- c(21, 22, 23, 24, 29, 30, 31, 32, 35, 36, 37, 38)

numIterations <- 1000

for (thisElectrode in 1:nlevels(cleanData$ElectrodeID)) {
  
  if (thisElectrode %in% badElectrodes) {
    next()
  }
  
  cat(paste0('\n\n\n\n\n\n\n\n\n\n\n\n', levels(cleanData$ElectrodeID)[thisElectrode], '\n\n\n\n\n\n\n\n\n\n\n\n'))

  thisData <- cleanData %>%
    filter(TimeBin == "Tele" &
             ElectrodeID == levels(ElectrodeID)[thisElectrode])
  
  numObservations <- thisData %>%
    group_by(TrialSpaceType) %>%
    summarise(Count = n() / nlevels(Frequency))
  
  numTimeObservations <- thisData %>%
    group_by(TrialTimeType, TrialSpaceType) %>%
    summarise(Count = n() / nlevels(Frequency))
  
  # percentage of the data set to use to train the classifier
  trainingPercent <- 0.75 
  trainingSize    <- floor(trainingPercent * min(numObservations$Count))
  testingSize     <- min(numObservations$Count) - trainingSize
  trainingTimeSize <- floor(trainingPercent * min(numTimeObservations$Count))
  testingTimeSize  <- min(numTimeObservations$Count) - trainingTimeSize
  
  for (thisIteration in 1:numIterations) {
    
    #### Run classification using both NT and FT ####
    accuracy <- classifyOnSpace(thisData, trainingSize, testingSize)
    thisResult <- data.frame(ElectrodeID = levels(cleanData$ElectrodeID)[thisElectrode],
                             Model = "Both",
                             Iteration = thisIteration,
                             Accuracy = accuracy)
    classificationResults <- rbind(classificationResults, thisResult)
    
    #### Run classification using NT only ####
    inputData <- thisData %>%
      filter(TrialTimeType == "NT")
    
    accuracy <- classifyOnSpace(inputData, trainingTimeSize, testingTimeSize)
    thisResult <- data.frame(ElectrodeID = levels(cleanData$ElectrodeID)[thisElectrode],
                             Model = "NT",
                             Iteration = thisIteration,
                             Accuracy = accuracy)
    classificationResults <- rbind(classificationResults, thisResult)
    
    #### Run classification using FT only ####
    inputData <- thisData %>%
      filter(TrialTimeType == "FT")
    
    accuracy <- classifyOnSpace(inputData, trainingTimeSize, testingTimeSize)
    thisResult <- data.frame(ElectrodeID = levels(cleanData$ElectrodeID)[thisElectrode],
                             Model = "FT",
                             Iteration = thisIteration,
                             Accuracy = accuracy)
    classificationResults <- rbind(classificationResults, thisResult)
    
    
  }
}

meanClassification <- classificationResults %>%
  group_by(ElectrodeID, Model) %>%
  summarise(Accuracy = mean(Accuracy))

ttestResults <- data.frame()

for (thisModel in 1:nlevels(meanClassification$Model)) {
  
  thisData <- meanClassification %>%
    filter(Model == levels(Model)[thisModel])
  thisTest <- t.test(thisData$Accuracy, alternative = "greater", mu = 0.5)
  thisResult <- data.frame(Model = levels(thisData$Model)[thisModel],
                           TStatistic = thisTest$statistic,
                           PValueUncorr = thisTest$p.value,
                           PValueBonf = thisTest$p.value * (nlevels(meanClassification$Model)))
  ttestResults <- rbind(ttestResults, thisResult)
  
}


save(file = 'Rda/allClassificationResults.Rda', list = c('classificationResults', 'ttestResults', 'meanClassification'))