# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         22 June 2015
# Purpose:      This script will analyze the data output by 2-Clean-Data.R. It
#               will test whether pepisode during the teleporter can be used to
#               classify the teleportation event based on its space type (NS/FS).

library(dplyr)
library(kernlab)

load('Rda/allCleanData.Rda')


# Functions ---------------------------------------------------------------

runClassifier <- function(trainingData, testingData) {
  
  tryCatch(fit <- ksvm(TrialSpaceType ~ ., data = trainingData),
           error = function(e) {
             fit <- ksvm(TrialSpaceType ~ ., data = trainingData)
           })
  
  predictions <- predict(fit, testingData[, 2], type = "response")
  confusionMatrix <- table(predictions, testingData$TrialSpaceType)
  accuracy <- sum(diag(confusionMatrix)) / sum(confusionMatrix)
  
}


# Perform classification --------------------------------------------------

classificationResults <- data.frame()

for (thisElectrode in 1:nlevels(cleanData$ElectrodeID)) {
  
  cat('\n\n-------------------------------------------------\n\n')
  cat(levels(cleanData$ElectrodeID)[thisElectrode])
  
  for (thisFreqBand in 1:nlevels(cleanData$FrequencyBand)) {
    
    cat('\n')
    cat(levels(cleanData$FrequencyBand)[thisFreqBand])
    cat('\n')
    
    thisData <- cleanData %>%
      filter(TimeBin == "Tele" &
               ElectrodeID == levels(ElectrodeID)[thisElectrode] &
               FrequencyBand == levels(FrequencyBand)[thisFreqBand])
    
    numObservations <- thisData %>%
      group_by(TrialSpaceType) %>%
      summarise(Count = n())
    
    numTimeObservations <- thisData %>%
      group_by(TrialTimeType, TrialSpaceType) %>%
      summarise(Count = n())
    
    # percentage of the data set to use to train the classifier
    trainingPercent <- 0.75 
    trainingSize    <- floor(trainingPercent * min(numObservations$Count))
    testingSize     <- min(numObservations$Count) - trainingSize
    trainingTimeSize <- floor(trainingPercent * min(numTimeObservations$Count))
    testingTimeSize  <- min(numTimeObservations$Count) - trainingTimeSize
    
    numIterations <- 1000
    
    
    for (thisIteration in 1:numIterations) {
      
      #### Run classification using both NT and FT ####
      
      # first sample so that both conditions have the same number of observations
      equalSizeData <- thisData %>%
        group_by(TrialSpaceType) %>%
        sample_n(min(numObservations$Count))
      
      # make training and testing sets
      trainingData <- equalSizeData %>%
        sample_n(trainingSize) 
      testingData  <- suppressMessages(anti_join(equalSizeData, trainingData))
      
      # remove unneeded columns
      trainingData <- trainingData %>%
        select(TrialSpaceType, Pepisode)
      testingData  <- testingData %>%
        select(TrialSpaceType, Pepisode)
      
      # run classifier
      accuracy <- runClassifier(trainingData, testingData)
      
      thisResult <- data.frame(ElectrodeID = levels(cleanData$ElectrodeID)[thisElectrode],
                               FrequencyBand = levels(cleanData$FrequencyBand)[thisFreqBand],
                               Model = "Both",
                               Iteration = thisIteration,
                               Accuracy = accuracy)
      classificationResults <- rbind(classificationResults, thisResult)
      
      #### Run classification using NT only ####
      
      # first sample so that both conditions have the same number of observations
      equalSizeData <- thisData %>%
        filter(TrialTimeType == "NT") %>%
        group_by(TrialSpaceType) %>%
        sample_n(min(numTimeObservations$Count))
      
      # make training and testing sets
      trainingData <- equalSizeData %>%
        sample_n(trainingTimeSize) 
      testingData  <- suppressMessages(anti_join(equalSizeData, trainingData))
      
      # remove unneeded columns
      trainingData <- trainingData %>%
        select(TrialSpaceType, Pepisode)
      testingData  <- testingData %>%
        select(TrialSpaceType, Pepisode)
      
      # run classifier
      accuracy <- runClassifier(trainingData, testingData)
      
      thisResult <- data.frame(ElectrodeID = levels(cleanData$ElectrodeID)[thisElectrode],
                               FrequencyBand = levels(cleanData$FrequencyBand)[thisFreqBand],
                               Model = "NT",
                               Iteration = thisIteration,
                               Accuracy = accuracy)
      classificationResults <- rbind(classificationResults, thisResult)
      
      #### Run classification using FT only ####
      
      # first sample so that both conditions have the same number of observations
      equalSizeData <- thisData %>%
        filter(TrialTimeType == "FT") %>%
        group_by(TrialSpaceType) %>%
        sample_n(min(numTimeObservations$Count))
      
      # make training and testing sets
      trainingData <- equalSizeData %>%
        sample_n(trainingTimeSize) 
      testingData  <- suppressMessages(anti_join(equalSizeData, trainingData))
      
      # remove unneeded columns
      trainingData <- trainingData %>%
        select(TrialSpaceType, Pepisode)
      testingData  <- testingData %>%
        select(TrialSpaceType, Pepisode)
      
      # run classifier
      accuracy <- runClassifier(trainingData, testingData)
      
      thisResult <- data.frame(ElectrodeID = levels(cleanData$ElectrodeID)[thisElectrode],
                               FrequencyBand = levels(cleanData$FrequencyBand)[thisFreqBand],
                               Model = "FT",
                               Iteration = thisIteration,
                               Accuracy = accuracy)
      classificationResults <- rbind(classificationResults, thisResult)
      
    }
  }
}