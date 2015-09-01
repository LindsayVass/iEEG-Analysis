# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         24 June 2015
# Purpose:      This script will analyze the data output by 2-Clean-Data.R. It
#               will test whether pepisode during the teleporter can be used to
#               classify the teleportation event based on its space type (NS/FS).
#               In contrast to the scripts in "Pepisode_Space_Classification",
#               this analysis will retain the frequency-specific information and
#               use all frequencies <= 8 Hz to create the pattern. In contrast to
#               3-Analyze-Data.R, it will perform all statistics within-electrode
#               and determine whether there are more significant electrodes than
#               would be expected by chance.

library(dplyr)
library(kernlab)
library(reshape2)
library(permute)

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

runClassification <- function(inputData, trainingSize, testingSize, electrode, model, iteration) {
  
  accuracy <- classifyOnSpace(inputData, trainingSize, testingSize)
  thisResult <- data.frame(ElectrodeID = electrode,
                           Model = model,
                           Iteration = iteration,
                           Accuracy = accuracy)
}

shuffleData <- function(inputData) {
  permOrder <- shuffle(nrow(inputData))
  inputData$Pepisode <- inputData$Pepisode[permOrder]
  return(inputData)
}

getCorrectedP <- function(permData, allBootstrapClassification, electrode, model) {
  bootstrapData <- allBootstrapClassification %>%
    filter(ElectrodeID == electrode,
           Model == model)
  permData <- permData %>%
    filter(ElectrodeID == electrode,
           Model == model)
  trueP <- min(which(bootstrapData$TVal == permData$TVal)) / nrow(bootstrapData)
  thisResult <- data.frame(ElectrodeID = electrode,
                           Model = model,
                           TVal = permData$TVal,
                           CorrP = trueP)
}

sortPerms <- function(permData) {
  rowPos <- data.frame(Row = seq(1, nrow(permData)))
  permData <- permData %>%
    arrange(desc(NSigElectrodes))
  nVals <- data.frame(NSigElectrodes = unique(permData$NSigElectrodes))
  pVals <- data.frame()
  for (thisN in 1:nrow(nVals)) {
    nInd  <- min(which(permData == nVals$NSigElectrodes[thisN]))
    pVal  <- data.frame(P = nInd / nrow(permData))
    pVals <- rbind(pVals, pVal)
  }
  npVals <- cbind(nVals, pVals)
  output <- inner_join(permData, npVals, by = "NSigElectrodes")
}

# Perform classification --------------------------------------------------

minFreq <- 1000 / (1830 / 3)

cleanData <- cleanData %>%
  filter(Frequency >= minFreq & Frequency <= 8) %>%
  mutate(TrialID = paste(ElectrodeID, TrialNumber, sep = "_"))
cleanData$Frequency <- factor(cleanData$Frequency)


numIterations <- 1000

allClassificationResults <- data.frame()
allMeanClassificationResults <- data.frame()
allBootstrapClassification <- data.frame()
for (thisElectrode in 1:nlevels(cleanData$ElectrodeID)) {
  
  classificationResults <- data.frame()
  permutationResults    <- data.frame()
  
  cat(paste0('\n\n\n\n\n\n\n\n\n\n\n\n', levels(cleanData$ElectrodeID)[thisElectrode], '\n\n\n\n\n\n\n\n\n\n\n\n'))
  
  electrode = levels(cleanData$ElectrodeID)[thisElectrode]

  thisData <- cleanData %>%
    filter(TimeBin == "Tele" &
             ElectrodeID == electrode)
  
  numObservations <- thisData %>%
    group_by(TrialSpaceType) %>%
    summarise(Count = n() / nlevels(Frequency))
  
  numTimeObservations <- thisData %>%
    group_by(TrialTimeType, TrialSpaceType) %>%
    summarise(Count = n() / nlevels(Frequency))
  
  # percentage of the data set to use to train the classifier
  trainingPercent  <- 0.75 
  trainingSize     <- floor(trainingPercent * min(numObservations$Count))
  testingSize      <- min(numObservations$Count) - trainingSize
  trainingTimeSize <- floor(trainingPercent * min(numTimeObservations$Count))
  testingTimeSize  <- min(numTimeObservations$Count) - trainingTimeSize
  
  if (testingSize < 3 | testingTimeSize < 3) {
    next
  }
  
  for (thisIteration in 1:numIterations) {
    
    #### Run classification using both NT and FT ####
    thisResult <- runClassification(thisData, trainingSize, testingSize, electrode, "Both", thisIteration)
    classificationResults <- rbind(classificationResults, thisResult)
    
    permData <- shuffleData(thisData)
    thisPermResult <- runClassification(permData, trainingSize, testingSize, electrode, "Both", thisIteration)
    permutationResults <- rbind(permutationResults, thisPermResult)
    
    #### Run classification using NT only ####
    inputData <- thisData %>%
      filter(TrialTimeType == "NT")
    thisResult <- runClassification(inputData, trainingTimeSize, testingTimeSize, electrode, "NT", thisIteration)
    classificationResults <- rbind(classificationResults, thisResult)
    
    permData <- shuffleData(inputData)
    thisPermResult <- runClassification(permData, trainingTimeSize, testingTimeSize, electrode, "NT", thisIteration)
    permutationResults <- rbind(permutationResults, thisPermResult)
    
    #### Run classification using FT only ####
    inputData <- thisData %>%
      filter(TrialTimeType == "FT")
    thisResult <- runClassification(inputData, trainingTimeSize, testingTimeSize, electrode, "FT", thisIteration)
    classificationResults <- rbind(classificationResults, thisResult)
    
    permData <- shuffleData(inputData)
    thisPermResult <- runClassification(permData, trainingTimeSize, testingTimeSize, electrode, "FT", thisIteration)
    permutationResults <- rbind(permutationResults, thisPermResult)
    
  }
  
  allClassificationResults <- rbind(allClassificationResults, classificationResults)
  
  meanClassification <- classificationResults %>%
    group_by(Model) %>%
    summarise(TVal = t.test(Accuracy, mu = 0.5, alternative = "greater")$statistic,
              UncorrP = t.test(Accuracy, mu = 0.5, alternative = "greater")$p.value,
              SEM = sd(Accuracy) / sqrt(n()),
              Accuracy = mean(Accuracy))
  
  # bootstrap permutation values
  bootstrapClassification <- data.frame()
  for (thisIteration in 1:numIterations) {
    permData <- permutationResults %>%
      group_by(Model) %>%
      sample_n(numIterations, replace = TRUE) %>%
      summarise(TVal = t.test(Accuracy, mu = 0.5, alternative = "greater")$statistic)
    bootstrapClassification <- rbind(bootstrapClassification, permData)
  }
  bootstrapClassification <- bootstrapClassification %>%
    mutate(ElectrodeID = electrode)
  allBootstrapClassification <- rbind(allBootstrapClassification, bootstrapClassification)
  
  # get corrected P values
  corrPResults <- data.frame()
  for (thisModel in 1:nlevels(classificationResults$Model)) {
    model <- levels(classificationResults$Model)[thisModel]
    permData <- bootstrapClassification %>%
      filter(Model == model) %>%
      select(TVal)
    trueData <- meanClassification %>%
      filter(Model == model) %>%
      select(TVal)
    permData <- permData %>%
      rbind(trueData) %>%
      arrange(desc(TVal))
    trueP <- min(which(permData$TVal == trueData$TVal)) / nrow(permData)
    thisResult <- data.frame(Model = model, CorrP = trueP, ElectrodeID = electrode)
    corrPResults <- rbind(corrPResults, thisResult)
  }
  meanClassification <- inner_join(meanClassification, corrPResults, by = "Model")
  allMeanClassificationResults <- rbind(allMeanClassificationResults, meanClassification)
}

allBootstrapClassification <- allBootstrapClassification %>%
  group_by(Model, ElectrodeID) %>%
  arrange(desc(TVal))
allBootstrapClassification$ElectrodeID <- factor(allBootstrapClassification$ElectrodeID)

# determine how many electrodes are significant by chance
maxElectrodePermutations <- data.frame()
for (thisIteration in 1:numIterations) {
  permData <- allBootstrapClassification %>%
    sample_n(1)
  permPVals <- data.frame()
  for (thisElectrode in 1:nlevels(permData$ElectrodeID)) {
    electrode <- levels(permData$ElectrodeID)[thisElectrode]
    for (thisModel in 1:nlevels(permData$Model)) {
      model <- levels(permData$Model)[thisModel]
      thisResult <- getCorrectedP(permData, allBootstrapClassification, electrode, model)
      permPVals <- rbind(permPVals, thisResult)
    }
  }
  permPVals <- permPVals %>%
    group_by(Model) %>%
    filter(CorrP < 0.05) %>%
    summarise(NSigElectrodes = n())
  thisMaxElec <- data.frame(NSigElectrodes = max(permPVals$NSigElectrodes))
  maxElectrodePermutations <- rbind(maxElectrodePermutations, thisMaxElec)
}

maxElectrodeP <- sortPerms(maxElectrodePermutations) %>%
  unique()

sigClassification <- allMeanClassificationResults %>%
  group_by(Model) %>%
  filter(CorrP < 0.05) %>%
  summarise(NSigElectrodes = n())

save(file = 'Rda/allClassificationResults.Rda', list = c('allClassificationResults', 'allMeanClassificationResults', 'maxElectrodeP', 'sigClassification', 'maxElectrodePermutations', 'allBootstrapClassification'))