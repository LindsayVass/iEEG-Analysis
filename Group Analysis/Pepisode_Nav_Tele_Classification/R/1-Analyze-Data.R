# Script Name:  1-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         30 September 2015
# Purpose:      This script will test whether navigation and teleportation can 
#               be distinguished using the frequency-specific pattern of 
#               pepisode values, using all frequencies <= 8 Hz. It will load the
#               data from the Pepisode_Sustainedness folder and run analyses as
#               in Pepisode_Space_Classification_Frequencies.

library(dplyr)
library(kernlab)
library(reshape2)
library(permute)

load('../Pepisode_Sustainedness/Rda/allCleanData.Rda')
remove(navSustain, teleSustain, sessionInfo)

# Functions ---------------------------------------------------------------

runClassifier <- function(trainingData, testingData) {
  
  tryCatch(fit <- ksvm(Condition ~ ., data = trainingData),
           error = function(e) {
             fit <- ksvm(Condition ~ ., data = trainingData)
           })
  predictions <- predict(fit, testingData[, 2:ncol(testingData)], type = "response")
  confusionMatrix <- table(predictions, testingData$Condition)
  accuracy <- sum(diag(confusionMatrix)) / sum(confusionMatrix)
}

filterSample <- function(thisData, condition, trainSize, testSize) {
 
  trialList <- thisData %>%
    filter(Condition == condition) %>%
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
    dcast(RealTrialNumber + Condition ~ Frequency, value.var = "Pepisode") %>%
    select(-(RealTrialNumber))
  
}

classifyOnCondition <- function(thisData, trainingSize, testingSize) {
  
  # Extract training and testing lists of trials
  navLists <- filterSample(thisData, "Navigation", trainingSize, testingSize)
  teleLists <- filterSample(thisData, "Teleportation", trainingSize, testingSize)
  
  navTrain <- getSampledData(thisData, navLists$train)
  navTest  <- getSampledData(thisData, navLists$test)
  teleTrain <- getSampledData(thisData, teleLists$train)
  teleTest  <- getSampledData(thisData, teleLists$test)
  
  # cast to wide
  navTrainWide <- castWide(navTrain)
  navTestWide  <- castWide(navTest)
  teleTrainWide <- castWide(teleTrain)
  teleTestWide  <- castWide(teleTest)
  
  # concatenate
  trainingData <- rbind(navTrainWide, teleTrainWide)
  testingData  <- rbind(navTestWide, teleTestWide)
  
  accuracy <- runClassifier(trainingData, testingData)
}

runClassification <- function(inputData, trainingSize, testingSize, electrode, model, iteration) {
  
  accuracy <- classifyOnCondition(inputData, trainingSize, testingSize)
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

cleanData <- allPepisode %>%
  filter(Frequency >= minFreq & Frequency <= 8,
         TimePoint == "Tele") %>%
  ungroup() %>%
  select(-c(TrialType:TimePoint, FrequencyBand))
cleanData$Frequency <- factor(cleanData$Frequency)

cleanData <- melt(cleanData, id.vars = c('ElectrodeID', 'RealTrialNumber', 'TrialSpaceType', 'TrialTimeType', 'Frequency'),
                  measure.vars = c('TelePepisode', 'NavPepisode'),
                  variable.name = 'Condition',
                  value.name = 'Pepisode') %>%
  mutate(TrialID = paste(ElectrodeID, RealTrialNumber, Condition, sep = "_")) 
cleanData$Condition <- plyr::mapvalues(cleanData$Condition, from = c('NavPepisode', 'TelePepisode'), to = c('Navigation', 'Teleportation'))


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
    filter(ElectrodeID == electrode)
  
  numObservations <- thisData %>%
    group_by(Condition) %>%
    summarise(Count = n() / nlevels(Frequency))
  
  numTimeObservations <- thisData %>%
    group_by(TrialTimeType, Condition) %>%
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