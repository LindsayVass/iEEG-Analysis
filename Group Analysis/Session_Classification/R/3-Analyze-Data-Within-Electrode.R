# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         25 November 2015
# Purpose:      This script will analyze the data output by 2-Clean-Data.R. It 
#               will restrict analysis to electrodes with recordings from 2 
#               test sessions and test whether PEpisode during Pre, Tele, or
#               Post can be used to classify trials based on session.
#               This analysis will retain the frequency-specific information and
#               use all frequencies <= 8 Hz to create the pattern. It will perform 
#               all statistics within-electrode and determine whether there are 
#               more significant electrodes than would be expected by chance.

library(dplyr)
library(kernlab)
library(reshape2)
library(permute)

load('Rda/allCleanData.Rda')


# Functions ---------------------------------------------------------------

runClassifier <- function(trainingData, testingData) {
  
  tryCatch(fit <- ksvm(Session ~ ., data = trainingData),
           error = function(e) {
             fit <- ksvm(Session ~ ., data = trainingData)
           })
  predictions <- predict(fit, testingData[, 2:ncol(testingData)], type = "response")
  confusionMatrix <- table(predictions, testingData$Session)
  accuracy <- sum(diag(confusionMatrix)) / sum(confusionMatrix)
}

filterSample <- function(thisData, session, trainSize, testSize) {
  
  trialList <- thisData %>%
    filter(Session == session) %>%
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
    dcast(TrialNumber + Session ~ Frequency, value.var = "Pepisode") %>%
    select(-(TrialNumber))
  
}

classifyOnSession <- function(thisData, trainingSize, testingSize) {
  
  # Extract training and testing lists of trials
  ALists <- filterSample(thisData, "TeleporterA", trainingSize, testingSize)
  BLists <- filterSample(thisData, "TeleporterB", trainingSize, testingSize)
  
  ATrain <- getSampledData(thisData, ALists$train)
  ATest  <- getSampledData(thisData, ALists$test)
  BTrain <- getSampledData(thisData, BLists$train)
  BTest  <- getSampledData(thisData, BLists$test)
  
  # cast to wide
  ATrainWide <- castWide(ATrain)
  ATestWide  <- castWide(ATest)
  BTrainWide <- castWide(BTrain)
  BTestWide  <- castWide(BTest)
  
  # concatenate
  trainingData <- rbind(ATrainWide, BTrainWide)
  testingData  <- rbind(ATestWide, BTestWide)
  
  accuracy <- runClassifier(trainingData, testingData)
}

runClassification <- function(inputData, trainingSize, testingSize, electrode, model, iteration) {
  
  accuracy <- classifyOnSession(inputData, trainingSize, testingSize)
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

# keep only electrodes with 2 sessions of data
cleanData <- tidyr::separate(cleanData, ElectrodeID, c('Subject', 'Session', 'Electrode')) %>%
  mutate(ElectrodeName = paste(Subject, Electrode, sep = "_"))
elecList <- cleanData %>%
  select(ElectrodeName, Session) %>%
  unique() %>%
  group_by(ElectrodeName) %>%
  summarise(Count = n()) %>%
  filter(Count == 2) %>%
  select(-Count)
cleanData <- inner_join(cleanData, elecList)
cleanData$ElectrodeName <- factor(cleanData$ElectrodeName)

numIterations <- 1000

allClassificationResults <- data.frame()
allMeanClassificationResults <- data.frame()
allBootstrapClassification <- data.frame()
#imeLevels  <- c('Both', 'NT', 'FT')
#spaceLevels <- c('Both', 'NS', 'FS')
#timeBins  <- c('Pre', 'Tele', 'Post')
timeLevels <- 'Both'
spaceLevels <- 'Both'
timeBins <- 'Tele'
for (thisElectrode in 1:nlevels(cleanData$ElectrodeName)) {
  classificationResults <- data.frame()
  permutationResults    <- data.frame()
  
  electrode <- levels(cleanData$ElectrodeName)[thisElectrode]
  
  cat(paste0('\n\nElectrode ', thisElectrode, ' of ', nlevels(cleanData$ElectrodeName)))
  
  for (thisTime in 1:length(timeLevels)) {
    
    cat(paste0('\nTime Type ', thisTime, ' of ', length(timeLevels)))
    
    timeType <- timeLevels[thisTime]
    
    for (thisSpace in 1:length(spaceLevels)) {
      
      cat(paste0('\nSpace Type ', thisSpace, ' of ', length(spaceLevels)))
      
      spaceType <- spaceLevels[thisSpace]
      
      for (thisBin in 1:length(timeBins)) {
        
        cat(paste0('\nTime Bin ', thisBin, ' of ', length(timeBins)))
        
        timeBin <- timeBins[thisBin]
        
        thisData <- cleanData %>%
          filter(ElectrodeName == electrode)
        
        if (timeType != 'Both') {
          thisData <- thisData %>%
            filter(TrialTimeType == timeType)
        }
        
        if (spaceType != 'Both') {
          thisData <- thisData %>%
            filter(TrialSpaceType == spaceType)
        }
        
        thisData <- thisData %>%
          filter(TimeBin == timeBin)
        
        numObservations <- thisData %>%
          group_by(Session) %>%
          summarise(Count = n() / nlevels(Frequency))
        
        # percentage of the data set to use to train the classifier
        trainingPercent  <- 0.75 
        trainingSize     <- floor(trainingPercent * min(numObservations$Count))
        testingSize      <- min(numObservations$Count) - trainingSize
        
        if (testingSize < 3) {
          next
        }
        
        thisModel <- paste0('S-', spaceType, '_T-', timeType, '_', timeBin)
        
        for (thisIteration in 1:numIterations) {
          log <- capture.output({
            thisResult <- suppressWarnings(runClassification(thisData, trainingSize, testingSize, electrode, thisModel, thisIteration))
          })
          classificationResults <- rbind(classificationResults, thisResult)
          
          permData <- shuffleData(thisData)
          log <- capture.output({
            thisPermResult <- suppressWarnings(runClassification(permData, trainingSize, testingSize, electrode, thisModel, thisIteration))
          })
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
      }
    }
  }
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