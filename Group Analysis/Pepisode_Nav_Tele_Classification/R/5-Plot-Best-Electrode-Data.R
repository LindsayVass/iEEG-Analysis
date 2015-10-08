# Script Name:  5-Plot-Best-Electrode-Data.R
# Author:       Lindsay Vass
# Date:         8 October 2015
# Purpose:      This script will produce plots of the raw data for the trials 
#               identified in 4-Select-Best-Trials.R
#               and exported by ../m/GetRawDataBestPepisodeTrials.m
#               It will also plot mean pepisode at each frequency, and bar plots
#               showing the classification for each trial across iterations.

library(dplyr)
library(ggplot2)
library(ggthemes)
library(R.matlab)
library(RColorBrewer)
library(reshape2)
library(grid)
library(kernlab)
library(gridExtra)

# Functions ---------------------------------------------------------------
runClassIterations <- function (thisData, electrode, model = "Both", numIterations = 1000, trainingPercent = 0.75) {
  numObservations <- thisData %>%
    group_by(Condition) %>%
    summarise(Count = n() / nlevels(Frequency))
  trainingSize     <- floor(trainingPercent * min(numObservations$Count))
  testingSize      <- min(numObservations$Count) - trainingSize
  
  classificationResults <- vector(mode = "list", length = numIterations)
  for (thisIteration in 1:numIterations) {
    #### Run classification using both NT and FT ####
    thisResult <- classifyOnCondition(thisData, trainingSize, testingSize, electrode)
    classificationResults[[thisIteration]] <- thisResult
  }
  classificationResults <- data.table::rbindlist(classificationResults)
  
  meanTrialResults <- classificationResults %>% 
    select(-ElectrodeID) %>%
    group_by(TrialNumber, Condition, Classification) %>%
    summarise(Count = n()) %>%
    as.data.frame() %>%
    arrange(TrialNumber) %>%
    ungroup() %>%
    group_by(TrialNumber) 
  zeroResult <- meanTrialResults %>%
    select(Condition, TrialNumber) %>%
    unique()
  zeroNav <- zeroResult %>%
    mutate(Classification = "Navigation")
  zeroTele <- zeroResult %>%
    mutate(Classification = "Teleportation")
  zeroResult <- rbind(zeroNav, zeroTele) %>%
    anti_join(meanTrialResults) %>%
    mutate(Count = 0) %>%
    as.data.frame()
  meanTrialResults <- as.data.frame(meanTrialResults) %>%
    rbind(zeroResult) %>%
    arrange(TrialNumber, Classification)
  trialInd <- strsplit(as.character(meanTrialResults$TrialNumber), "_")
  trialIndNum <- vector(length = length(trialInd))
  for (i in 1:length(trialInd)) {
    trialIndNum[i] <- trialInd[[i]][1]
  }
  trialIndNum <- data.frame(TrialInd = trialIndNum)
  meanTrialResults <- cbind(meanTrialResults, trialIndNum) %>%
    group_by(TrialInd, Condition) %>%
    mutate(TotalCount = sum(Count),
           PropNav = Count / TotalCount - 0.5) %>%
    filter(Classification == "Navigation")

  
  trialAccuracy <- meanTrialResults %>%
    select(TrialNumber, Condition, PropNav) %>%
    unique() %>%
    mutate(FinalClass = ifelse(PropNav > 0, "Navigation", "Teleportation"),
           Accurate = ifelse(Condition == FinalClass, TRUE, FALSE)) %>%
    group_by(Accurate) %>%
    summarise(Count = n())
  accuracyLabel <- paste0(trialAccuracy$Count[which(trialAccuracy$Accurate == TRUE)], '/', sum(trialAccuracy$Count), ' trials')
  
  return(list(meanTrialResults = meanTrialResults, accuracyLabel = accuracyLabel))
}

classifyOnCondition <- function(thisData, trainingSize, testingSize, electrode) {
  
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
  
  testTrials <- rbind(navTest, teleTest) %>%
    mutate(TrialName = paste(RealTrialNumber, Condition, sep = "_")) %>%
    select(TrialName, Condition) %>%
    unique()
  
  accuracy <- runClassifier(trainingData, testingData, testTrials, electrode)
}

runClassifier <- function(trainingData, testingData, testTrials, electrode) {
  
  tryCatch(fit <- ksvm(Condition ~ ., data = trainingData),
           error = function(e) {
             fit <- ksvm(Condition ~ ., data = trainingData)
           })
  predictions <- data.frame(Classification = predict(fit, testingData[, 2:ncol(testingData)], type = "response"),
                            ElectrodeID = electrode,
                            TrialNumber = testTrials$TrialName,
                            Condition = testTrials$Condition)
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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# Prep Raw Data ----------------------------------------------------------------
load('Rda/bestPepisodeTrials.Rda')
rawData <- readMat('mat/RawDataBestPepisodeTrials_1-8Hz.mat')

allTeleEEGData <- rawData$allTeleEEGData
allNavEEGData  <- rawData$allNavEEGData
allTeleBinaryData <- rawData$allTeleBinaryData 
allNavBinaryData  <- rawData$allNavBinaryData 
eegNtTimes <- rawData$eegNtTimes
eegFtTimes <- rawData$eegFtTimes
frequencies <- rawData$frequencies

bestTele$ElectrodeID <- factor(bestTele$ElectrodeID)
bestNav$ElectrodeID  <- factor(bestNav$ElectrodeID)

# parameters for plotting lines below EEG showing pepisode 
linePlotScaleMin <- 1.2
linePlotInterval <- 0.1
linePlotScaleMax <- linePlotScaleMin + linePlotInterval * (length(frequencies) - 1)
linePlotScales   <- seq(from = linePlotScaleMax, to = linePlotScaleMin, by = -linePlotInterval)

# Prep electrode data -----------------------------------------------------

load('Rda/allClassificationResults.Rda')
load('../Pepisode_Sustainedness/Rda/allCleanData.Rda')

minFreq <- 1000 / (1830 / 3)
maxFreq <- 8

classData <- allPepisode %>%
  inner_join(topClass, by = c("ElectrodeID")) %>%
  filter(TimePoint == 'Tele',
         Frequency > minFreq,
         Frequency <= maxFreq) %>%
  group_by(ElectrodeID)
classData$ElectrodeID <- factor(classData$ElectrodeID)

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

# Make plots --------------------------------------------------------------

stripColors <- brewer.pal(5, "Set1")
stripColors <- c(stripColors[2], stripColors[5])

for (thisElec in 1:nlevels(bestNav$ElectrodeID)) {
  
  electrode <- levels(bestNav$ElectrodeID)[thisElec]
  
  thisElectrodeNavData <- bestNav %>%
    filter(ElectrodeID == electrode)
  thisElectrodeTeleData <- bestTele %>%
    filter(ElectrodeID == electrode)
  
  ######## Make scatterplot of mean pepisode at each frequency #############
  meanFreqPepisode <- cleanData %>%
    filter(ElectrodeID == electrode) %>%
    group_by(Condition, Frequency) %>%
    summarise(MeanPepisode = mean(Pepisode),
              SEM = sd(Pepisode) / sqrt(n()))
  meanFreqPepisode$Frequency <- as.numeric(levels(meanFreqPepisode$Frequency))[meanFreqPepisode$Frequency]
  scatterP <- meanFreqPepisode %>%
    ggplot(aes(x = Frequency,
               y = MeanPepisode,
               ymin = MeanPepisode - SEM,
               ymax = MeanPepisode + SEM,
               color = Condition)) +
    geom_point(size = 5) +
    geom_pointrange() +
    scale_colour_manual(values = stripColors) +
    scale_x_log10(breaks = frequencies, labels = round(frequencies, digits = 2)) + 
    theme_few() +
    theme(text = element_text(size = 18),
          axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = -0.25),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.position = c(0.2, 0.85)) +
    labs(x = 'Frequency (Hz)',
         y = expression('Mean P'['Episode']),
         colour = 'Condition')
  scatterP
  
  ######### Make bar plot of trialwise classification #############
  
  thisData <- cleanData %>%
    filter(ElectrodeID == electrode)
  
  classResults <- runClassIterations(thisData, electrode)
  
  trialClassP <- classResults$meanTrialResults %>%
    transform(TrialNumber = reorder(TrialNumber, PropNav)) %>%
    ggplot(aes(x = factor(TrialNumber),
               y = PropNav,
               fill = Condition)) +
    scale_fill_manual(values = stripColors) +
    geom_bar(stat = 'identity', position = 'identity') +
    scale_y_continuous(limits = c(-0.5, 0.5),
                       breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                       labels = c(0, 25, 50, 75, 100)) +
    annotation_custom(grob = textGrob(classResults$accuracyLabel, gp = gpar(fontsize = 18)), xmin = -15, xmax = Inf, ymin = 0, ymax = Inf) +
    coord_flip() +
    theme_few() +
    theme(text = element_text(size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = c(0.2, 0.85),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(vjust = -0.5)) +
    labs(y = '% Iterations Classified as "Navigation"',
         x = 'Trial',
         fill = 'Condition') 
  trialClassP
  
  for (thisTrial in 1:nrow(thisElectrodeNavData)) {
    
    realTrialNumber <- thisElectrodeNavData$RealTrialNumber[thisTrial]
    
    trialTimeType <- thisElectrodeNavData$TrialTimeType[thisTrial]
    
    if (trialTimeType == "NT") {
      theTimes = eegNtTimes
    } else {
      theTimes = eegFtTimes
    }
    
    ########### RAW DATA PLOT #############
    
    navInd <- which(bestNav$ElectrodeID == electrode & bestNav$RealTrialNumber == realTrialNumber)
    teleInd <- which(bestTele$ElectrodeID == electrode & bestTele$RealTrialNumber == realTrialNumber)
    
    if (length(navInd) == 0 | length(teleInd) == 0)
      next
    
    navBinary <- t(as.data.frame(allNavBinaryData[[navInd]])) %>%
      as.data.frame()
    names(navBinary) <- round(frequencies, digits = 2)
    
    teleBinary <- t(as.data.frame(allTeleBinaryData[[teleInd]])) %>%
      as.data.frame()
    names(teleBinary) <- round(frequencies, digits = 2)
    
    allBinary <- rbind(navBinary, teleBinary)
    
    allEEG <- data.frame(Condition = c(rep('Navigation', length(unlist(allNavEEGData[[navInd]]))), rep('Teleportation', length(unlist(allTeleEEGData[[teleInd]])))),
                         EEG = c(unlist(allNavEEGData[[navInd]]), unlist(allTeleEEGData[[teleInd]])),
                         Time = as.vector(theTimes))
    
    minEEG <- min(allEEG$EEG)
    scaleMatrix <- t(replicate(nrow(allBinary), linePlotScales)) * minEEG
    scaledBinary <- allBinary * scaleMatrix
    
    allEEG$Condition <- factor(allEEG$Condition, levels = c("Navigation", "Teleportation"))
    
    scaledBinary[scaledBinary == 0] <- NA  
    scaledBinary <- cbind(Time = allEEG$Time, Condition = allEEG$Condition,  scaledBinary) %>%
      melt(id = c("Time", "Condition"), variable.name = "Frequency",  value.name = "y") 
    scaledBinary$Frequency <- as.numeric(levels(scaledBinary$Frequency))[scaledBinary$Frequency]
    
    yVal <- data.frame(Frequency = t(round(frequencies, digits = 2)), y = minEEG * linePlotScales)
    
    pepisodeVals <- scaledBinary %>%
      group_by(Frequency, Condition) %>%
      mutate(Pepisode = ifelse(is.na(y) == TRUE, 0, 1)) %>%
      summarise(MeanPepisode = round(mean(Pepisode), digits = 2)) %>%
      mutate(Label = paste0(sprintf(Frequency, fmt = "%#.3g"), " Hz (", MeanPepisode, ")"),
             x = 1.01 * max(allEEG$Time)) %>%
      inner_join(yVal)
    
    rawTrace <- ggplot(allEEG, aes(x = Time, y = EEG)) +
      geom_line(colour = "black") +
      facet_grid(Condition ~ .) +
      theme_few() +
      guides(colour = FALSE, fill = FALSE) +
      labs(x = 'Time (ms)',
           y = expression(paste('Voltage (', mu, 'V)'))) +
      scale_shape_identity() +
      geom_point(data = scaledBinary,
                 colour = NA,
                 aes(x = Time,
                     y = y,
                     shape = 22,
                     fill = Frequency)) +
      expand_limits(x = c(0, 1.2 * max(allEEG$Time))) +
      scale_fill_continuous(limits = c(signif(minFreq, 2), signif(maxFreq,2)), breaks =  c(signif(minFreq, 2), signif(maxFreq,2))) +
      guides(fill = guide_colorbar(label.vjust = -0.05)) + 
      theme(text = element_text(size = 18),
            axis.title.x = element_text(vjust = -0.35),
            axis.ticks = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA),
            panel.margin = unit(1, "lines"),
            strip.background = element_rect(colour = "black", fill = stripColors[1]),
            strip.text = element_text(colour = "white"),
            legend.position = c(0.9, 0.08),
            legend.title = element_blank(),
            legend.key.height = unit(0.2, 'cm'),
            legend.background = element_blank())
    rawTrace
    
    
    
    ######## Plot all subplots together #########
    fileName <- paste0('Figures/Single Electrode/MultiPlot/', electrode, '_', realTrialNumber, '.pdf')
    
    pdf(file = fileName, width = 16, height = 4, useDingbats = FALSE)
    grid.arrange(rawTrace, scatterP, trialClassP, ncol = 3)
    dev.off()
  }
}