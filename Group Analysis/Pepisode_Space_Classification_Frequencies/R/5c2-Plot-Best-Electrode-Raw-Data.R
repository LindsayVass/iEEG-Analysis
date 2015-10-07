# Script Name:  5c-Plot-Best-Electrode-Raw-Data.R
# Author:       Lindsay Vass
# Date:         1 September 2015
# Purpose:      This script will produce plots of the raw data for the trials 
#               identified in 4c-Plot-Best-Electrode-Data.R
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
plotLabeller <- function(variable, value) {
  variable <- spaceNames[value]
}

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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

runClassification <- function(classData, thisResult, numPerm = 1000, trainingPercent = 0.75) {
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
  oneElecResults <- data.table::rbindlist(oneElecResults)
  
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
  accuracyLabel <- paste0(trialAccuracy$Count[which(trialAccuracy$Accurate == TRUE)], '/', sum(trialAccuracy$Count), ' trials')
  
  return(list(meanTrialResults = meanTrialResults, accuracyLabel = accuracyLabel))
}

# Prep Raw Data ----------------------------------------------------------------
load('Rda/bestPepisodeTrials.Rda')

#rawData <- readMat('mat/RawDataBestPepisodeTrials.mat')
rawData <- readMat('mat/RawDataBestPepisodeTrials_1-8Hz.mat')
allEEGData <- rawData$allEEGData
allBinaryData <- rawData$allBinaryData
eegTimes <- rawData$eegTimes
frequencies <- rawData$frequencies

bestData$ElectrodeID <- factor(bestData$ElectrodeID)
bestData$TrialSpaceType <- factor(bestData$TrialSpaceType)

spaceNames <- c('Short', 'Long')

# parameters for plotting lines below EEG showing pepisode for the non-plotted
# frequencies. 
linePlotScaleMin <- 1.2
linePlotInterval <- 0.1
linePlotScaleMax <- linePlotScaleMin + linePlotInterval * (length(frequencies) - 1)
linePlotScales   <- seq(from = linePlotScaleMax, to = linePlotScaleMin, by = -linePlotInterval)


# Prep electrode data -----------------------------------------------------

load('Rda/allClassificationResults.Rda')
load('Rda/allCleanData.Rda')

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


# Make plots --------------------------------------------------------------

stripColors <- brewer.pal(4, "Set1")
stripColors <- stripColors[3:4]

for (thisElec in 1:nlevels(bestData$ElectrodeID)) {
  
  ########### RAW DATA PLOT #############
  electrode <- levels(bestData$ElectrodeID)[thisElec]
  
  nsInd <- which(bestData$ElectrodeID == electrode & bestData$TrialSpaceType == "NS")
  fsInd <- which(bestData$ElectrodeID == electrode & bestData$TrialSpaceType == "FS")
  
  nsBinary <- t(as.data.frame(allBinaryData[[nsInd]])) %>%
    as.data.frame()
  names(nsBinary) <- round(frequencies, digits = 2)
  fsBinary <- t(as.data.frame(allBinaryData[[fsInd]])) %>%
    as.data.frame()
  names(fsBinary) <- round(frequencies, digits = 2)
  allBinary <- rbind(nsBinary, fsBinary)
  
  plotColorFreq <- round(bestData$Frequency[nsInd], digits = 2)
  plotColorInd  <- which(names(nsBinary) == plotColorFreq)
  
  allEEG <- data.frame(TrialSpaceType = c(rep('NS', length(unlist(allEEGData[[nsInd]]))), rep('FS', length(unlist(allEEGData[[fsInd]])))),
                       EEG = c(unlist(allEEGData[[nsInd]]), unlist(allEEGData[[fsInd]])),
                       Time = as.vector(eegTimes))
  
  minEEG <- min(allEEG$EEG)
  scaleMatrix <- t(replicate(nrow(allBinary), linePlotScales)) * minEEG
  scaledBinary <- allBinary * scaleMatrix
  
  allEEG$TrialSpaceType <- factor(allEEG$TrialSpaceType, levels = c("NS", "FS"))
  allEEG <- cbind(allEEG, Pepisode = allBinary[, plotColorInd])
  
  scaledBinary[scaledBinary == 0] <- NA  
  scaledBinary <- cbind(Time = allEEG$Time, TrialSpaceType = allEEG$TrialSpaceType,  scaledBinary) %>%
    melt(id = c("Time", "TrialSpaceType"), variable.name = "Frequency",  value.name = "y") 
  scaledBinary$Frequency <- as.numeric(levels(scaledBinary$Frequency))[scaledBinary$Frequency]
   
  yVal <- data.frame(Frequency = t(round(frequencies, digits = 2)), y = minEEG * linePlotScales)
  
  pepisodeVals <- scaledBinary %>%
    group_by(Frequency, TrialSpaceType) %>%
    mutate(Pepisode = ifelse(is.na(y) == TRUE, 0, 1)) %>%
    summarise(MeanPepisode = round(mean(Pepisode), digits = 2)) %>%
    mutate(Label = paste0(sprintf(Frequency, fmt = "%#.3g"), " Hz (", MeanPepisode, ")"),
           x = 1.01 * max(allEEG$Time)) %>%
    inner_join(yVal)
  
  rawTrace <- ggplot(allEEG, aes(x = Time, y = EEG)) +
    geom_line(colour = "black") +
    facet_grid(TrialSpaceType ~ ., labeller = plotLabeller) +
    theme_few() +
    guides(colour = FALSE, fill = FALSE) +
    labs(x = 'Time in Teleporter (ms)',
         y = expression(paste('Voltage (', mu, 'V)'))) +
    scale_shape_identity() +
    geom_point(data = scaledBinary,
               colour = NA,
               aes(x = Time,
                   y = y,
                   shape = 22,
                   fill = Frequency)) +
    expand_limits(x = c(0, 1.07 * max(allEEG$Time))) +
    scale_fill_continuous(limits = c(signif(minFreq, 2), signif(maxFreq,2)), breaks =  c(signif(minFreq, 2), signif(maxFreq,2))) +
    guides(fill = guide_colorbar(label.vjust = -0.05)) + 
    theme(text = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.35),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.margin = unit(1, "lines"),
          strip.background = element_rect(colour = "black", fill = stripColors[1]),
          strip.text = element_text(colour = "white"),
          legend.position = c(0.95, 0.1),
          legend.title = element_blank(),
          legend.key.height = unit(0.22, 'cm'),
          legend.background = element_blank())
  rawTrace
  
  ######## Make scatterplot of mean pepisode at each frequency #############
  meanFreqPepisode <- classData %>%
    filter(ElectrodeID == as.character(topClass$ElectrodeID[thisElec])) %>%
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
    scale_colour_manual(values = stripColors) +
    scale_x_log10(breaks = unique(meanFreqPepisode$Frequency), labels = round(unique(meanFreqPepisode$Frequency), digits = 2)) + 
    theme_few() +
    theme(text = element_text(size = 18),
          axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = -0.25),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.position = c(0.27, 0.85)) +
    labs(x = 'Frequency (Hz)',
         y = expression('Mean P'['Episode']),
         colour = 'Distance Teleported')
  scatterP
  
  ######### Make bar plot of trialwise classification #############
  
  classResults <- runClassification(classData, thisElec)
  
  trialClassP <- classResults$meanTrialResults %>%
    transform(TrialNumber = reorder(TrialNumber, PropNS)) %>%
    ggplot(aes(x = factor(TrialNumber),
               y = PropNS,
               fill = TrialSpaceType)) +
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
          legend.position = c(0.25, 0.85),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(vjust = -0.5)) +
    labs(y = '% Iterations Classified as "Short"',
         x = 'Trial',
         fill = 'Distance Teleported') 
  trialClassP
  
  ######## Plot all subplots together #########
  fileName <- paste0('Figures/Single Electrode/MultiPlot/', electrode, '.pdf')
  pdf(file = fileName, width = 10, height = 8, useDingbats = FALSE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(rawTrace, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  print(scatterP, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(trialClassP, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  dev.off()
}