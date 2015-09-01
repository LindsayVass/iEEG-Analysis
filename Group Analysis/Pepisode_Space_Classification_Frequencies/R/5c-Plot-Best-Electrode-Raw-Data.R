# Script Name:  5c-Plot-Best-Electrode-Raw-Data.R
# Author:       Lindsay Vass
# Date:         1 September 2015
# Purpose:      This script will produce plots of the raw data, with highlighted
#               pepisode, for the trials identified in 4c-Plot-Best-Electrode-Data.R
#               and exported by ../m/GetRawDataBestPepisodeTrials.m

library(dplyr)
library(ggplot2)
library(ggthemes)
library(R.matlab)


# Functions ---------------------------------------------------------------
plotLabeller <- function(variable, value) {
  variable <- spaceNames[value]
}

# Analysis ----------------------------------------------------------------
load('Rda/bestPepisodeTrials.Rda')
rawData <- readMat('mat/RawDataBestPepisodeTrials.mat')
allEEGData <- rawData$allEEGData
allBinaryData <- rawData$allBinaryData
eegTimes <- rawData$eegTimes

bestData$ElectrodeID <- factor(bestData$ElectrodeID)
bestData$TrialSpaceType <- factor(bestData$TrialSpaceType)

spaceNames <- c('Short', 'Long')

for (thisElec in 1:nlevels(bestData$ElectrodeID)) {
  electrode <- levels(bestData$ElectrodeID)[thisElec]
  
  nsInd <- which(bestData$ElectrodeID == electrode & bestData$TrialSpaceType == "NS")
  fsInd <- which(bestData$ElectrodeID == electrode & bestData$TrialSpaceType == "FS")
  
  allEEG <- data.frame(TrialSpaceType = c(rep('NS', length(unlist(allEEGData[[nsInd]]))), rep('FS', length(unlist(allEEGData[[fsInd]])))),
                       EEG = c(unlist(allEEGData[[nsInd]]), unlist(allEEGData[[fsInd]])),
                       Pepisode = c(unlist(allBinaryData[[nsInd]]), unlist(allBinaryData[[fsInd]])),
                       Time = as.vector(eegTimes))
   
  
  allEEG$TrialSpaceType <- factor(allEEG$TrialSpaceType, levels = c("NS", "FS"))
  
  ggplot(allEEG, aes(x = Time, y = EEG)) +
    geom_line(aes(colour = Pepisode)) +
    facet_grid(TrialSpaceType ~ ., labeller = plotLabeller) +
    scale_colour_gradient(low = "black", high = "red") +
    theme_few() +
    guides(colour = FALSE) +
    labs(x = 'Time in Teleporter (ms)',
         y = expression(paste('Voltage (', mu, 'V)'))) +
    theme(text = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.35),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA)) 
  
  ggsave(paste0('Figures/Single Electrode/RawTrace_', electrode, '.png'))
}