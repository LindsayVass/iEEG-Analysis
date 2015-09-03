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
library(RColorBrewer)
library(reshape2)


# Functions ---------------------------------------------------------------
plotLabeller <- function(variable, value) {
  variable <- spaceNames[value]
}

# Analysis ----------------------------------------------------------------
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

for (thisElec in 1:nlevels(bestData$ElectrodeID)) {
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
 # scaleMatrix <- cbind(scaleMatrix[, 1:plotColorInd - 1], rep(1, nrow(allBinary)), scaleMatrix[, plotColorInd:ncol(scaleMatrix)])
  scaledBinary <- allBinary * scaleMatrix

  
  allEEG$TrialSpaceType <- factor(allEEG$TrialSpaceType, levels = c("NS", "FS"))
  allEEG <- cbind(allEEG, Pepisode = allBinary[, plotColorInd])
  
  scaledBinary[scaledBinary == 0] <- NA  
  scaledBinary <- cbind(Time = allEEG$Time, TrialSpaceType = allEEG$TrialSpaceType,  scaledBinary) %>%
    melt(id = c("Time", "TrialSpaceType"), variable.name = "Frequency",  value.name = "y") 
  scaledBinary$Frequency <- as.numeric(levels(scaledBinary$Frequency))[scaledBinary$Frequency]
 
 #   pepisodeLabels <- bestData %>%
 #     filter(ElectrodeID == electrode) %>%
 #     mutate(Label = paste("P[Episode](", round(Frequency, digits = 2), "~Hz", ") == ", round(Pepisode, digits = 2)),
 #            x = 0,
 #            y = min(allEEG$EEG))
 
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
#     scale_colour_gradient(low = "black", high = "red") +
    theme_few() +
    guides(colour = FALSE, fill = FALSE) +
    labs(x = 'Time in Teleporter (ms)',
         y = expression(paste('Voltage (', mu, 'V)'))) +
    theme(text = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.35),
          axis.ticks = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.margin = unit(1, "lines")) +
     scale_shape_identity() +
     geom_point(data = scaledBinary,
                colour = NA,
                aes(x = Time,
                    y = y,
                    shape = 22,
                    fill = Frequency)) +
    geom_text(data = pepisodeVals, aes(x = x, y = y, label = Label, hjust = 0), size = 2) +
    expand_limits(x = c(0, 1.07 * max(allEEG$Time)))
   

    
    ggsave(paste0('Figures/Single Electrode/RawTrace_', electrode, '.pdf'), width = 10.7)
}