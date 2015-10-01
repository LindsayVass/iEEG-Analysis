# Script Name:  5b-Plot-Single-Trial-Data.R
# Author:       Lindsay Vass
# Date:         25 September
# Purpose:      This script will load the raw data from matlab and plot single
#               trial data, showing the raw traces, highlighted with pepisode,
#               a bar chart to the side showing the mean pepisode for each trial 
#               segment, and PSDs for navigation and teleportation time points.


library(R.matlab)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)


# Functions ---------------------------------------------------------------

moveStripToTop <- function(plot) {
  
  # Get the gtable
  gt <- ggplotGrob(plot)
  
  # Get the position of the panels in the layout
  panels <-c(subset(gt$layout, name=="panel", se=t:r))
  
  # Add a row above each panel
  for(i in rev(panels$t-1)) gt = gtable_add_rows(gt, unit(1, "lines"), i)
  
  # Get the positions of the panels and the strips in the revised layout
  panels <-c(subset(gt$layout, name=="panel", se=t:r))
  strips <- c(subset(gt$layout, name=="strip-right", se=t:r))
  
  # Get the strip grobs
  stripText = gtable_filter(gt, "strip-right")
  
  # Insert the strip grobs into the new rows
  for(i in 1:length(strips$t)) gt = gtable_add_grob(gt, stripText$grobs[[i]],  t=panels$t[i]-1, l=4, r=4)
  
  # Remove the old strips
  gt = gt[,-5]
  
  # For this plot - adjust the heights of the strips and the empty row above the strips
  for(i in panels$t) {
    gt$heights[i-1] = list(unit(2, "lines"))
    gt$heights[i-2] = list(unit(0.2, "lines"))
  }
  
  return(gt)
}

# Import and prepare data -------------------------------------------------

matData           <- readMat('mat/bestSingleTrialsRawData.mat')
allNavEEGData     <- matData$allNavEEGData
allTeleEEGData    <- matData$allTeleEEGData
allNavBinaryData  <- matData$allNavBinaryData
allTeleBinaryData <- matData$allTeleBinaryData
allNavEEGTimes    <- matData$allNavEEGTimes
allTeleEEGTimes   <- matData$allTeleEEGTimes

load('Rda/bestSingleTrials.Rda')
load('Rda/allCleanPowerData.Rda')


# Filter power data -------------------------------------------------------

plotSubs <- manualBestTrials %>%
  select(ElectrodeID) %>%
  unique()
allPower <- allPower %>%
  inner_join(plotSubs) %>%
  group_by(ElectrodeID, TimePoint, Frequency) %>%
  summarise(MeanTelePower = mean(telePower),
            MeanNavPower  = mean(navPower))
allPower$ElectrodeID <- factor(allPower$ElectrodeID)

# Clean up raw data -------------------------------------------------------

allNavData  <- vector(mode = "list", length = nrow(manualBestTrials))
allTeleData <- vector(mode = "list", length = nrow(manualBestTrials))
for (i in 1:nrow(manualBestTrials)) {
  allNavData[[i]] <- data.frame(ElectrodeID = manualBestTrials$ElectrodeID[i],
                                RealTrialNumber = manualBestTrials$RealTrialNumber[i],
                                Frequency = manualBestTrials$Frequency[i],
                                TrialSpaceType = manualBestTrials$TrialSpaceType[i],
                                TrialTimeType = manualBestTrials$TrialTimeType[i],
                                TimePoint = manualBestTrials$TimePoint[i],
                                Condition = "Navigation",
                                Time = unlist(allNavEEGTimes[[i]]),
                                EEG = unlist(allNavEEGData[[i]]),
                                Pepisode = unlist(allNavBinaryData[[i]]))
  allTeleData[[i]] <- data.frame(ElectrodeID = manualBestTrials$ElectrodeID[i],
                                 RealTrialNumber = manualBestTrials$RealTrialNumber[i],
                                 Frequency = manualBestTrials$Frequency[i],
                                 TrialSpaceType = manualBestTrials$TrialSpaceType[i],
                                 TrialTimeType = manualBestTrials$TrialTimeType[i],
                                 TimePoint = manualBestTrials$TimePoint[i],
                                 Condition = "Teleporter",
                                 Time = unlist(allTeleEEGTimes[[i]]),
                                 EEG = unlist(allTeleEEGData[[i]]),
                                 Pepisode = unlist(allTeleBinaryData[[i]]))
}

allNavData <- rbindlist(allNavData)
allTeleData <- rbindlist(allTeleData)
allData <- rbind(allNavData, allTeleData) %>%
  mutate(ObservationID = paste(ElectrodeID, Frequency, RealTrialNumber, sep = "_"))
allData$ObservationID <- factor(allData$ObservationID)

allData$TimePoint <- mapvalues(allData$TimePoint, from = c("Pre1", "Tele", "Post1"), to = c("Pre", "Tele", "Post"))
allData$Condition <- mapvalues(allData$Condition, from = c("Navigation", "Teleporter"), to = c("Navigation", "Teleportation"))



# Make plots --------------------------------------------------------------

for (i in 1:nlevels(allData$ObservationID)) {
  thisData <- allData %>%
    filter(ObservationID == levels(allData$ObservationID)[i])
  thisPowerData <- allPower %>%
    filter(ElectrodeID == thisData$ElectrodeID[1]) %>%
    melt(id.vars = c("ElectrodeID", "TimePoint", "Frequency"), measure.vars = c("MeanTelePower", "MeanNavPower"), variable.name = "Condition", value.name = "MeanPower")
  thisPowerData$Condition <- mapvalues(thisPowerData$Condition, from = c("MeanNavPower", "MeanTelePower"), to = c("Navigation", "Teleportation"))
  thisPowerData$Condition <- factor(thisPowerData$Condition, levels = c("Navigation", "Teleportation"))
  thisPowerData$TimePoint <- mapvalues(thisPowerData$TimePoint, from = c("Pre1", "Tele", "Post1"), to = c("Pre", "Tele", "Post"))
  
  if (thisData$TrialTimeType[1] == "NT") {
    x = c(0, 1830)
  }  else {
    x = c(0, 2830)
  }
  
  # make plots of raw data
  pLine <- thisData %>%
    ggplot(aes(x = Time, y = EEG)) +
    geom_line(aes(colour = Pepisode)) +
    scale_colour_gradient(low = "black", high = "red") +
    facet_grid(Condition ~ ., space = "free") + 
    geom_vline(xintercept = x, colour = "blue", size = 2) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(0.5, "lines"),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey", colour = "black"),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(face = "plain"),
          axis.title.y = element_text(vjust = 1)) +
    labs(x = "Time (ms)",
         y = expression(paste('Voltage (', mu, 'V)'))) +
    guides(colour = FALSE)
  pLine
  
  gLine <- moveStripToTop(pLine)
  
  # Draw it
  grid.newpage()
  grid.draw(gLine)
  
  
  # plot bar graph showing mean pepisode at each timepoint
  meanData <- thisData %>%
    group_by(TimePoint, Condition) %>%
    summarise(Pepisode = mean(Pepisode))
  
  pBar <- meanData %>%
    ggplot(aes(x = TimePoint, y = Pepisode)) +
    geom_bar(stat = "identity") +
    facet_grid(Condition ~ .) +
    labs(y = expression(P[Episode]),
         x = "Time Point") +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(0.5, "lines"),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey", colour = "black"),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(face = "plain"),
          axis.title.y = element_text(vjust = 1)) +
    scale_y_continuous(limits = c(0,1))
  
  gBar <- moveStripToTop(pBar)
  grid.newpage()
  grid.draw(gBar)
  
  # plot PSDs
  freqList <- unique(thisPowerData$Frequency)
  freqInds <- seq(1, length(freqList), 5)
  freqList <- freqList[freqInds]
  
  pPSD <- thisPowerData %>%
    ggplot(aes(x = Frequency, y = MeanPower)) +
    geom_line(aes(color = TimePoint), size = 1) +
    facet_grid(Condition ~ .) +
    scale_x_log10(breaks = freqList, labels = signif(freqList, 2)) +
    scale_y_log10() +
    labs(x = "Frequency (Hz)",
         y = expression(paste('Power (', mu, 'V'^'2',' / Hz)'))) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(0.5, "lines"),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey", colour = "black"),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(face = "plain"),
          axis.title.y = element_text(vjust = 1))
    
  gPSD <- moveStripToTop(pPSD)
  grid.newpage()
  grid.draw(gPSD)
  
  # save combined figure
  electrodeID <- unique(thisData$ElectrodeID)
  freq <- unique(thisData$Frequency)
  trial <- unique(thisData$RealTrialNumber)
  fileName <- paste0('Figures/RawTracePlusPepisode/', electrodeID, '_', freq, 'Hz_Trial', trial, '.pdf')
  pdf(file = fileName, width = 20, height = 6)
  #grid.arrange(gt, pBar, ncol = 2, widths = c(6, 1.5))
  grid.arrange(gt, gBar, gPSD, ncol = 3, widths = c(50, 20, 30))
  dev.off()
}


