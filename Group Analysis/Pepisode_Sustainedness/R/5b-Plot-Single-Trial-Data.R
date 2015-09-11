# Script Name:  5b-Plot-Single-Trial-Data.R
# Author:       Lindsay Vass
# Date:         9 September
# Purpose:      This script will load the raw data from matlab and plot single
#               trial data, showing the raw traces, highlighted with pepisode,
#               and a bar chart to the side showing the mean for each trial 
#               segment.


library(R.matlab)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(gridExtra)
library(cowplot)


# Import and prepare data -------------------------------------------------

matData           <- readMat('mat/bestSingleTrialsRawData.mat')
allNavEEGData     <- matData$allNavEEGData
allTeleEEGData    <- matData$allTeleEEGData
allNavBinaryData  <- matData$allNavBinaryData
allTeleBinaryData <- matData$allTeleBinaryData
allNavEEGTimes    <- matData$allNavEEGTimes
allTeleEEGTimes   <- matData$allTeleEEGTimes

load('Rda/bestSingleTrials.Rda')

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

for (i in 1:nlevels(allData$ObservationID)) {
  thisData <- allData %>%
    filter(ObservationID == levels(ObservationID)[i])
  
  if (thisData$TrialTimeType[1] == "NT") {
    x = c(0, 1830)
  }  else {
    x = c(0, 2830)
  }
  
  pLine <- thisData %>%
    ggplot(aes(x = Time, y = EEG)) +
    geom_line(aes(colour = Pepisode)) +
    scale_colour_gradient(low = "black", high = "red") +
    facet_grid(Condition ~ ., space = "free") + 
    geom_vline(xintercept = x, colour = "blue", size = 2) +
    theme_few() +
    theme(text = element_text(size = 18),
          panel.margin = unit(0.5, "lines"),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey", colour = "black"),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(size = 18, face = "plain")) +
    labs(x = "Time (ms)",
         y = expression(paste('Voltage (', mu, 'V)'))) +
    guides(colour = FALSE)
  pLine
  
  ## code below is so we can move the facet labels from the right to the top
  # Get the gtable
  gt <- ggplotGrob(pLine)
  
  # Get the position of the panels in the layout
  panels <-c(subset(gt$layout, name=="panel", se=t:r))
  
  # Add a row above each panel
  for(i in rev(panels$t-1)) gt = gtable_add_rows(gt, unit(.5, "lines"), i)
  
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
    gt$heights[i-1] = list(unit(1, "lines"))
    gt$heights[i-2] = list(unit(0.2, "lines"))
  }
  
  # Draw it
  grid.newpage()
  grid.draw(gt)
  ## end raw trace plot
  
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
    theme(text = element_text(size = 18),
          strip.text = element_blank(),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(size = 18, face = "plain")) +
    scale_y_continuous(limits = c(0,1))
    
  # save combined figure
  electrodeID <- unique(thisData$ElectrodeID)
  freq <- unique(thisData$Frequency)
  trial <- unique(thisData$RealTrialNumber)
  fileName <- paste0('Figures/RawTracePlusPepisode/', electrodeID, '_', freq, 'Hz_Trial', trial, '.pdf')
  pdf(file = fileName, width = 16, height = 4)
  grid.arrange(gt, pBar, ncol = 2, widths = c(5, 1.5))
  dev.off()
}


