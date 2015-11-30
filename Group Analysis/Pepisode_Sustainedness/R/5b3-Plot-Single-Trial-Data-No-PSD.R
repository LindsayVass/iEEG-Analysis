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

#source('../../functions/ggplot2_formatter.R')

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

matData           <- readMat('mat/allSingleTrialsRawData.mat')
allNavEEGData     <- matData$allNavEEGData
allTeleEEGData    <- matData$allTeleEEGData
allNavBinaryData  <- matData$allNavBinaryData
allTeleBinaryData <- matData$allTeleBinaryData
allNavEEGTimes    <- matData$allNavEEGTimes
allTeleEEGTimes   <- matData$allTeleEEGTimes

load('Rda/allSingleTrials.Rda')
load('Rda/allCleanData.Rda')

#dir.create('Figures/RawTracePlusPepisodeNoPSD/')

# Filter power and pepisode data ------------------------------------------

plotSubs <- allTrials %>%
  select(ElectrodeID) %>%
  unique()

plotSubsTrials <- allTrials %>%
  select(ElectrodeID, RealTrialNumber)
origPepisode <- allPepisode
allPepisode <- allPepisode %>%
  inner_join(plotSubsTrials) %>%
  filter(FrequencyBand == "Delta-Theta") %>%
  group_by(ElectrodeID, TimePoint, FrequencyBand, RealTrialNumber) %>%
  summarise(MeanTelePepisode = mean(TelePepisode),
            MeanNavPepisode = mean(NavPepisode))

# Clean up raw data -------------------------------------------------------

allNavData  <- vector(mode = "list", length = nrow(allTrials))
allTeleData <- vector(mode = "list", length = nrow(allTrials))
for (i in 1:nrow(allTrials)) {
  allNavData[[i]] <- data.frame(ElectrodeID = allTrials$ElectrodeID[i],
                                RealTrialNumber = allTrials$RealTrialNumber[i],
                                Frequency = allTrials$Frequency[i],
                                TrialSpaceType = allTrials$TrialSpaceType[i],
                                TrialTimeType = allTrials$TrialTimeType[i],
                                TimePoint = allTrials$TimePoint[i],
                                Condition = "Navigation",
                                Time = unlist(allNavEEGTimes[[i]]),
                                EEG = unlist(allNavEEGData[[i]]),
                                Pepisode = unlist(allNavBinaryData[[i]]))
  allTeleData[[i]] <- data.frame(ElectrodeID = allTrials$ElectrodeID[i],
                                 RealTrialNumber = allTrials$RealTrialNumber[i],
                                 Frequency = allTrials$Frequency[i],
                                 TrialSpaceType = allTrials$TrialSpaceType[i],
                                 TrialTimeType = allTrials$TrialTimeType[i],
                                 TimePoint = allTrials$TimePoint[i],
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
fontSize <- 30

for (i in 1:nlevels(allData$ObservationID)) {
  thisData <- allData %>%
    filter(ObservationID == levels(allData$ObservationID)[i])
  thisElecData <- origPepisode %>%
    filter(ElectrodeID == as.character(thisData$ElectrodeID[1]),
           FrequencyBand == 'Delta-Theta') %>%
    group_by(ElectrodeID, TimePoint) %>%
    summarise(MeanTelePepisode = mean(TelePepisode),
              MeanNavPepisode = mean(NavPepisode),
              SEMTelePepisode = sd(TelePepisode) / sqrt(n()),
              SEMNavPepisode = sd(NavPepisode) / sqrt(n()))
  firstMelt <- thisElecData %>%
    melt(id.vars = c("ElectrodeID", "TimePoint"), measure.vars = c("MeanTelePepisode", "MeanNavPepisode"), variable.name = "Condition", value.name = "MeanPepisode")
  secondMelt <- thisElecData %>%
    melt(id.vars = c("ElectrodeID", "TimePoint"), measure.vars = c("SEMTelePepisode", "SEMNavPepisode"), variable.name = "Condition", value.name = "SEMPepisode")
  firstMelt$Condition <- mapvalues(firstMelt$Condition, from = c("MeanNavPepisode", "MeanTelePepisode"), to = c("Navigation", "Teleportation"))
  secondMelt$Condition <- mapvalues(secondMelt$Condition, from = c("SEMNavPepisode", "SEMTelePepisode"), to = c("Navigation", "Teleportation"))
  thisElecData <- inner_join(firstMelt, secondMelt)
  
  thisPepisodeData <- allPepisode %>%
    filter(ElectrodeID == thisData$ElectrodeID[1],
           RealTrialNumber == thisData$RealTrialNumber[1]) %>%
    melt(id.vars = c("ElectrodeID", "TimePoint", "RealTrialNumber"), measure.vars = c("MeanTelePepisode", "MeanNavPepisode"), variable.name = "Condition", value.name = "MeanPepisode")
  thisPepisodeData$Condition <- mapvalues(thisPepisodeData$Condition, from = c("MeanNavPepisode", "MeanTelePepisode"), to = c("Navigation", "Teleportation"))
  thisPepisodeData$Condition <- factor(thisPepisodeData$Condition, levels = c("Navigation", "Teleportation"))
  thisPepisodeData$TimePoint <- mapvalues(thisPepisodeData$TimePoint, from = c("Pre1", "Tele", "Post1"), to = c("Pre", "Tele", "Post"))
  
  thisElecData$Condition <- factor(thisElecData$Condition, levels = c("Navigation", "Teleportation"))
  thisElecData$TimePoint <- mapvalues(thisElecData$TimePoint, from = c("Pre1", "Tele", "Post1"), to = c("Pre", "Tele", "Post"))
  
  
  if (thisData$TrialTimeType[1] == "NT") {
    x = c(0, 1830)
  }  else {
    x = c(0, 2830)
  }
  
  # make plots of raw data
  pLine <- thisData %>%
    ggplot(aes(x = Time, y = EEG)) +
    geom_line() +
    facet_grid(Condition ~ ., space = "free") + 
    geom_vline(xintercept = x, colour = "blue", size = 2, linetype = "longdash") +
    theme_few() +
    theme(text = element_text(size = fontSize),
          panel.margin = unit(0.5, "lines"),
          strip.text.y = element_text(angle = 0),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x = element_text(face = "plain"),
          axis.title.y = element_text(vjust = 1),
          strip.background = element_rect(colour = "black", fill = "deepskyblue"),
          plot.margin = unit(c(0.3, 1, 0.3, 0.3), "cm")) +
    labs(x = "Time (ms)",
         y = expression(paste('Voltage (', mu, 'V)'))) +
    guides(colour = FALSE)
  pLine
  
  gLine <- moveStripToTop(pLine)
  grid.newpage()
  grid.draw(gLine)
  
  # plot bar graph showing mean pepisode at each timepoint for ALL TRIALS
  pBar2 <- thisElecData %>%
    ggplot(aes(x = TimePoint, y = MeanPepisode, ymin = MeanPepisode - SEMPepisode, ymax = MeanPepisode + SEMPepisode, fill = TimePoint)) +
    geom_bar(stat = "identity") +
    geom_errorbar() +
    facet_grid(Condition ~ .) +
    labs(y = expression(paste("All Trials Delta-Theta ", P[Episode])),
         x = "  ") +
    theme_few() +
    theme(text = element_text(size = fontSize),
          panel.margin = unit(0.5, "lines"),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "deepskyblue", colour = "black"),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          #  axis.title.x = element_blank(),
          axis.title.x = element_text(face = "plain"),
          axis.text.x = element_text(colour = "white"),
          axis.ticks.x = element_line(colour = "white"),
          axis.title.y = element_text(vjust = 1),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    guides(fill = FALSE) 
  
  gBar2 <- moveStripToTop(pBar2)
  grid.newpage()
  grid.draw(gBar2)
  
  tmp <- ggplot_build(pBar2)
  tmp <- tmp$panel$ranges
  ylim2 <- tmp[[1]]$y.range
  ylabels2 <- tmp[[1]]$y.labels
  ybreaks2 <- as.numeric(ylabels2)
  
  # plot bar graph showing mean pepisode at each timepoint
  pBar <- thisPepisodeData %>%
    ggplot(aes(x = TimePoint, y = MeanPepisode, fill = TimePoint)) +
    geom_bar(stat = "identity") +
    facet_grid(Condition ~ .) +
    labs(y = expression(paste("Mean Delta-Theta ", P[Episode])),
         x = "  ") +
    theme_few() +
    theme(text = element_text(size = fontSize),
          panel.margin = unit(0.5, "lines"),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "deepskyblue", colour = "black"),
          line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          #  axis.title.x = element_blank(),
          axis.title.x = element_text(face = "plain"),
          axis.text.x = element_text(colour = "white"),
          axis.ticks.x = element_line(colour = "white"),
          axis.title.y = element_text(vjust = 1),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    guides(fill = FALSE) +
    scale_y_continuous(lim = c(0, ylim2[2]), breaks = ybreaks2)
  
  gBar <- moveStripToTop(pBar)
  grid.newpage()
  grid.draw(gBar)
  
  tmp <- ggplot_build(pBar)
  tmp <- tmp$panel$ranges
  ylim <- tmp[[1]]$y.range
  ylabels <- tmp[[1]]$y.labels
  ybreaks <- as.numeric(ylabels)
  
  
  
  # save combined figure
  electrodeID <- unique(thisData$ElectrodeID)
  freq <- unique(thisData$Frequency)
  trial <- unique(thisData$RealTrialNumber)
  fileName <- paste0('Figures/RawTracePlusPepisodeNoPSD/', electrodeID, '_', freq, 'Hz_Trial', trial, '.pdf')
  pdf(file = fileName, width = 16, height = 6)
  #grid.arrange(gt, pBar, ncol = 2, widths = c(6, 1.5))
  grid.arrange(gLine, gBar, gBar2, ncol = 3, widths = c(50, 25, 25), padding = 10)
  dev.off()
}


