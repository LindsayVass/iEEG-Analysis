# Script Name:  3-Plot-Individual-Electrodes.R
# Author:       Lindsay Vass
# Date:         14 August 2015
# Purpose:      This script will make plots of data from individual electrodes 
#               showing how log(power) varies across conditions and frequencies
#               according to landmark richness and speed. For landmark richness
#               it will separately plot power for rich and poor landmarks. For
#               speed, it will separately plot power for higher and lower speed,
#               divided by a median split. This will produce plots similar to 
#               Figure 2 in Watrous et al (2011) J Neurophys.


library(dplyr)
library(ggplot2)



# Functions ---------------------------------------------------------------

plotSpeed <- function(inputData, frequencyThresh) {
  labels <- inputData$Frequency %>%
    unique()
  labels <- labels[labels < frequencyThresh]
  
  p <- inputData %>%
    filter(Frequency < frequencyThresh) %>%
    ggplot(aes(x = Frequency, 
               y = MeanPower, 
               ymin = MeanPower - SEMPower, 
               ymax = MeanPower + SEMPower, 
               color = MedSpeedSplit)) +
    geom_line(size = 1) +
    geom_linerange(size = 1) +
    scale_x_log10(breaks = labels, labels = signif(labels, digits = 3)) +
    ylab('Mean log(Power)') +
    xlab('Frequency (Hz)') +
    theme(text = element_text(size = 24)) +
    guides(color = guide_legend(title = "Speed"))
  ggsave(filename = paste0('Figures/Speed/', inputData$ElectrodeID, '_LinePlot_Power_by_Frequency_and_Speed.png'))
  return(p)
}

plotLandmarks <- function(inputData, frequencyThresh) {
  labels <- inputData$Frequency %>%
    unique()
  labels <- labels[labels < frequencyThresh]
  
  p <- inputData %>%
    filter(Frequency < frequencyThresh) %>%
    ggplot(aes(x = Frequency, 
               y = MeanPower, 
               ymin = MeanPower - SEMPower, 
               ymax = MeanPower + SEMPower, 
               color = Landmarks)) +
    geom_line(size = 1) +
    geom_linerange(size = 1) +
    scale_x_log10(breaks = labels, labels = signif(labels, digits = 3)) +
    ylab('Mean log(Power)') +
    xlab('Frequency (Hz)') +
    theme(text = element_text(size = 24))
  ggsave(filename = paste0('Figures/Landmarks/', inputData$ElectrodeID, '_LinePlot_Power_by_Frequency_and_Landmarks.png'))
  return(p)
}

# prepare data for plotting -----------------------------------------------

# select electrodes with significant results < 10 Hz
noCentralStatsCorrected$Frequency <- as.numeric(noCentralStatsCorrected$Frequency)
noCentralSigStats <- noCentralStatsCorrected %>%
  filter(CorrP < 0.05,
         Frequency < 10) %>%
  ungroup() %>%
  group_by(term) %>%
  arrange(CorrP)
noCentralSigSpeed <- noCentralSigStats %>%
  filter(term == "Speed") %>%
  ungroup() %>%
  select(ElectrodeID) %>%
  unique()
noCentralSigLandmarks <- noCentralSigStats %>%
  filter(term == "LandmarksRICH") %>%
  ungroup() %>%
  select(ElectrodeID) %>%
  unique()

# filter and summarize raw data
noCentralRawData <- allRawData %>%
  filter(Landmarks != "CENTRAL") %>%
  select(Frequency:ElectrodeID) %>%
  group_by(ElectrodeID, Frequency)
medianSpeed <- noCentralRawData %>%
  ungroup() %>%
  group_by(ElectrodeID) %>%
  summarise(MedianSpeed = median(Speed))
noCentralRawData <- noCentralRawData %>%
  inner_join(medianSpeed) %>%
  mutate(MedSpeedSplit = ifelse(Speed < MedianSpeed, "Slow", "Fast"))
powerByMedSplitSpeed <- noCentralRawData %>%
  group_by(ElectrodeID, Frequency, MedSpeedSplit) %>%
  summarise(MeanPower = mean(Power),
            SEMPower = sd(Power) / sqrt(n())) %>%
  inner_join(noCentralSigSpeed)
powerByLandmarks <- noCentralRawData %>%
  group_by(ElectrodeID, Frequency, Landmarks) %>%
  summarise(MeanPower = mean(Power),
            SEMPower = sd(Power) / sqrt(n())) %>%
  inner_join(noCentralSigLandmarks)
powerByLandmarks$Landmarks <- factor(powerByLandmarks$Landmarks, levels = rev(levels(powerByLandmarks$Landmarks)))

# plot power by frequency line plots --------------------------------------

dir.create('Figures/Speed')
dir.create('Figures/Landmarks')

freqThresh <- 10

powerByMedSplitSpeed %>%
  ungroup() %>%
  group_by(ElectrodeID) %>%
  do(plotSpeed(., freqThresh))

powerByLandmarks %>%
  ungroup() %>%
  group_by(ElectrodeID) %>%
  do(plotLandmarks(., freqThresh))



