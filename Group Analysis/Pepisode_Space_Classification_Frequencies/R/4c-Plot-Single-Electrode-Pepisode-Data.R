# Script Name:  4c-Plot-Single-Electrode-Pepisode-Data.R
# Author:       Lindsay Vass
# Date:         11 September
# Purpose:      This script will make a scatterplot showing the difference between
#               navigation/teleportation pepisode at each timepoint, with 
#               separate plots for each electrode. This will allow us to visualize
#               the differences on a trialwise basis. 

library(dplyr)
library(ggplot2)

load('Rda/allCleanData.Rda')
load('Rda/allAnalyzedData.Rda')
remove(navSustain, teleSustain, sessionInfo, durationPermResults, durationTrueDataCorrected, durationTrueNSigElectrodesCorrected, pepisodePermResults, permData, permDurationData)

#dir.create('Figures/SingleElectrodePepisode/')

# separate plot for each frequency band
for (thisElec in 1:nlevels(allPepisode$ElectrodeID)) {
  
  thisData <- allPepisode %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElec]) %>%
    ungroup() %>%
    group_by(TimePoint, FrequencyBand) %>%
    mutate(NavMinusTele = NavPepisode - TelePepisode)
  
  pLabels <- pepisodeTrueDataCorrected %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElec]) %>%
    ungroup() %>%
    group_by(TimePoint, FrequencyBand)
  
  p <- thisData %>%
    ggplot(aes(x = TimePoint, y = NavMinusTele)) +
    geom_violin() +
    facet_wrap(~FrequencyBand) +
    geom_text(data = pLabels, aes(label = paste0("P = ", round(pLabels$CorrP, digits = 2)), x = TimePoint, y = 1.2)) +
    labs(y = 'Navigation Pepisode - Teleporter Pepisode')
  
  ggsave(paste0('Figures/SingleElectrodePepisode/', levels(allPepisode$ElectrodeID)[thisElec], '_pepisode_by_timepoint_violin.png'))
  
}

# separate plot for each frequency within delta-theta
dtData <- allPepisode %>%
  filter(FrequencyBand == "Delta-Theta")
for (thisElec in 1:nlevels(allPepisode$ElectrodeID)) {
  
  thisData <- dtData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElec]) %>%
    ungroup() %>%
    group_by(TimePoint, Frequency) %>%
    mutate(NavMinusTele = NavPepisode - TelePepisode)
  
  # violin plot doesn't like constant values (e.g., all 0) so add a tiny jitter
  zeroSD <- thisData %>%
    summarise(SD = sd(NavMinusTele)) %>%
    filter(SD == 0) %>%
    select(-SD)
  for (i in 1:nrow(zeroSD)) {
    ind <- min(which(thisData$TimePoint == zeroSD$TimePoint[i] & thisData$Frequency == zeroSD$Frequency[i]))
    thisData$NavMinusTele[ind] <- thisData$NavMinusTele[ind] + 0.00001
  }
  
  pLabels <- pepisodeTrueDataCorrected %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElec],
           FrequencyBand == "Delta-Theta") %>%
    ungroup() %>%
    group_by(TimePoint)
  
  freqs <- unique(thisData$Frequency) %>%
    rep(each = nrow(pLabels))
  
  pLabels <- cbind(pLabels, freqs) %>%
    rename(Frequency = freqs)
  
  p <- thisData %>%
    ggplot(aes(x = TimePoint, y = NavMinusTele)) +
    geom_violin() +
    facet_wrap(~Frequency) +
    geom_text(data = pLabels, aes(label = paste0("P = ", round(pLabels$CorrP, digits = 2)), x = TimePoint, y = 1.2)) +
    labs(y = 'Navigation Pepisode - Teleporter Pepisode')
  
  ggsave(paste0('Figures/SingleElectrodePepisode/', levels(allPepisode$ElectrodeID)[thisElec], '_delta-theta_pepisode_by_timepoint_violin.png'))
  
}