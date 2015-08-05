# Script Name:  4-Plot-Data.R
# Author:       Lindsay Vass
# Date:         23 July 2015
# Purpose:      This script will plot the data from 3-Analyze-Data.R

library(dplyr)
library(ggplot2)

load('Rda/allAnalyzedData.Rda')


# Plot pepisode histogram -------------------------------------------------

pepisodeElectrodeCounts <- pepisodeResults %>%
  select(FrequencyBand, TimePoint) %>%
  unique() %>%
  anti_join(pepisodeSigElectrodes) %>%
  mutate(NSigElectrodes = 0, P = 1) %>%
  rbind(pepisodeSigElectrodes) 

freqBandOrder <- c("Delta-Theta", "Alpha", "Beta", "Gamma")
pepisodeElectrodeCounts$TimePoint <- sub('Pre1', 'Pre', pepisodeElectrodeCounts$TimePoint)
pepisodeElectrodeCounts$TimePoint <- sub('Post1', 'Post', pepisodeElectrodeCounts$TimePoint)
pepisodeElectrodeCounts$FrequencyBand <- factor(pepisodeElectrodeCounts$FrequencyBand, levels = freqBandOrder)
pepisodeElectrodeCounts$TimePoint <- factor(pepisodeElectrodeCounts$TimePoint, levels = c("Pre", "Tele", "Post"))

pepisodeHist <- pepisodeElectrodeCounts %>%
  ggplot(aes(x = TimePoint, y = NSigElectrodes)) +
  geom_histogram(stat = "identity") +
  facet_grid(~FrequencyBand) +
  ggtitle("Number of Electrodes with Lower Pepisode During Teleportation Relative to Navigation") +
  ylab('# of Significant Electrodes')
pepisodeHist
ggsave(filename = 'Figures/Histogram_Pepisode_Tele_lt_Nav.png')

# Plot episode duration histogram -----------------------------------------

episodeElectrodeCounts <- meanEpisodeDuration %>%
  select(FrequencyBand, TimeType, BoundaryType) %>%
  unique() %>%
  anti_join(episodeDurSigElectrodes) %>%
  mutate(NSigElectrodes = 0, P = 1) %>%
  rbind(episodeDurSigElectrodes)

timeTypeOrder <- c("NT", "FT")
boundaryTypeOrder <- c("Entry", "Exit")
episodeElectrodeCounts$TimeType <- factor(episodeElectrodeCounts$TimeType, levels = timeTypeOrder)
episodeElectrodeCounts$FrequencyBand <- factor(episodeElectrodeCounts$FrequencyBand, levels = freqBandOrder)
episodeElectrodeCounts$BoundaryType <- factor(episodeElectrodeCounts$BoundaryType, levels = boundaryTypeOrder)

episodeDurHist <- episodeElectrodeCounts %>%
  ggplot(aes(x = BoundaryType, y = NSigElectrodes)) +
  geom_histogram(stat = "identity") +
  facet_grid(TimeType ~ FrequencyBand) +
  ggtitle("Number of Electrodes Whose Oscillations are Shorter During Teleportation Relative to Navigation") +
  ylab("# of Significant Electrodes")
episodeDurHist
ggsave(filename = 'Figures/Histogram_Episode_Duration_Tele_lt_Nav.png')

# Plot post-event episode duration histogram ------------------------------

postEpisodeElectrodeCounts <- meanPostEpisodeDuration %>%
  select(FrequencyBand, TimeType, BoundaryType) %>%
  unique() %>%
  anti_join(postEpisodeDurSigElectrodes) %>%
  mutate(NSigElectrodes = 0, P = 1) %>%
  rbind(postEpisodeDurSigElectrodes)

postEpisodeElectrodeCounts$TimeType <- factor(postEpisodeElectrodeCounts$TimeType, levels = timeTypeOrder)
postEpisodeElectrodeCounts$FrequencyBand <- factor(postEpisodeElectrodeCounts$FrequencyBand, levels = freqBandOrder)
postEpisodeElectrodeCounts$BoundaryType <- factor(postEpisodeElectrodeCounts$BoundaryType, levels = boundaryTypeOrder)

postEpisodeDurHist <- postEpisodeElectrodeCounts %>%
  ggplot(aes(x = BoundaryType, y = NSigElectrodes)) +
  geom_histogram(stat = "identity") +
  facet_grid(TimeType ~ FrequencyBand) +
  ggtitle("Number of Electrodes Whose Post-Event Oscillations (Teleporter Entry/Exit)\nare Shorter During Teleportation Relative to Navigation") +
  ylab("# of Significant Electrodes")
postEpisodeDurHist
ggsave(filename = 'Figures/Histogram_Post_Event_Episode_Duration_Tele_lt_Nav.png')
