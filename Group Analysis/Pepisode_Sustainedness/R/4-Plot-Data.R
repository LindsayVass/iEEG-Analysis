# Script Name:  4-Plot-Data.R
# Author:       Lindsay Vass
# Date:         23 July 2015
# Purpose:      This script will plot the data from 3-Analyze-Data.R

library(dplyr)
library(ggplot2)
library(ggthemes)

load('Rda/allAnalyzedData.Rda')
load('Rda/postEventDurationData.Rda')

# Functions ---------------------------------------------------------------

facetLabeller <- function(variable, value) {
  label <- facetLabels[value]
}


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

# Make plot of individual electrode effects for post-event episode duration ----

validData$TimeType <- factor(validData$TimeType, c('NT', 'FT'))
facetLabels <- c('\nShort Teleport Time\n', '\nLong Teleport Time\n')

for (i in 1:nlevels(validData$FrequencyBand)) {
  p <- validData %>%
    filter(FrequencyBand == levels(FrequencyBand)[i]) %>%
    ggplot(aes(x = Condition, 
               y = MeanDuration, 
               ymin = MeanDuration - SEM, 
               ymax = MeanDuration + SEM,
               group = ElectrodeID,
               colour = Difference)) +
    geom_point(size = 5) +
    geom_pointrange() +
    geom_line() +
    scale_color_gradientn(colours = rainbow(4)) +
    theme_stata() +
    theme(plot.background = element_rect(fill = "white"),
          text = element_text(size = 30),
          legend.text = element_blank(),
          legend.title = element_text(size = 24),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1.5),
          strip.background = element_rect(colour = "black", size = 0.75),
          panel.border = element_rect(colour = "black", size = 0.75, fill = NA)) +
    ylab("Mean Duration (ms)") +
    facet_grid(~ TimeType, labeller = facetLabeller)
  p
  
  ggsave(paste0('Figures/PostEntryEpisodeDuration_LinePlot_', levels(validData$FrequencyBand)[i], '.png'))
}
