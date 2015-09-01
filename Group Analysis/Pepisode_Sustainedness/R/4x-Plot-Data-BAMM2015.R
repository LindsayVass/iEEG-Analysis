# Date:     19 August 2015
# Purpose:  This script will plot data for the BAMM 2015 talk. It will plot
#           the number of significant electrodes showing a longer post-event
#           duration for navigation relative to teleportation, specifically for
#           Delta-Theta and teleporter entry. It will also plot the difference
#           between navigation/teleportation for each electrode.

library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

load('Rda/allAnalyzedData.Rda')
load('Rda/BAMM2015.Rda')
#dir.create('Figures/BAMM2015/')

# Functions ---------------------------------------------------------------

facetLabeller <- function(variable, value) {
  label <- facetLabels[value]
}

# Make plot of individual electrode effects -------------------------------

validData$TimeType <- factor(validData$TimeType, c('NT', 'FT'))
facetLabels <- c('\nShort Teleport Time\n', '\nLong Teleport Time\n')

p <- validData %>%
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
        panel.border = element_rect(colour = "black", size = 0.75, fill = NA),
        panel.background = element_rect(fill = "black"),
        panel.grid.major.y = element_line(colour = "dimgray")) +
  ylab("Mean Duration (ms)") +
  facet_grid(~ TimeType, labeller = facetLabeller)
p
ggsave('Figures/BAMM2015/PostEntryEpisodeDuration_DeltaThetaOnly_LinePlot.png')

pBlank <- validData %>%
  ggplot(aes(x = Condition, 
             y = MeanDuration, 
             ymin = MeanDuration - SEM, 
             ymax = MeanDuration + SEM,
             group = ElectrodeID,
             colour = Difference)) +
  geom_blank() +
  scale_color_gradientn(colours = rainbow(4)) +
  theme_stata() +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 30),
        legend.text = element_blank(),
        legend.title = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1.5),
        strip.background = element_rect(colour = "black", size = 0.75),
        panel.border = element_rect(colour = "black", size = 0.75, fill = NA),
        panel.background = element_rect(fill = "black"),
        panel.grid.major.y = element_line(colour = "dimgray")) +
  ylab("Mean Duration (ms)") +
  facet_grid(~ TimeType, labeller = facetLabeller)
pBlank
ggsave('Figures/BAMM2015/PostEntryEpisodeDuration_DeltaThetaOnly_LinePlot_BLANK.png')

# Make plot of electrode counts -------------------------------------------

postEpisodeElectrodeCounts <- meanPostEpisodeDuration %>%
  select(FrequencyBand, TimeType, BoundaryType) %>%
  unique() %>%
  anti_join(postEpisodeDurSigElectrodes) %>%
  mutate(NSigElectrodes = 0, P = 1) %>%
  rbind(postEpisodeDurSigElectrodes) %>%
  filter(FrequencyBand == "Delta-Theta",
         BoundaryType == "Entry")

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
