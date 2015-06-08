# Script Name:  4-Plot-Data.R
# Author:       Lindsay Vass
# Date:         8 June 2015
# Purpose:      This script will plot the data produced by 3-Analyze-Data.R. It
#               will produce line plots that show when pepisodes are present
#               over the course of the trial as well as histograms that show
#               when episodes are most likely to begin and end.

library(ggplot2)
library(dplyr)
library(ggthemes)
library(grid)

# Make line plots ---------------------------------------------------------
timepointMarkers <- data.frame(x = c(0, 1830, 0, 2830), TrialTimeType = c("NT", "NT", "FT", "FT"))
lineData <- episodeData %>%
  group_by(FrequencyBand, TrialTimeType, Time) %>%
  summarise(GroupMean = mean(Mean), GroupSEM = sd(Mean) / sqrt(n())) 
linePlot <- lineData %>%
  ggplot(aes(x = Time, y = GroupMean, ymin = GroupMean - GroupSEM, ymax = GroupMean + GroupSEM)) + 
  geom_ribbon(color = "lightskyblue", fill = "lightskyblue") +
  geom_line(color = "steelblue4") +
  geom_vline(aes(xintercept = x, linetype = "dashed"), timepointMarkers) +
  facet_grid(FrequencyBand ~ TrialTimeType) +
  theme_few() +
  labs(y = "Mean Pepisode", title = "Pepisode Over Time") +
  theme(text = element_text(size = 24),
        axis.text = element_text(size = 18),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 1),
        panel.margin = unit(1, "lines"))
linePlot

dir.create("Figures")
ggsave('Figures/Pepisode_Over_Time.png', linePlot)

# Make histograms ---------------------------------------------------------

onsetData <- onOffData %>%
  ungroup() %>%
  filter(ObservationType == "Onset") %>%
  ggplot(aes(x = Time)) +
  geom_histogram(colour = "steelblue4", fill = "lightskyblue") +
  facet_grid(FrequencyBand ~ TrialTimeType, scales = "free") +
  geom_vline(aes(xintercept = x), timepointMarkers) +
  theme_few() +
  theme(text = element_text(size = 24),
        panel.margin = unit(1, "lines")) +
  labs(y = "# of Episodes", title = "Pepisode Onset Times")
onsetData
ggsave('Figures/Pepisode_Onset_Times_Histogram.png', onsetData)

offsetData <- onOffData %>%
  ungroup() %>%
  filter(ObservationType == "Offset") %>%
  ggplot(aes(x = Time)) +
  geom_histogram(colour = "steelblue4", fill = "lightskyblue") +
  facet_grid(FrequencyBand ~ TrialTimeType, scales = "free") +
  geom_vline(aes(xintercept = x), timepointMarkers) +
  theme_few() +
  theme(text = element_text(size = 24),
        panel.margin = unit(1, "lines")) +
  labs(y = "# of Episodes", title = "Pepisode Offset Times")
offsetData
ggsave('Figures/Pepisode_Offset_Times_Histogram.png', offsetData)