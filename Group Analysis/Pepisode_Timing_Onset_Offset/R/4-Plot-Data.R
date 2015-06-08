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

# group data
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
dir.create("Figures")
ggsave('Figures/Pepisode_Over_Time.png', linePlot)

# individual electrode data
for (thisFrequency in 1:nlevels(episodeData$FrequencyBand)) {
  lineNtData <- episodeData %>%
    filter(FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency] &
             TrialTimeType == "NT")%>%
    ggplot(aes(x = Time, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM)) + 
    geom_ribbon(color = "lightskyblue", fill = "lightskyblue") +
    geom_line(color = "steelblue4") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 1830) +
    facet_wrap(~ElectrodeID) +
    theme_few() +
    labs(y = "Mean Pepisode", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Over Time")) +
    theme(text = element_text(size = 24),
          axis.text = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 1),
          panel.margin = unit(1, "lines"))
  
    ggsave(filename = paste0('Figures/', levels(episodeData$FrequencyBand)[thisFrequency], '_NT_Electrodewise_Pepisode_Over_Time.png'),
           plot = lineNtData,
           width = 30,
           height = 15)

  lineFtData <- episodeData %>%
    filter(FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency] &
             TrialTimeType == "FT")%>%
    ggplot(aes(x = Time, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM)) + 
    geom_ribbon(color = "lightskyblue", fill = "lightskyblue") +
    geom_line(color = "steelblue4") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 2830) +
    facet_wrap(~ElectrodeID) +
    theme_few() +
    labs(y = "Mean Pepisode", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Over Time")) +
    theme(text = element_text(size = 24),
          axis.text = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 1),
          panel.margin = unit(1, "lines"))
  
  ggsave(filename = paste0('Figures/', levels(episodeData$FrequencyBand)[thisFrequency], '_FT_Electrodewise_Pepisode_Over_Time.png'),
         plot = lineFtData,
         width = 30,
         height = 15)
}

# Make onset histograms ------------------------------------------------------

# group data
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
ggsave('Figures/Pepisode_Onset_Times_Histogram.png', onsetData)

# electrodewise data
for (thisFrequency in 1:nlevels(episodeData$FrequencyBand)) {
  onsetNtData <- onOffData %>%
    ungroup() %>%
    filter(ObservationType == "Onset" & 
             TrialTimeType == "NT" & 
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) %>%
    ggplot(aes(x = Time)) +
    geom_histogram(colour = "steelblue4", fill = "lightskyblue") +
    facet_wrap(~ElectrodeID, scales = "free_y") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 1830) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Onset Times"))
  ggsave(filename =  paste0('Figures/', levels(episodeData$FrequencyBand)[thisFrequency], '_NT_Electrodewise_Pepisode_Onset_Times_Histogram.png'),
         plot = onsetNtData,
         width = 30,
         height = 15)
  
  onsetFtData <- onOffData %>%
    ungroup() %>%
    filter(ObservationType == "Onset" & 
             TrialTimeType == "FT" & 
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) %>%
    ggplot(aes(x = Time)) +
    geom_histogram(colour = "steelblue4", fill = "lightskyblue") +
    facet_wrap(~ElectrodeID, scales = "free_y") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 2830) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Onset Times"))
  ggsave(filename =  paste0('Figures/', levels(episodeData$FrequencyBand)[thisFrequency], '_FT_Electrodewise_Pepisode_Onset_Times_Histogram.png'),
         plot = onsetFtData,
         width = 30,
         height = 15)
         
}
 

# Make offset histograms ------------------------------------------------------

# group data
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
ggsave('Figures/Pepisode_Offset_Times_Histogram.png', offsetData)

# electrodewise data
for (thisFrequency in 1:nlevels(episodeData$FrequencyBand)) {
  offsetNtData <- onOffData %>%
    ungroup() %>%
    filter(ObservationType == "Offset" & 
             TrialTimeType == "NT" & 
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) %>%
    ggplot(aes(x = Time)) +
    geom_histogram(colour = "steelblue4", fill = "lightskyblue") +
    facet_wrap(~ElectrodeID, scales = "free_y") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 1830) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Offset Times"))
  ggsave(filename =  paste0('Figures/', levels(episodeData$FrequencyBand)[thisFrequency], '_NT_Electrodewise_Pepisode_Offset_Times_Histogram.png'),
         plot = offsetNtData,
         width = 30,
         height = 15)
  
  offsetFtData <- onOffData %>%
    ungroup() %>%
    filter(ObservationType == "Offset" & 
             TrialTimeType == "FT" & 
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) %>%
    ggplot(aes(x = Time)) +
    geom_histogram(colour = "steelblue4", fill = "lightskyblue") +
    facet_wrap(~ElectrodeID, scales = "free_y") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 2830) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Offset Times"))
  ggsave(filename =  paste0('Figures/', levels(episodeData$FrequencyBand)[thisFrequency], '_FT_Electrodewise_Pepisode_Offset_Times_Histogram.png'),
         plot = onsetFtData,
         width = 30,
         height = 15)
  
}