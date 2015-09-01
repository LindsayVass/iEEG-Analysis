# Script Name:  4x-Plot-BAMM-Data.R
# Author:       Lindsay Vass
# Date:         17 August 2015
# Purpose:      This script will plot the results from 3-Analyze-Data.R. This 
#               includes histograms of significant electrodes counts for each, 
#               frequency band as well as scatterplots of individual electrodes' 
#               effects.

library(plyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(dplyr)

load('Rda/allAnalyzedData.Rda')
load('Rda/allCleanData.Rda')


# Histogram ---------------------------------------------------------------

allAnalyses <- permResults %>%
  as.data.frame() %>%
  select(FrequencyBand, Contrast) %>%
  unique() %>%
  anti_join(trueNSigElectrodesCorrected) %>%
  mutate(Count = 0,
         CorrP = 1) %>%
  rbind(trueNSigElectrodesCorrected)
allAnalyses$Contrast <- factor(allAnalyses$Contrast, levels = c("Pre < Tele", "Pre > Tele", "Tele < Post", "Tele > Post"))

sigLevel <- 7 # number of electrodes at which results are significant  

histogram <- allAnalyses %>%
  ggplot(aes(x = Contrast, y = Count)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = sigLevel, 
                 color = "red")) +
  facet_grid(. ~ FrequencyBand) +
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 24, vjust = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text.x  = element_text(size = 18, angle = 90, hjust = 1)) +
  ylab("# of Significant Electrodes") +
  scale_y_continuous(limits = c(0, nlevels(cleanData$ElectrodeID))) +
  ggtitle("Frequency of Significant Differences in Pepisode by Frequency Band")
histogram

ggsave('Figures/Pepisode_by_TimePoint_Histogram.png')


# Individual electrode data -----------------------------------------------

validData <- cleanData %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint) %>%
  summarise(SEM = sd(Pepisode) / sqrt(n()),
            Pepisode = mean(Pepisode))
differenceData <- validData %>%
  dcast(ElectrodeID + FrequencyBand ~ TimePoint, value.var = 'Pepisode') %>%
  mutate(Difference = Tele - ((Pre1 + Post1) / 2)) %>%
  select(ElectrodeID, FrequencyBand, Difference) %>%
  inner_join(validData)
differenceData$TimePoint <- revalue(differenceData$TimePoint, c("Pre1" = "Pre", "Tele" = "Teleport", "Post1" = "Post"))

for (thisFreq in 1:nlevels(cleanData$FrequencyBand)) {
  p <- differenceData %>%
    filter(FrequencyBand == levels(cleanData$FrequencyBand)[thisFreq]) %>%
    ggplot(aes(x = TimePoint, 
               y = Pepisode, 
               ymin = Pepisode - SEM, 
               ymax = Pepisode + SEM,
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
    ylab(expression("Mean P"["Episode"])) +
    facet_grid(~ FrequencyBand) +
    scale_y_continuous(limits = c(0,1))
  p
  ggsave(paste0('Figures/Pepisode_by_TimePoint_Scatterplot', levels(cleanData$FrequencyBand)[thisFreq], '.png'))
}