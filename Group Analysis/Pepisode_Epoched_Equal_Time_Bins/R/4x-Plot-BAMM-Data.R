# Script Name:  4x-Plot-BAMM-Data.R
# Author:       Lindsay Vass
# Date:         17 August 2015
# Purpose:      This script will plot the results for the 2015 BAMM talk. This 
#               includes histograms of significant electrodes for delta and 
#               theta only, as well as scatterplots of individual electrodes' 
#               effects.

library(plyr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(reshape2)


dir.create('Figures/BAMM2015/')

load('Rda/allAnalyzedData_SepDeltaTheta.Rda')
load('Rda/allCleanData_SepDeltaTheta.Rda')


# Histogram ---------------------------------------------------------------

allAnalyses <- permResults %>%
  select(FrequencyBand, Contrast) %>%
  unique() %>%
  anti_join(trueNSigElectrodesCorrected) %>%
  mutate(Count = 0,
         CorrP = 1) %>%
  rbind(trueNSigElectrodesCorrected)
allAnalyses$Contrast <- factor(allAnalyses$Contrast, levels = c("Pre < Tele", "Pre > Tele", "Tele < Post", "Tele > Post"))

sigLevel <- 6 # number of electrodes at which results are significant  

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

ggsave('Figures/BAMM2015/Pepisode_by_TimePoint_Histogram.png')

blankHistogram <- allAnalyses %>%
  filter(FrequencyBand == "Delta" | FrequencyBand == "Theta") %>%
  ggplot(aes(x = Contrast, y = Count)) +
  geom_blank() +
  geom_hline(aes(yintercept = sigLevel,
                 color = "red")) +
  facet_grid(. ~ FrequencyBand) +
  ylab("# of Significant Electrodes") +
  theme_stata() +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 18),
        axis.title.y = element_text(vjust = 1),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
        strip.background = element_rect(fill = "grey", colour = "black", size = 1.5)) +
  scale_y_continuous(limits = c(0, nlevels(cleanData$ElectrodeID)))
blankHistogram
ggsave('Figures/BAMM2015/BLANK_Pepisode_by_Timepoint_Histogram_DeltaThetaOnly.png')

deltaThetaHistogram <- allAnalyses %>%
  filter(FrequencyBand == "Delta" | FrequencyBand == "Theta") %>%
  ggplot(aes(x = Contrast, y = Count)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = sigLevel,
                 color = "red")) +
  facet_grid(. ~ FrequencyBand) +
  ylab("# of Significant Electrodes") +
  theme_stata() +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 18),
        axis.title.y = element_text(vjust = 1),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
        strip.background = element_rect(fill = "grey", colour = "black", size = 1.5)) +
  scale_y_continuous(limits = c(0, nlevels(cleanData$ElectrodeID))) 
deltaThetaHistogram

ggsave('Figures/BAMM2015/Pepisode_by_Timepoint_Histogram_DeltaThetaOnly.png')


# Individual electrode data -----------------------------------------------

validData <- cleanData %>%
  filter(FrequencyBand == "Delta" | FrequencyBand == "Theta") %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint) %>%
  summarise(SEM = sd(Pepisode) / sqrt(n()),
            Pepisode = mean(Pepisode))
differenceData <- validData %>%
  dcast(ElectrodeID + FrequencyBand ~ TimePoint, value.var = 'Pepisode') %>%
  mutate(Difference = Tele - ((Pre1 + Post1) / 2)) %>%
  select(ElectrodeID, FrequencyBand, Difference) %>%
  inner_join(validData)
differenceData$TimePoint <- revalue(differenceData$TimePoint, c("Pre1" = "Pre", "Tele" = "Teleport", "Post1" = "Post"))

p <- differenceData %>%
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
        panel.border = element_rect(colour = "black", size = 0.75, fill = NA),
        panel.background = element_rect(fill = "black"),
        panel.grid.major.y = element_line(colour = "dimgray")) +
  ylab(expression("Mean P"["Episode"])) +
  facet_grid(~ FrequencyBand) +
  scale_y_continuous(limits = c(0,1))
p
ggsave('Figures/BAMM2015/Pepisode_by_TimePoint_Scatterplot.png')

blankP <- differenceData %>%
  ggplot(aes(x = TimePoint, 
             y = Pepisode, 
             ymin = Pepisode - SEM, 
             ymax = Pepisode + SEM,
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
  ylab(expression("Mean P"["Episode"])) +
  facet_grid(~ FrequencyBand) +
  scale_y_continuous(limits = c(0,1))
blankP
ggsave('Figures/BAMM2015/Pepisode_by_TimePoint_Scatterplot_BLANK.png')