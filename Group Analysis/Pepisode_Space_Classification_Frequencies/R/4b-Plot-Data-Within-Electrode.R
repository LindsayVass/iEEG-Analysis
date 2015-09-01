# Script Name:  4b-Plot-Data-Within-Electrode.R
# Author:       Lindsay Vass
# Date:         19 August 2015
# Purpose:      This script will plot the classification results from 
#               3b-Analyze-Data-Within-Electrode.R

library(dplyr)
library(ggplot2)
library(ggthemes)
load('Rda/allClassificationResults.Rda')

theData <- allClassificationResults %>%
  group_by(ElectrodeID, Model) %>%
  summarise(SEM = sd(Accuracy) / sqrt(n()), Accuracy = mean(Accuracy))

for (i in 1:nlevels(theData$Model)) {
  
  p <- theData %>%
    filter(Model == levels(Model)[i]) %>%
    transform(ElectrodeID = reorder(ElectrodeID, -Accuracy)) %>%
    ggplot(aes(x = ElectrodeID, y = Accuracy, ymin = Accuracy - SEM, ymax = Accuracy + SEM)) +
    geom_bar(stat = "identity") +
    geom_errorbar() +
    theme_stata() +
    scale_y_continuous(limits = c(0, 1)) +
    theme(text = element_text(size = 24),
          axis.text.x = element_blank(),
          axis.title.y = element_text(vjust = 1.5),
          axis.ticks.x = element_blank(),
          plot.background = element_rect(fill = "white")) +
    geom_hline(yintercept = 0.5, color = "red") +
    xlab("Electrode")
  p
  ggsave(paste0('Figures/SpaceClassification_WithinElectrode_BarChart_', levels(theData$Model)[i], '.png'))
}

p <- theData %>%
  filter(Model == levels(Model)[i]) %>%
  transform(ElectrodeID = reorder(ElectrodeID, -Accuracy)) %>%
  ggplot(aes(x = ElectrodeID, y = Accuracy, ymin = Accuracy - SEM, ymax = Accuracy + SEM)) +
  geom_blank() +
  theme_stata() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(vjust = 1.5),
        axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0.5, color = "red") +
  xlab("Electrode")
p
ggsave('Figures/SpaceClassification_WithinElectrode_BarChart_BLANK.png')