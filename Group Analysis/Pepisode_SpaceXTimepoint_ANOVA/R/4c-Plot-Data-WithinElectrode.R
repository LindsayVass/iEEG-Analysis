# Script Name:  4c-Plot-Data-WithinElectrode.R
# Author:       Lindsay Vass
# Date:         29 July 2015
# Purpose:      This script will plot the data from 3c-Analyze-Data-WithinElectrode.R

load('Rda/allAnalyzedDataWithinElectrode.Rda')

library(dplyr)
library(ggplot2)


# Plot histograms ---------------------------------------------------------

missingAnalyses <- allTrueAnovaResults %>%
  select(Contrast, FrequencyBand, TrialTimeType) %>%
  unique() %>%
  anti_join(allSigElectrodes) %>%
  mutate(NSigElectrodes = 0,
         P = 1)
allSigElectrodes <- rbind(allSigElectrodes, missingAnalyses)
allSigElectrodes$Contrast      <- factor(allSigElectrodes$Contrast, levels = c("Space", "Timepoint", "Interaction"))
allSigElectrodes$FrequencyBand <- factor(allSigElectrodes$FrequencyBand, levels = c("Delta-Theta", "Alpha", "Beta", "Gamma"))
allSigElectrodes$TrialTimeType <- factor(allSigElectrodes$TrialTimeType, levels = c("NT", "FT"))

p <- allSigElectrodes %>%
  ggplot(aes(x = Contrast, y = NSigElectrodes)) + 
  geom_bar(stat = "identity") +
  facet_grid(TrialTimeType ~ FrequencyBand) +
  theme(text = element_text(size = 24)) +
  ylab("# Significant Electrodes")
p

ggsave('Figures/SpacexTimepointANOVA_WithinElectrode_Histogram.png')


# Plot scatter plots with condition means for each significant ANOVA --------

for (thisEffect in 1:nrow(sigAnovas)) {
  
  plotTitle <- paste("Electrodes with a significant", sigAnovas$Contrast[thisEffect], "effect for", sigAnovas$FrequencyBand[thisEffect], sigAnovas$TrialTimeType[thisEffect])
  
  thisPlot <- allConditionMeans %>%
    filter(FrequencyBand == sigAnovas$FrequencyBand[thisEffect],
           TrialTimeType == sigAnovas$TrialTimeType[thisEffect]) %>%
    ggplot(aes(x = TimeBin, y = Mean, ymin = Mean - SE, ymax = Mean + SE, color = TrialSpaceType)) +
    geom_point(size = 5) +
    geom_pointrange() +
    facet_wrap(~ElectrodeID, scales = "free_y") +
    ylab("Mean Pepisode") +
    xlab("") +
    theme(text = element_text(size = 24)) +
    ggtitle(plotTitle)
  ggsave(paste0('Figures/SpacexTimepointANOVA_WithinElectrode_Scatter_', sigAnovas$FrequencyBand[thisEffect], '_', sigAnovas$TrialTimeType[thisEffect], '_', sigAnovas$Contrast[thisEffect], '.png'))
}

