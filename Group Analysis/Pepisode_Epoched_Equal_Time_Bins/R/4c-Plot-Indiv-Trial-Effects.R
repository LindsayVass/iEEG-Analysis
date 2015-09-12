# Script Name:  4c-Plot-Indiv-Trial-Effects.R
# Author:       Lindsay Vass
# Date:         11 September 2015
# Purpose:      This script will plot violin plots for each timepoint for each
#               frequency band for each electrode. 

library(dplyr)
library(ggplot2)

load('Rda/allAnalyzedData.Rda')

#dir.create('Figures/SingleElectrodeData/Violin/')

diffData <- cleanData %>%
  dcast(ElectrodeID + FrequencyBand + TrialNumber ~ TimePoint, value.var = "MeanPepisode") %>%
  mutate(TeleMinusPre = Tele - Pre1,
         TeleMinusPost = Tele - Post1) %>%
  select(-c(Pre1:Post1)) %>%
  melt(id.vars = c("ElectrodeID", "FrequencyBand", "TrialNumber"), variable.name = "Contrast", value.name = "TeleMinus")
diffData$Contrast <- sub("TeleMinus", "", diffData$Contrast)
diffData$Contrast <- factor(diffData$Contrast, levels = c("Pre", "Post"))

plotLabels <- moltenTrueData %>%
  mutate(Type = ifelse(Contrast == "TeleLtPre" | Contrast == "TeleLtPost", "Neg", "Pos"))
plotLabels$Contrast <- sub("TeleLt", "", plotLabels$Contrast)
plotLabels$Contrast <- sub("TeleGt", "", plotLabels$Contrast)
plotLabels <- dcast(plotLabels, ElectrodeID + FrequencyBand + Contrast ~ Type, value.var = "PValue")

for (i in 1:nlevels(diffData$ElectrodeID)) {
  
  pLabel <- plotLabels %>%
    filter(ElectrodeID == levels(ElectrodeID)[i])
  
  p <- diffData %>%
    filter(ElectrodeID == levels(ElectrodeID)[i]) %>%
    ggplot(aes(x = Contrast, y = TeleMinus)) +
    geom_violin() +
    facet_wrap(~FrequencyBand) +
    scale_y_continuous(limits = c(-1, 1.2)) + 
    theme(text = element_text(size = 18),
          axis.title.x = element_blank()) +
    labs(y = "Teleporter Pepisode - Other Pepisode") +
    geom_text(data = pLabel, aes(y = 1.2, label = paste0("Pos P = ", round(Pos, digits = 2)))) +
    geom_text(data = pLabel, aes(y = 1.08, label = paste0("Neg P = ", round(Neg, digits = 2))))
  ggsave(paste0('Figures/SingleElectrodeData/Violin/', levels(diffData$ElectrodeID)[i], '_Pepisode_Timepoint_Difference.png'))
}