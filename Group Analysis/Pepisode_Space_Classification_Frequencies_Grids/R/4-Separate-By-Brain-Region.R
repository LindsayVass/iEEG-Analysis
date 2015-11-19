library(tidyr)

load('Rda/allClassificationResults.Rda')

parietal <- c('LPA3', 'LPA4', 'LPA5', 'LPA6', 'LPA7')
occipital <- c('LPA8', 'LMO6', 'LMO7', 'LMO8', 'LPO5', 'LPO6', 'LPO7', 'LPO8')

allMeanClassificationResults <- separate(allMeanClassificationResults, ElectrodeID, c('Subject', 'Session', 'Electrode'), sep = '_', remove = FALSE)

parietalResults <- allMeanClassificationResults %>%
  filter(Electrode %in% parietal)

parietalSigClass <- parietalResults %>%
  group_by(Model) %>%
  filter(CorrP < 0.05) %>%
  summarise(NSigElectrodes = n())

occipitalResults <- allMeanClassificationResults %>%
  filter(Electrode %in% occipital)

occipitalSigClass <- occipitalResults %>%
  group_by(Model) %>%
  filter(CorrP < 0.05) %>%
  summarise(NSigElectrodes = n())