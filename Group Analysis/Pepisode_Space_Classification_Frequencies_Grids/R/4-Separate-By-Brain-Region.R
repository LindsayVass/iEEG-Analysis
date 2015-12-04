library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)

load('Rda/allClassificationResults.Rda')


# Functions ---------------------------------------------------------------
getCorrectedP <- function(permData, allBootstrapClassification, electrode, model) {
  bootstrapData <- allBootstrapClassification %>%
    filter(ElectrodeID == electrode,
           Model == model)
  permData <- permData %>%
    filter(ElectrodeID == electrode,
           Model == model)
  trueP <- min(which(bootstrapData$TVal == permData$TVal)) / nrow(bootstrapData)
  thisResult <- data.frame(ElectrodeID = electrode,
                           Model = model,
                           TVal = permData$TVal,
                           CorrP = trueP)
}

sortPerms <- function(permData) {
  rowPos <- data.frame(Row = seq(1, nrow(permData)))
  permData <- permData %>%
    arrange(desc(NSigElectrodes))
  nVals <- data.frame(NSigElectrodes = unique(permData$NSigElectrodes))
  pVals <- data.frame()
  for (thisN in 1:nrow(nVals)) {
    nInd  <- min(which(permData == nVals$NSigElectrodes[thisN]))
    pVal  <- data.frame(P = nInd / nrow(permData))
    pVals <- rbind(pVals, pVal)
  }
  npVals <- cbind(nVals, pVals)
  output <- inner_join(permData, npVals, by = "NSigElectrodes")
}

# Analysis ----------------------------------------------------------------

parietal <- c('LPA3', 'LPA4', 'LPA5', 'LPA6', 'LPA7')
occipital <- c('LPA8', 'LMO6', 'LMO7', 'LMO8', 'LPO5', 'LPO6', 'LPO7', 'LPO8')

allMeanClassificationResults <- separate(allMeanClassificationResults, ElectrodeID, c('Subject', 'Session', 'Electrode'), sep = '_', remove = FALSE)
allBootstrapClassification   <- separate(allBootstrapClassification, ElectrodeID, c('Subject', 'Session', 'Electrode'), sep = '_', remove = FALSE)

parietalResults <- allMeanClassificationResults %>%
  filter(Electrode %in% parietal)
parietalResults$ElectrodeID <- factor(parietalResults$ElectrodeID)
parietalBootstrap <- allBootstrapClassification %>%
  filter(Electrode %in% parietal)
parietalBootstrap$ElectrodeID <- factor(parietalBootstrap$ElectrodeID) 
parietalBootstrap <- as.data.frame(parietalBootstrap) %>%
  group_by(ElectrodeID, Model)

parietalSigClass <- parietalResults %>%
  group_by(Model) %>%
  filter(CorrP < 0.05) %>%
  summarise(NSigElectrodes = n())

occipitalResults <- allMeanClassificationResults %>%
  filter(Electrode %in% occipital)
occipitalResults$ElectrodeID <- factor(occipitalResults$ElectrodeID)
occipitalBootstrap <- allBootstrapClassification %>%
  filter(Electrode %in% occipital)
occipitalBootstrap$ElectrodeID <- factor(occipitalBootstrap$ElectrodeID) 
occipitalBootstrap <- as.data.frame(occipitalBootstrap) %>%
  group_by(ElectrodeID, Model)

occipitalSigClass <- occipitalResults %>%
  group_by(Model) %>%
  filter(CorrP < 0.05) %>%
  summarise(NSigElectrodes = n())

# determine how many electrodes are significant by chance
numIterations <- 1000
parietalMaxElectrodePermutations <- data.frame()
for (thisIteration in 1:numIterations) {
  permData <- parietalBootstrap %>%
    sample_n(1)
  permPVals <- data.frame()
  for (thisElectrode in 1:nlevels(permData$ElectrodeID)) {
    electrode <- levels(permData$ElectrodeID)[thisElectrode]
    for (thisModel in 1:nlevels(permData$Model)) {
      model <- levels(permData$Model)[thisModel]
      thisResult <- getCorrectedP(permData, parietalBootstrap, electrode, model)
      permPVals <- rbind(permPVals, thisResult)
    }
  }
  permPVals <- permPVals %>%
    group_by(Model) %>%
    filter(CorrP < 0.05) %>%
    summarise(NSigElectrodes = n())
  thisMaxElec <- data.frame(NSigElectrodes = max(permPVals$NSigElectrodes))
  parietalMaxElectrodePermutations <- rbind(parietalMaxElectrodePermutations, thisMaxElec)
}

parietalMaxElectrodeP <- sortPerms(parietalMaxElectrodePermutations) %>%
  unique()

occipitalMaxElectrodePermutations <- data.frame()
for (thisIteration in 1:numIterations) {
  permData <- occipitalBootstrap %>%
    sample_n(1)
  permPVals <- data.frame()
  for (thisElectrode in 1:nlevels(permData$ElectrodeID)) {
    electrode <- levels(permData$ElectrodeID)[thisElectrode]
    for (thisModel in 1:nlevels(permData$Model)) {
      model <- levels(permData$Model)[thisModel]
      thisResult <- getCorrectedP(permData, occipitalBootstrap, electrode, model)
      permPVals <- rbind(permPVals, thisResult)
    }
  }
  permPVals <- permPVals %>%
    group_by(Model) %>%
    filter(CorrP < 0.05) %>%
    summarise(NSigElectrodes = n())
  thisMaxElec <- data.frame(NSigElectrodes = max(permPVals$NSigElectrodes))
  occipitalMaxElectrodePermutations <- rbind(occipitalMaxElectrodePermutations, thisMaxElec)
}

occipitalMaxElectrodeP <- sortPerms(occipitalMaxElectrodePermutations) %>%
  unique()


# Plot --------------------------------------------------------------------

occipP <- occipitalSigClass
occipP$Model <- plyr::mapvalues(occipP$Model, c('Both', 'NT', 'FT'), c('All Trials', 'Short Duration', 'Long Duration'))
occipP <- occipP %>%
  mutate(PropSigElectrodes = NSigElectrodes / nlevels(occipitalResults$ElectrodeID)) %>%
  ggplot(aes(x = Model, y = PropSigElectrodes)) +
  geom_bar(stat = 'identity') +
  geom_hline(aes(yintercept = 2 / nlevels(occipitalResults$ElectrodeID)), color = 'red') +
  theme_few() +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank()) +
  ylab('Proportion of Significant Electrodes') +
  ggtitle('Distance Classification by Occipital Electrodes')
ggsave('Figures/Occipital_Classification.pdf', width = 10, height = 8)

parietalP <- parietalSigClass
parietalP$Model <- plyr::mapvalues(parietalP$Model, c('Both', 'NT', 'FT'), c('All Trials', 'Short Duration', 'Long Duration'))
parietalP <- parietalP %>%
  mutate(PropSigElectrodes = NSigElectrodes / nlevels(parietalResults$ElectrodeID)) %>%
  ggplot(aes(x = Model, y = PropSigElectrodes)) +
  geom_bar(stat = 'identity') +
  geom_hline(aes(yintercept = 1 / nlevels(parietalResults$ElectrodeID)), color = 'red') +
  theme_few() +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank()) +
  ylab('Proportion of Significant Electrodes') +
  ggtitle('Distance Classification by Parietal Electrodes')
ggsave('Figures/Parietal_Classification.pdf', width = 10, height = 8)

save(file = 'Rda/allClassificationByRegion.Rda', list = c('allBootstrapClassification', 'allClassificationResults', 'allMeanClassificationResults', 'occipitalBootstrap', 'occipitalMaxElectrodeP', 'occipitalResults', 'occipitalSigClass', 'parietalBootstrap', 'parietalMaxElectrodeP', 'parietalResults', 'parietalSigClass'))