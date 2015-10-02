
library(dplyr)
library(reshape2)
library(ggplot2)

load('Rda/allCleanData.Rda')

deltaThetaPepisode <- allPepisode %>%
  filter(FrequencyBand == "Delta-Theta",
         TimePoint == "Tele") %>%
  melt(id.vars = c('ElectrodeID', 'RealTrialNumber', 'TrialSpaceType', 'TrialTimeType', 'Frequency'), measure.vars = c('TelePepisode', 'NavPepisode'), variable.name = "Condition", value.name = "Pepisode")
deltaThetaPepisode$Condition <- plyr::mapvalues(deltaThetaPepisode$Condition, from = c("TelePepisode", "NavPepisode"), to = c("Teleportation", "Navigation"))

binaryPepisode <- deltaThetaPepisode %>%
  mutate(PepisodeBin = ifelse(Pepisode > 0, 1, 0)) %>%
  group_by(Condition, Frequency) %>%
  summarise(TotalPepisodeBin = sum(PepisodeBin))

hist <- ggplot(binaryPepisode, aes(x = Frequency, y = TotalPepisodeBin, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge")
  
groupKsResult <- ks.test(binaryPepisode$TotalPepisodeBin[which(binaryPepisode$Condition == "Teleportation")], binaryPepisode$TotalPepisodeBin[which(binaryPepisode$Condition == "Navigation")])

electrodeKsResults <- vector(mode = "list", length = nlevels(deltaThetaPepisode$ElectrodeID))
for (i in 1:nlevels(deltaThetaPepisode$ElectrodeID)) {
  thisData <- deltaThetaPepisode %>%
    filter(ElectrodeID == levels(ElectrodeID)[i]) %>%
    mutate(PepisodeBin = ifelse(Pepisode > 0, 1, 0)) %>%
    group_by(Condition, Frequency) %>%
    summarise(TotalPepisodeBin = sum(PepisodeBin))
  thisKsResult <- ks.test(thisData$TotalPepisodeBin[which(thisData$Condition == "Teleportation")], thisData$TotalPepisodeBin[which(thisData$Condition == "Navigation")])
  thisKsResult <- data.frame(ElectrodeID = levels(deltaThetaPepisode$ElectrodeID)[i],
                                        D = thisKsResult$statistic,
                                        PValue = thisKsResult$p.value)
  electrodeKsResults[[i]] <- thisKsResult
}
electrodeKsResults <- data.table::rbindlist(electrodeKsResults)