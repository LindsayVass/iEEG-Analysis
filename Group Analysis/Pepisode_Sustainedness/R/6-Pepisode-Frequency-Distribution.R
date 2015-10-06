
library(dplyr)
library(reshape2)
library(ggplot2)
library(coin)
library(permute)
library(car)

load('Rda/allCleanData.Rda')
#dir.create('Figures/FrequencyDist/')

deltaThetaPepisode <- allPepisode %>%
  filter(FrequencyBand == "Delta-Theta",
         TimePoint == "Tele") %>%
  melt(id.vars = c('ElectrodeID', 'RealTrialNumber', 'TrialSpaceType', 'TrialTimeType', 'Frequency'), measure.vars = c('TelePepisode', 'NavPepisode'), variable.name = "Condition", value.name = "Pepisode")
deltaThetaPepisode$Condition <- plyr::mapvalues(deltaThetaPepisode$Condition, from = c("TelePepisode", "NavPepisode"), to = c("Teleportation", "Navigation"))

binaryPepisode <- deltaThetaPepisode %>%
  mutate(PepisodeBin = ifelse(Pepisode > 0, 1, 0)) %>%
  group_by(Condition, Frequency) %>%
  summarise(TotalPepisodeBin = sum(PepisodeBin))

meanPepisode <- deltaThetaPepisode %>%
  group_by(Condition, Frequency) %>%
  summarise(MeanPepisode = mean(Pepisode),
            SEMPepisode = sd(Pepisode)/sqrt(n()))

binHist <- ggplot(binaryPepisode, aes(x = Frequency, y = TotalPepisodeBin, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge")
meanHist <- ggplot(meanPepisode, aes(x = Frequency, y = MeanPepisode, ymin = MeanPepisode - SEMPepisode, ymax = MeanPepisode + SEMPepisode, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(position = "dodge")

ggsave(filename = 'Figures/FrequencyDist/Mean_Pepisode_Distribution.png', plot = meanHist)
  
binGroupKsResult <- ks.test(binaryPepisode$TotalPepisodeBin[which(binaryPepisode$Condition == "Teleportation")], binaryPepisode$TotalPepisodeBin[which(binaryPepisode$Condition == "Navigation")])
meanGroupKsResult <- ks.test(meanPepisode$MeanPepisode[which(meanPepisode$Condition == "Teleportation")], meanPepisode$MeanPepisode[which(meanPepisode$Condition == "Navigation")])

wilcox_results <- vector(mode = "list", length = length(unique(deltaThetaPepisode$Frequency)))
for (i in 1:length(unique(deltaThetaPepisode$Frequency))) {
  thisData <- deltaThetaPepisode %>%
    filter(Frequency == unique(deltaThetaPepisode$Frequency)[i]) %>%
    group_by(ElectrodeID, Condition) %>%
    summarise(MeanPepisode = mean(Pepisode)) 
  thisTest <- wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = thisData, distribution = approximate(B = 10000))
  thisResult <- data.frame(Frequency = unique(deltaThetaPepisode$Frequency)[i],
                           Statistic = statistic(thisTest)[[1]],
                           PValue = pvalue(thisTest)[[1]])
  wilcox_results[[i]] <- thisResult
}
wilcox_results <- data.table::rbindlist(wilcox_results)

# run 2-way nonparametric anova
wideData <- deltaThetaPepisode %>%
  group_by(ElectrodeID, Frequency, Condition) %>%
  summarise(MeanPepisode = mean(Pepisode)) %>%
  dcast(ElectrodeID ~ Frequency + Condition, value.var = "MeanPepisode")

# Convert to matrix form and remove ElectrodeID
wideMatrix <- data.matrix(wideData[2:dim(wideData)[2]])

# Use lm() to generate a linear model on all columns
firstLM <- lm(wideMatrix ~ 1)

# iData for two-way ANOVA
freqFactor <- deltaThetaPepisode %>%
  select(Frequency) %>%
  unique()
freqFactor$Frequency <- as.factor(freqFactor$Frequency)

condition  <- c(rep(c("Teleportation", "Navigation"), nlevels(freqFactor$Frequency)))
frequency  <- rep(levels(freqFactor$Frequency), each = 2)
iData <- data.frame(condition, frequency)

twoWayAnova <- Anova(firstLM, 
                     idata = iData, 
                     idesign = ~condition * frequency)
twoWayAnova <- summary(twoWayAnova)$univariate.tests