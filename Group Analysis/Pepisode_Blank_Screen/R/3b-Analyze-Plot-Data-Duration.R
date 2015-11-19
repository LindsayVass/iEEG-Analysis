# Script Name:  3-Analyze-Plot-Data.R
# Author:       Lindsay Vass
# Date:         23 October 2015
# Purpose:      This script will test whether Pepisode differs between Blank &
#               Tele and between Blank & Nav using Wilcoxon signed rank tests.
#               It will then plot the mean pepisode across electrodes in a bar chart.

library(dplyr)
library(ggplot2)
library(coin)
library(ggthemes)

load('Rda/allCleanData.Rda')

# get mean value for each electrode/condition
meanPepisode <- allPepisode %>%
  filter(Condition == "Blank" | Condition == "Teleportation_Tele" | Condition == "Navigation_Tele") %>%
  group_by(ElectrodeID, Condition, TrialTimeType) %>%
  summarise(MeanPepisode = mean(Pepisode),
            SEMPepisode = sd(Pepisode) / n())

# exclude electrodes that don't have data for "blank"
goodElec <- meanPepisode %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Teleportation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "Blank") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

goodPepisode <- inner_join(meanPepisode, goodElec)

##### wilcoxon signed rank test for Nav vs Blank (SHORT)
navBlankShort <- meanPepisode %>%
  filter(Condition != "Teleportation_Tele",
         TrialTimeType == 'Short') 

# exclude electrodes that don't have data for "blank"
goodElec <- navBlankShort %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Navigation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "Blank") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

navBlankShort <- inner_join(navBlankShort, goodElec)

navBlankShort$ElectrodeID <- factor(navBlankShort$ElectrodeID)
navBlankShort$Condition <- factor(navBlankShort$Condition)
navBlankShortStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = navBlankShort, distribution = approximate(B = 10000))


##### wilcoxon signed rank test for Nav vs Blank (LONG)
navBlankLong <- meanPepisode %>%
  filter(Condition != "Teleportation_Tele",
         TrialTimeType == 'Long') 

# exclude electrodes that don't have data for "blank"
goodElec <- navBlankLong %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Navigation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "Blank") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

navBlankLong <- inner_join(navBlankLong, goodElec)

navBlankLong$ElectrodeID <- factor(navBlankLong$ElectrodeID)
navBlankLong$Condition <- factor(navBlankLong$Condition)
navBlankLongStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = navBlankLong, distribution = approximate(B = 10000))

##### wilcoxon signed rank test for Tele vs Blank (SHORT)
teleBlankShort <- meanPepisode %>%
  filter(Condition != "Navigation_Tele",
         TrialTimeType == 'Short') 

# exclude electrodes that don't have data for "blank"
goodElec <- teleBlankShort %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Teleportation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "Blank") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

teleBlankShort <- inner_join(teleBlankShort, goodElec)

teleBlankShort$ElectrodeID <- factor(teleBlankShort$ElectrodeID)
teleBlankShort$Condition <- factor(teleBlankShort$Condition)
teleBlankShortStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = teleBlankShort, distribution = approximate(B = 10000))

##### wilcoxon signed rank test for Tele vs Blank (LONG)
teleBlankLong <- meanPepisode %>%
  filter(Condition != "Navigation_Tele",
         TrialTimeType == 'Long') 

# exclude electrodes that don't have data for "blank"
goodElec <- teleBlankLong %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Teleportation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "Blank") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

teleBlankLong <- inner_join(teleBlankLong, goodElec)

teleBlankLong$ElectrodeID <- factor(teleBlankLong$ElectrodeID)
teleBlankLong$Condition <- factor(teleBlankLong$Condition)
teleBlankLongStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = teleBlankLong, distribution = approximate(B = 10000))




# plot the results
goodPepisode$Condition <- plyr::mapvalues(goodPepisode$Condition, from = c("Blank", "Teleportation_Tele", "Navigation_Tele"), to = c("Blank", "Teleportation", "Navigation"))
p <- goodPepisode %>%
  group_by(Condition) %>%
  summarise(GroupMeanPepisode = mean(MeanPepisode),
            GroupSEMPepisode = sd(MeanPepisode) / n()) %>%
  ggplot(aes(x = Condition, y = GroupMeanPepisode, ymin = GroupMeanPepisode - GroupSEMPepisode, ymax = GroupMeanPepisode + GroupSEMPepisode)) +
  geom_bar(stat = "identity") +
  geom_errorbar() +
  theme_few() +
  ylab(expression(paste(Mean~P[Episode]))) +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank())
ggsave('Figures/Bar_Mean_Pepisode.pdf')

# save data
save(file = 'Rda/allAnalyzedData.Rda', list = c('meanPepisode', 'goodPepisode', 'navBlankShort', 'teleBlank', 'navBlankShortStats', 'teleBlankStats'))