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
  group_by(ElectrodeID, Condition) %>%
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

# wilcoxon signed rank test for Nav vs Blank
navBlank <- goodPepisode %>%
  filter(Condition != "Teleportation_Tele")
navBlank$ElectrodeID <- factor(navBlank$ElectrodeID)
navBlank$Condition <- factor(navBlank$Condition)
navBlankStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = navBlank, distribution = approximate(B = 10000))

# wilcoxon signed rank test for Tele vs Blank
teleBlank <- goodPepisode %>%
  filter(Condition != "Navigation_Tele")
teleBlank$ElectrodeID <- factor(teleBlank$ElectrodeID)
teleBlank$Condition <- factor(teleBlank$Condition)
teleBlankStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = teleBlank, distribution = approximate(B = 10000))

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
save(file = 'Rda/allAnalyzedData.Rda', list = c('meanPepisode', 'goodPepisode', 'navBlank', 'teleBlank', 'navBlankStats', 'teleBlankStats'))