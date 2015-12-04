# Script Name:  3-Analyze-Plot-Data.R
# Author:       Lindsay Vass
# Date:         30 November 2015
# Purpose:      This script will test whether Pepisode differs between Blank &
#               Tele and between Blank & Nav using Wilcoxon signed rank tests.
#               It will then plot the mean pepisode across electrodes in a bar chart.

library(dplyr)
library(ggplot2)
library(coin)
library(ggthemes)

load('Rda/allCleanData_sepBlank.Rda')

# get mean value for each electrode/condition
meanPepisode <- allPepisode %>%
  filter(Condition == "BlankFree" | Condition == "BlankNav" | Condition == "Teleportation_Tele" | Condition == "Navigation_Tele") %>%
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
  mutate(Condition = "BlankFree") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

goodPepisode <- inner_join(meanPepisode, goodElec)
goodPepisode$ElectrodeID <- factor(goodPepisode$ElectrodeID)

##### wilcoxon signed rank test for Nav vs Blank (SHORT) (FREE EXPLORE)
navBlankShortFree <- meanPepisode %>%
  filter(Condition != "Teleportation_Tele",
         TrialTimeType == 'Short') 

# exclude electrodes that don't have data for "blank"
goodElec <- navBlankShortFree %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Navigation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankFree") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

navBlankShortFree <- inner_join(navBlankShortFree, goodElec) %>%
  filter(Condition != "BlankNav")

navBlankShortFree$ElectrodeID <- factor(navBlankShortFree$ElectrodeID)
navBlankShortFree$Condition <- factor(navBlankShortFree$Condition)
navBlankShortFreeStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = navBlankShortFree, distribution = approximate(B = 10000))

##### wilcoxon signed rank test for Nav vs Blank (SHORT) (NAVIGATION)
navBlankShortNav <- meanPepisode %>%
  filter(Condition != "Teleportation_Tele",
         TrialTimeType == 'Short') 

# exclude electrodes that don't have data for "blank"
goodElec <- navBlankShortNav %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Navigation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankNav") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

navBlankShortNav <- inner_join(navBlankShortNav, goodElec) %>%
  filter(Condition != "BlankFree")

navBlankShortNav$ElectrodeID <- factor(navBlankShortNav$ElectrodeID)
navBlankShortNav$Condition <- factor(navBlankShortNav$Condition)
navBlankShortNavStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = navBlankShortNav, distribution = approximate(B = 10000))


##### wilcoxon signed rank test for Nav vs Blank (LONG) (FREE EXPLORE)
navBlankLongFree <- meanPepisode %>%
  filter(Condition != "Teleportation_Tele",
         TrialTimeType == 'Long') 

# exclude electrodes that don't have data for "blank"
goodElec <- navBlankLongFree %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Navigation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankFree") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

navBlankLongFree <- inner_join(navBlankLongFree, goodElec) %>%
  filter(Condition != "BlankNav")

navBlankLongFree$ElectrodeID <- factor(navBlankLongFree$ElectrodeID)
navBlankLongFree$Condition <- factor(navBlankLongFree$Condition)
navBlankLongFreeStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = navBlankLongFree, distribution = approximate(B = 10000))

##### wilcoxon signed rank test for Nav vs Blank (LONG) (NAVIGATION)
navBlankLongNav <- meanPepisode %>%
  filter(Condition != "Teleportation_Tele",
         TrialTimeType == 'Long') 

# exclude electrodes that don't have data for "blank"
goodElec <- navBlankLongNav %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Navigation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankNav") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

navBlankLongNav <- inner_join(navBlankLongNav, goodElec) %>%
  filter(Condition != "BlankFree")

navBlankLongNav$ElectrodeID <- factor(navBlankLongNav$ElectrodeID)
navBlankLongNav$Condition <- factor(navBlankLongNav$Condition)
navBlankLongNavStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = navBlankLongNav, distribution = approximate(B = 10000))


##### wilcoxon signed rank test for Tele vs Blank (SHORT) (FREE EXPLORE)
teleBlankShortFree <- meanPepisode %>%
  filter(Condition != "Navigation_Tele",
         TrialTimeType == 'Short') 

# exclude electrodes that don't have data for "blank"
goodElec <- teleBlankShortFree %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Teleportation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankFree") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

teleBlankShortFree <- inner_join(teleBlankShortFree, goodElec) %>%
  filter(Condition != "BlankNav")

teleBlankShortFree$ElectrodeID <- factor(teleBlankShortFree$ElectrodeID)
teleBlankShortFree$Condition <- factor(teleBlankShortFree$Condition)
teleBlankShortFreeStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = teleBlankShortFree, distribution = approximate(B = 10000))

##### wilcoxon signed rank test for Tele vs Blank (SHORT) (NAVIGATION)
teleBlankShortNav <- meanPepisode %>%
  filter(Condition != "Navigation_Tele",
         TrialTimeType == 'Short') 

# exclude electrodes that don't have data for "blank"
goodElec <- teleBlankShortNav %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Teleportation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankNav") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

teleBlankShortNav <- inner_join(teleBlankShortNav, goodElec) %>%
  filter(Condition != "BlankFree")

teleBlankShortNav$ElectrodeID <- factor(teleBlankShortNav$ElectrodeID)
teleBlankShortNav$Condition <- factor(teleBlankShortNav$Condition)
teleBlankShortNavStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = teleBlankShortNav, distribution = approximate(B = 10000))


##### wilcoxon signed rank test for Tele vs Blank (LONG) (FREE EXPLORE)
teleBlankLongFree <- meanPepisode %>%
  filter(Condition != "Navigation_Tele",
         TrialTimeType == 'Long') 

# exclude electrodes that don't have data for "blank"
goodElec <- teleBlankLongFree %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Teleportation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankFree") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

teleBlankLongFree <- inner_join(teleBlankLongFree, goodElec) %>%
  filter(Condition != "BlankNav")

teleBlankLongFree$ElectrodeID <- factor(teleBlankLongFree$ElectrodeID)
teleBlankLongFree$Condition <- factor(teleBlankLongFree$Condition)
teleBlankLongFreeStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = teleBlankLongFree, distribution = approximate(B = 10000))

##### wilcoxon signed rank test for Tele vs Blank (LONG) (Nav)
teleBlankLongNav <- meanPepisode %>%
  filter(Condition != "Navigation_Tele",
         TrialTimeType == 'Long') 

# exclude electrodes that don't have data for "blank"
goodElec <- teleBlankLongNav %>%
  ungroup() %>%
  select(ElectrodeID, Condition)  
badElec <- goodElec %>%
  filter(Condition == "Teleportation_Tele") %>%
  select(ElectrodeID) %>%
  unique() %>%
  mutate(Condition = "BlankNav") %>%
  anti_join(goodElec) %>%
  select(-Condition)
goodElec <- anti_join(goodElec, badElec)

teleBlankLongNav <- inner_join(teleBlankLongNav, goodElec) %>%
  filter(Condition != "BlankFree")

teleBlankLongNav$ElectrodeID <- factor(teleBlankLongNav$ElectrodeID)
teleBlankLongNav$Condition <- factor(teleBlankLongNav$Condition)
teleBlankLongNavStats <-  wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = teleBlankLongNav, distribution = approximate(B = 10000))



# plot the results
goodPepisode$Condition <- plyr::mapvalues(goodPepisode$Condition, from = c("BlankFree", "BlankNav", "Teleportation_Tele", "Navigation_Tele"), to = c("Blank-AfterNav", "Blank-AfterTele", "Teleportation", "Navigation"))
p <- goodPepisode %>%
  group_by(Condition, TrialTimeType) %>%
  summarise(GroupMeanPepisode = mean(MeanPepisode),
            GroupSEMPepisode = sd(MeanPepisode) / n()) %>%
  ggplot(aes(x = Condition, y = GroupMeanPepisode, ymin = GroupMeanPepisode - GroupSEMPepisode, ymax = GroupMeanPepisode + GroupSEMPepisode)) +
  facet_wrap(~TrialTimeType) +
  geom_bar(stat = "identity") +
  geom_errorbar() +
  theme_few() +
  ylab(expression(paste(Mean~P[Episode]))) +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank())
ggsave('Figures/Bar_Mean_Pepisode_BlankType_TimeType.pdf')

pElec <- goodPepisode %>%
  group_by(Condition, TrialTimeType, ElectrodeID) %>%
  ggplot(aes(x = Condition, y = MeanPepisode, ymin = MeanPepisode - SEMPepisode, ymax = MeanPepisode + SEMPepisode, group = ElectrodeID, color = ElectrodeID)) +
  facet_wrap(~TrialTimeType) +
  geom_line() +
  theme_few() +
  ylab(expression(paste(Mean~P[Episode]))) +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank()) +
  guides(color = FALSE)
ggsave('Figures/Line_Each_Elec_Mean_Pepisode_BlankType_TimeType.pdf')


# get mean value for each electrode/condition, including Post (still)
pepisode <- allPepisode %>%
  filter(Condition == "BlankFree" | Condition == "BlankNav" | Condition == "Teleportation_Tele" | Condition == "Navigation_Tele" | Condition == "Teleportation_Post1") %>%
  group_by(ElectrodeID, Condition, TrialTimeType) %>%
  summarise(MeanPepisode = mean(Pepisode),
            SEMPepisode = sd(Pepisode) / n())
pepisode$Condition <- factor(pepisode$Condition)
pepisode$Condition <- plyr::mapvalues(pepisode$Condition, 
                                      from = c("BlankFree", "BlankNav", "Teleportation_Tele", "Navigation_Tele", "Teleportation_Post1"),
                                      to = c("Blank-AfterNav", "Blank-AfterTele", "Teleportation", "Navigation", "Still/Slow"))
p2 <- pepisode %>%
  group_by(Condition, TrialTimeType) %>%
  summarise(GroupMeanPepisode = mean(MeanPepisode),
            GroupSEMPepisode = sd(MeanPepisode) / n()) %>%
  ggplot(aes(x = Condition, y = GroupMeanPepisode, ymin = GroupMeanPepisode - GroupSEMPepisode, ymax = GroupMeanPepisode + GroupSEMPepisode)) +
  facet_wrap(~TrialTimeType) +
  geom_bar(stat = "identity") +
  geom_errorbar() +
  theme_few() +
  ylab(expression(paste(Mean~P[Episode]))) +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank())
ggsave('Figures/Bar_Mean_Pepisode_BlankType_TimeType_Still.pdf')

# save data
save(file = 'Rda/allAnalyzedData_BlankType.Rda', list = c('meanPepisode', 'goodPepisode', 'navBlankShortFreeStats', 'navBlankShortNavStats', 'navBlankLongNavStats', 'navBlankLongFreeStats', 'teleBlankShortFreeStats', 'teleBlankShortNavStats', 'teleBlankLongFreeStats', 'teleBlankLongNavStats'))