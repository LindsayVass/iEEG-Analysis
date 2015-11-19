library(tidyr)
library(dplyr)
library(reshape2)

load('Rda/allAnalyzedData_COIN.Rda')

pepisodeTrueData <- separate(pepisodeTrueData, ElectrodeID, c('Subject', 'Session', 'Electrode'), sep = '_', remove = FALSE)

pepisodeBySubject <- pepisodeTrueData %>%
  mutate(Signif = ifelse(PValue < 0.05, 'TRUE', 'FALSE')) %>%
  filter(FrequencyBand == 'Delta-Theta',
         TimePoint == 'Post1') %>%
  group_by(Subject, Signif) %>%
  summarise(Count = n())

sessionBreakdown <- pepisodeTrueData %>%
  mutate(Signif = ifelse(PValue < 0.05, 'TRUE', 'FALSE'),
         Subj_Elec = paste(Subject, Electrode, sep = '_')) %>%
  filter(FrequencyBand == 'Delta-Theta') %>%
  dcast(Subj_Elec + TimePoint ~ Session, value.var = 'Signif') %>%
  mutate(SignifMatch = ifelse(TeleporterA == TeleporterB, 'TRUE', 'FALSE'))
sessionBreakdown$Subj_Elec <- as.factor(sessionBreakdown$Subj_Elec)

sessionBreakdownSummary <- sessionBreakdown %>%
  group_by(TimePoint, SignifMatch) %>%
  summarise(Count = n())