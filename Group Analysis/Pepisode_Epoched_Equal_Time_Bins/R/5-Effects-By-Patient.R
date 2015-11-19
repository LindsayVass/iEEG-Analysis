library(tidyr)
library(dplyr)
library(reshape2)

load('Rda/allAnalyzedData.Rda')

moltenTrueData <- separate(moltenTrueData, ElectrodeID, c('Subject', 'Session', 'Electrode'), sep = '_', remove = FALSE)

subjectEffectBreakdown <- moltenTrueData %>%
  mutate(Signif = ifelse(PValue < 0.05, 'TRUE', 'FALSE')) %>%
  group_by(Subject, FrequencyBand, Contrast, Signif) %>%
  summarise(Count = n())

sessionBreakdown <- moltenTrueData %>%
  mutate(Signif = ifelse(PValue < 0.05, 'TRUE', 'FALSE'),
         Subj_Elec = paste(Subject, Electrode, sep = '_')) %>%
  filter(FrequencyBand == 'Delta-Theta') %>%
  dcast(Subj_Elec + Contrast ~ Session, value.var = 'Signif') %>%
  mutate(SignifMatch = ifelse(TeleporterA == TeleporterB, 'TRUE', 'FALSE'))
sessionBreakdown$Subj_Elec <- as.factor(sessionBreakdown$Subj_Elec)

sessionBreakdownSummary <- sessionBreakdown %>%
  group_by(Contrast, SignifMatch) %>%
  summarise(Count = n())