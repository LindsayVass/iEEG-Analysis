library(dplyr)
library(tidyr)
library(reshape2)

load('Rda/allClassificationResults.Rda')

allMeanClassificationResults <- separate(allMeanClassificationResults, ElectrodeID, c('Subject', 'Session', 'Electrode'), sep = '_', remove = FALSE)

patientBreakdown <- allMeanClassificationResults %>%
  mutate(Signif = ifelse(CorrP < 0.05, 'TRUE', 'FALSE')) %>%
  group_by(Subject, Model, Signif) %>%
  summarise(Count = n())


sessionBreakdown <- allMeanClassificationResults %>%
  mutate(Signif = ifelse(CorrP < 0.05, 'TRUE', 'FALSE'),
         Subj_Elec = paste(Subject, Electrode, sep = '_')) %>%
  dcast(Subj_Elec + Subject + Model ~ Session, value.var = 'Signif') %>%
  mutate(SignifMatch = ifelse(TeleporterA == TeleporterB, 'TRUE', 'FALSE'))
sessionBreakdown$Subj_Elec <- as.factor(sessionBreakdown$Subj_Elec)

sessionBreakdownSummary <- sessionBreakdown %>%
  filter(Subject != 'UCDMC13') %>%
  group_by(Model, SignifMatch) %>%
  summarise(Count = n())
