library(dplyr)
library(tidyr)
library(reshape2)

load('Rda/allClassificationResults.Rda')

allClassificationResults <- separate(allClassificationResults, ElectrodeID, c('Subject', 'Session', 'Electrode'), sep = '_', remove = FALSE)
allClassificationResults <- allClassificationResults %>% 
  mutate(Subj_Elec = paste(Subject, Electrode, sep = '_'))
allClassificationResults$Subj_Elec <- as.factor(allClassificationResults$Subj_Elec)

patientBreakdown <- allClassificationResults %>%
  mutate(Signif = ifelse(CorrP < 0.05, 'TRUE', 'FALSE')) %>%
  group_by(Subject, Model, Signif) %>%
  summarise(Count = n())

# examine whether same electrode showed same result across sessions
sessionBreakdown <- allClassificationResults %>%
  mutate(Signif = ifelse(CorrP < 0.05, 'TRUE', 'FALSE')) %>%
  filter(Subject != 'UCDMC13') %>%
  arrange(Model, Subject, Electrode, Session) %>%
  select(Model, Subject, Session, Electrode, Signif) %>%
  dcast(Model + Subject + Electrode ~ Session, value.var = 'Signif') %>%
  mutate(MatchSignif = ifelse(TeleporterA == TeleporterB, 'TRUE', 'FALSE'))
sessionBreakdownSummary <- sessionBreakdown %>%
  group_by(Model, MatchSignif) %>%
  summarise(Count = n())