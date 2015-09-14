# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         23 July 2015
# Purpose:      This script will analyze the data from 2-Clean-Data.R. It will
#               perform the following analyses:
#               1. Compare pepisode between teleporter and navigation epochs. In
#                  this analysis, we will simply compare pepisode at each time
#                  bin (pre/tele/post) for the teleporter epoch and the preceding
#                  active navigation epoch. 
#               2. Compare the duration of oscillatory episodes for teleporter
#                  and navigation epochs. This analysis will identify episodes
#                  that cross a boundary (teleporter entry) and
#                  will test whether the mean duration of these episodes AFTER 
#                  the boundary differs between teleporter epochs and their 
#                  matched navigation epochs.
#                  NB There is no true "boundary" event for the navigation epochs;
#                  we are simply using the same temporal event structure as the
#                  teleporter epochs.



library(permute)
library(reshape2)
library(dplyr)
library(broom)
library(data.table)
library(coin)
load('Rda/allCleanData.Rda')

# Functions ---------------------------------------------------------------

meanPostEventOscDuration <- function(dataFrame, eventTime, trialTimeType) {
  output <- dataFrame %>%
    filter(Onset < eventTime & Offset > eventTime & TrialTimeType == trialTimeType) %>%
    mutate(PostEventDuration = Offset - eventTime) %>%
    group_by(ElectrodeID, FrequencyBand, RealTrialNumber) %>%
    summarise(MeanPostEventDuration = mean(PostEventDuration))
}

# Pepisode Analysis (1) ---------------------------------------------------
numPerm = 10000

meanPepisode <- allPepisode %>%
  summarise(MeanTelePepisode = mean(TelePepisode), MeanNavPepisode = mean(NavPepisode)) %>%
  ungroup() %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint)

# get permutation-corrected wilcoxon results for true data
pepisodeTrueData <- meanPepisode %>%
  do(NavGtTele = wilcoxsign_test(.$MeanNavPepisode ~ .$MeanTelePepisode, distribution = approximate(B = numPerm), alternative = "greater")) %>%
  mutate(Statistic = statistic(NavGtTele),
         PValue = pvalue(NavGtTele))
pepisodeTrueData$Statistic <- as.numeric(pepisodeTrueData$Statistic)
pepisodeTrueData$PValue    <- as.numeric(pepisodeTrueData$PValue)

pepisodeTrueNSigElectrodes <- pepisodeTrueData %>%
  group_by(FrequencyBand, TimePoint) %>%
  filter(PValue < 0.05) %>%
  summarise(Count = n())

# Determine number of electrodes significant by chance 
elecPerm = 1000

pepisodePermData <- vector(mode = "list", length = elecPerm)
pb <- txtProgressBar(min = 1, max = elecPerm, style = 3)
for (i in 1:elecPerm) {
  setTxtProgressBar(pb, i)
  pepisodePermData[[i]] <- pepisodeTrueData %>%
    select(-(Statistic:PValue)) %>%
    mutate(PValue = pperm(NavGtTele, rperm(NavGtTele, 1)),
           Iteration = i) %>%
    select(-NavGtTele)
}

pepisodePermData <- as.data.frame(rbindlist(pepisodePermData))
pepisodePermData$PValue <- as.numeric(pepisodePermData$PValue)

pepisodePermNSigElectrodes <- pepisodePermData %>%
  group_by(FrequencyBand, TimePoint, Iteration) %>%
  filter(PValue < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Iteration) %>%
  summarise(MaxCount = max(Count))

electrodeP <- data.frame(Iteration = 0,
                         MaxCount = 1:max(pepisodeTrueNSigElectrodes$Count))
electrodePCorr <- electrodeP %>%
  select(-Iteration) %>%
  rename(Count = MaxCount) %>%
  mutate(PValueElec = NA)

for (i in 1:nrow(electrodeP)) {
  thisData <- rbind(pepisodePermNSigElectrodes, electrodeP[i, ]) %>%
    arrange(desc(MaxCount))
  electrodePCorr$PValueElec[i] <- min(which(thisData$MaxCount == electrodeP$MaxCount[i])) / nrow(thisData)
}
remove(electrodeP)

pepisodeTrueNSigElectrodes <- inner_join(pepisodeTrueNSigElectrodes, electrodePCorr)

# Oscillatory episode post-event duration analysis (2) ------------------------

minTrials = 5 # must have this many trials of each condition for analysis

# get oscillatory event durations post-teleporter-entry
navNtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Navigation")
navFtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Navigation")

teleNtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Teleportation")
teleFtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Teleportation")

meanPostEntryDur <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType, Condition) %>%
  filter(n() > minTrials) %>%
  summarise(MeanDuration = mean(MeanPostEventDuration),
            SEM = sd(MeanPostEventDuration) / sqrt(n()))

# get permutation-corrected wilcoxon results
postEntryObs <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType, Condition) %>%
  select(-RealTrialNumber) %>%
  dcast(ElectrodeID + FrequencyBand + TimeType ~ Condition, fun.aggregate = length, value.var = "MeanPostEventDuration") %>%
  filter(Navigation >= minTrials & Teleportation >= minTrials) %>%
  select(-c(Navigation, Teleportation))
postEntryTrueData <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  select(-RealTrialNumber) %>%
  inner_join(postEntryObs) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType) %>%
  mutate_each(funs(as.factor), TimeType, Condition) %>%
  do(NavGtTele = wilcox_test(MeanPostEventDuration ~ Condition, data = ., distribution = approximate(B = numPerm), alternative = "greater")) %>%
  mutate(Statistic = statistic(NavGtTele),
         PValue = pvalue(NavGtTele)) %>%
  mutate_each(funs(as.numeric), Statistic, PValue)

postEntryTrueNSigElectrodes <- postEntryTrueData %>%
  group_by(FrequencyBand, TimeType) %>%
  filter(PValue < 0.05) %>%
  summarise(Count = n())


# Determine number of electrodes significant by chance 

postEntryPermData <- vector(mode = "list", length = elecPerm)
pb <- txtProgressBar(min = 1, max = elecPerm, style = 3)
for (i in 1:elecPerm) {
  setTxtProgressBar(pb, i)
  postEntryPermData[[i]] <- postEntryTrueData %>%
    select(-(Statistic:PValue)) %>%
    mutate(PValue = pperm(NavGtTele, rperm(NavGtTele, 1)),
           Iteration = i) %>%
    select(-NavGtTele)
}

postEntryPermData <- as.data.frame(rbindlist(postEntryPermData))
postEntryPermData$PValue <- as.numeric(postEntryPermData$PValue)

postEntryPermNSigElectrodes <- postEntryPermData %>%
  group_by(FrequencyBand, TimeType, Iteration) %>%
  filter(PValue < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Iteration) %>%
  summarise(MaxCount = max(Count))

electrodeP <- data.frame(Iteration = 0,
                         MaxCount = 1:max(postEntryTrueNSigElectrodes$Count))
electrodePCorr <- electrodeP %>%
  select(-Iteration) %>%
  rename(Count = MaxCount) %>%
  mutate(PValueElec = NA)
for (i in 1:nrow(electrodeP)) {
  thisData <- rbind(postEntryPermNSigElectrodes, electrodeP[i, ]) %>%
    arrange(desc(MaxCount))
  electrodePCorr$PValueElec[i] <- min(which(thisData$MaxCount == electrodeP$MaxCount[i])) / nrow(thisData)
}
remove(electrodeP)

postEntryTrueNSigElectrodes <- inner_join(postEntryTrueNSigElectrodes, electrodePCorr)

remove(electrodePCorr)

# Save data ---------------------------------------------------------------

save(file = 'Rda/allAnalyzedData_COIN.Rda', list = c('meanPepisode', 'pepisodePermData', 'pepisodePermNSigElectrodes', 'pepisodeTrueData', 'pepisodeTrueNSigElectrodes', 'postEntryPermData', 'postEntryPermNSigElectrodes','postEntryTrueData', 'postEntryTrueNSigElectrodes'))
