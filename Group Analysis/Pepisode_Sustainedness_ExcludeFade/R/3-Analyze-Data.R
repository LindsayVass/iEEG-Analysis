# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         23 November 2015
# Purpose:      This script will analyze the data from 2-Clean-Data.R. It will
#               perform the following analyses:
#               1. Compare pepisode between teleporter and navigation epochs. In
#                  this analysis, we will simply compare pepisode at each time
#                  bin (pre/tele/post) for the teleporter epoch and the preceding
#                  active navigation epoch. 

library(permute)
library(reshape2)
library(dplyr)
library(broom)
library(data.table)
library(coin)
load('Rda/allCleanData.Rda')

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


# Save data ---------------------------------------------------------------

save(file = 'Rda/allAnalyzedData_COIN.Rda', list = c('meanPepisode', 'pepisodePermData', 'pepisodePermNSigElectrodes', 'pepisodeTrueData', 'pepisodeTrueNSigElectrodes'))
