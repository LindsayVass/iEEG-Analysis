# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         23 November 2015
# Purpose:      This script will analyze the data frame produced by 
#               "2-Clean-Data.R". For each electrode and frequency band, it will 
#               use a Wilcoxon signed-rank test to determine whether pepisode 
#               significantly differs between timepoints. 

library(dplyr)
library(reshape2)
library(data.table)
library(coin)

load('Rda/allCleanData.Rda')


# Functions ---------------------------------------------------------------
tidyWilcoxon <- function(inputData) {
  preGtTele <- tidy(inputData, PreGtTele) %>%
    mutate(Contrast = "Pre > Tele")
  preLtTele <- tidy(inputData, PreLtTele) %>%
    mutate(Contrast = "Pre < Tele")
  teleGtPost <- tidy(inputData, TeleGtPost) %>%
    mutate(Contrast = "Tele > Post")
  teleLtPost <- tidy(inputData, TeleLtPost) %>%
    mutate(Contrast = "Tele < Post")
  output <- rbind(preGtTele, preLtTele, teleGtPost, teleLtPost)
}

replicateData <- function(inputData, nrep) {
  repData <- inputData[rep(1:nrow(inputData), times = nrep), ] %>%
    mutate(Iteration = gl(nrep, nrow(inputData)))
}

shuffleData <- function(thisData, control) {
  shuffleOrder <- shuffle(nrow(thisData), control = control)
  thisData$MeanPepisode <- thisData$MeanPepisode[shuffleOrder]
  return(thisData)
}

getCorrectedP <- function(permData, trueData, trueOrPerm = "true") {
  
  if (tolower(trueOrPerm) == "perm") {
    origData <- trueData
    trueData <- trueData %>%
      ungroup() %>%
      select(-c(Iteration, ElectrodeIteration))
  }
  
  trueData <- trueData %>%
    mutate(Iteration = 0)
  trueData$Iteration <- factor(trueData$Iteration)
  
  trueDataCorrected <- rbind(permData, trueData)
  trueDataCorrected$Contrast <- factor(trueDataCorrected$Contrast)
  trueDataCorrected <- trueDataCorrected %>%
    ungroup() %>%
    group_by(ElectrodeID, FrequencyBand, Contrast) 
  
  trueDataCorrectedGt <- trueDataCorrected %>%
    arrange(desc(statistic)) %>%
    mutate(CorrP = row_number(desc(statistic)) / n()) %>%
    filter(Iteration == 0)
  trueDataCorrectedLt <- trueDataCorrected %>%
    arrange(statistic) %>%
    mutate(CorrP = row_number(statistic) / n()) %>%
    filter(Iteration == 0)
  
  if (tolower(trueOrPerm) == "perm") {
    trueDataGt <- trueData %>%
      filter(Contrast == "Pre > Tele" | Contrast == "Tele > Post") %>%
      inner_join(trueDataCorrectedGt, by = c("ElectrodeID", "FrequencyBand", "statistic", "p.value", "Contrast", "Iteration")) %>%
      select(-Iteration) %>%
      inner_join(origData, by = c("ElectrodeID", "FrequencyBand", "statistic", "p.value", "Contrast"))
    trueDataLt <- trueData %>%
      filter(Contrast == "Pre < Tele" | Contrast == "Tele < Post") %>%
      inner_join(trueDataCorrectedLt, by = c("ElectrodeID", "FrequencyBand", "statistic", "p.value", "Contrast", "Iteration")) %>%
      select(-Iteration) %>%
      inner_join(origData, by = c("ElectrodeID", "FrequencyBand", "statistic", "p.value", "Contrast"))
  } else {
    trueDataGt <- trueData %>%
      filter(Contrast == "Pre > Tele" | Contrast == "Tele > Post") %>%
      inner_join(trueDataCorrectedGt, by = c("ElectrodeID", "FrequencyBand", "statistic", "p.value", "Contrast", "Iteration"))
    trueDataLt <- trueData %>%
      filter(Contrast == "Pre < Tele" | Contrast == "Tele < Post") %>%
      inner_join(trueDataCorrectedLt, by = c("ElectrodeID", "FrequencyBand", "statistic", "p.value", "Contrast", "Iteration"))
  }
  
  trueDataCorrected <- rbind(trueDataGt, trueDataLt)
  return(trueDataCorrected)
}

getElectrodeCorrectedP <- function(trueNSigElectrodes, permNSigElectrodes) {
  permNSigElectrodes <- permNSigElectrodes %>%
    c(trueNSigElectrodes) %>%
    sort(decreasing = TRUE)
  CorrP <- min(which(permNSigElectrodes == trueNSigElectrodes)) / length(permNSigElectrodes)
}

# Wilcoxon signed-rank at each electrode and frequency band ---------------

# permutation testing will be automatically calculated using coin's version of the test
numPerm=10000

cleanData <- cleanData %>%
  group_by(ElectrodeID, FrequencyBand, TrialNumber, TimePoint) %>%
  summarise(MeanPepisode = mean(Pepisode)) 

trueData <- cleanData %>%
  group_by(ElectrodeID, FrequencyBand, TrialNumber, TimePoint) %>%
  dcast(TrialNumber + ElectrodeID + FrequencyBand ~ TimePoint, value.var = "MeanPepisode") %>%
  group_by(ElectrodeID, FrequencyBand) %>%
  do(TeleLtPre  = wilcoxsign_test(.$Tele ~ .$Pre1,  distribution = approximate(B=numPerm), alternative = "less"),
     TeleGtPre  = wilcoxsign_test(.$Tele ~ .$Pre1,  distribution = approximate(B=numPerm), alternative = "greater"),
     TeleGtPost = wilcoxsign_test(.$Tele ~ .$Post1, distribution = approximate(B=numPerm), alternative = "greater"),
     TeleLtPost = wilcoxsign_test(.$Tele ~ .$Post1, distribution = approximate(B=numPerm), alternative = "less")) %>%
  mutate_each(funs(statistic, pvalue), TeleLtPre, TeleGtPre, TeleLtPost, TeleGtPost)

# put the data into molten format
moltenP <- trueData %>%
  select(-c(TeleLtPre, TeleGtPre, TeleGtPost, TeleLtPost)) %>%
  select(-(ends_with("statistic"))) %>%
  melt(id.vars = c("ElectrodeID", "FrequencyBand"), variable.name = "Contrast", value.name = "PValue")
moltenP$Contrast <- sub("_pvalue", "", moltenP$Contrast)
moltenZ <- trueData %>%
  select(-c(TeleLtPre, TeleGtPre, TeleGtPost, TeleLtPost)) %>%
  select(-(ends_with("pvalue"))) %>%
  melt(id.vars = c("ElectrodeID", "FrequencyBand"), variable.name = "Contrast", value.name = "Statistic")
moltenZ$Contrast <- sub("_statistic", "", moltenZ$Contrast)
moltenTrueData <- inner_join(moltenZ, moltenP)
remove(moltenP, moltenZ)

# how many electrodes were significant for each contrast
trueNSigElectrodes <- moltenTrueData %>%
  group_by(FrequencyBand, Contrast) %>%
  filter(PValue < 0.05) %>%
  summarise(Count = n())


# Determine number of electrodes significant by chance --------------------
elecPerm = 1000

permData <- vector(mode = "list", length = elecPerm)
pb <- txtProgressBar(min = 1, max = elecPerm, style = 3)
for (i in 1:elecPerm) {
  setTxtProgressBar(pb, i)
  permData[[i]] <- trueData %>%
    select(-(7:14)) %>%
    mutate_each(funs(list(pperm(., rperm(., 1)))), TeleLtPre, TeleGtPre, TeleGtPost, TeleLtPost) %>%
    mutate(Iteration = i)
}

permData <- as.data.frame(rbindlist(permData)) %>%
  mutate_each(funs(as.numeric), TeleLtPre, TeleGtPre, TeleGtPost, TeleLtPost)

permNSigElectrodes <- permData %>%
  melt(id.vars = c("ElectrodeID", "FrequencyBand", "Iteration"), variable.name = "Contrast", value.name = "PValue") %>%
  group_by(FrequencyBand, Iteration, Contrast) %>%
  filter(PValue < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Iteration) %>%
  summarise(MaxCount = max(Count))

electrodeP <- data.frame(Iteration = 0,
                         MaxCount = 1:max(trueNSigElectrodes$Count))
electrodePCorr <- electrodeP %>%
  select(-Iteration) %>%
  rename(Count = MaxCount) %>%
  mutate(PValueElec = NA)
for (i in 1:nrow(electrodeP)) {
  thisData <- rbind(permNSigElectrodes, electrodeP[i, ]) %>%
    arrange(desc(MaxCount)) 
  electrodePCorr$PValueElec[i] <- min(which(thisData$MaxCount == electrodeP$MaxCount[i])) / nrow(thisData)
}
remove(electrodeP)

trueNSigElectrodes <- inner_join(trueNSigElectrodes, electrodePCorr)

save(file = 'Rda/allAnalyzedData.Rda', list = c('cleanData', 'electrodePCorr', 'moltenTrueData', 'permData', 'permNSigElectrodes', 'trueData', 'trueNSigElectrodes'))