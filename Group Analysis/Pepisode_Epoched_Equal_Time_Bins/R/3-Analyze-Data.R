# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         3 June 2015
# Purpose:      This script will analyze the data frame produced by 
#               "2-Clean-Data.R". For each electrode and frequency band, it will 
#               use a Wilcoxon signed-rank test to determine whether pepisode 
#               significantly differs between timepoints. 
#
# Updated 20 August 2015: Refactored and added permutation testing

library(dplyr)
library(reshape2)
library(broom)
library(permute)
library(data.table)

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

cleanData <- cleanData %>%
  group_by(ElectrodeID, FrequencyBand, TrialNumber, TimePoint) %>%
  summarise(MeanPepisode = mean(Pepisode)) 
trueData <- cleanData %>%
  group_by(ElectrodeID, FrequencyBand, TrialNumber, TimePoint) %>%
  dcast(TrialNumber + ElectrodeID + FrequencyBand ~ TimePoint, value.var = "MeanPepisode") %>%
  group_by(ElectrodeID, FrequencyBand) %>%
  do(PreGtTele  = wilcox.test(.$Pre1, .$Tele, alternative = "greater", paired = TRUE),
     PreLtTele  = wilcox.test(.$Pre1, .$Tele, alternative = "less", paired = TRUE),
     TeleGtPost = wilcox.test(.$Tele, .$Post1, alternative = "greater", paired = TRUE),
     TeleLtPost = wilcox.test(.$Tele, .$Post1, alternative = "less", paired = TRUE)) %>%
  tidyWilcoxon()  

# Run permuted wilcoxon signed-rank tests ---------------------------------

cleanData <- cleanData %>%
  mutate(Observation = paste(ElectrodeID, FrequencyBand, sep = "_"))
cleanData$Observation <- factor(cleanData$Observation)

nperm <- 1000
control <- how(plots = Plots(cleanData$TrialNumber), within = Within(type = "free"), blocks = cleanData$Observation)

permData <- cleanData %>%
  ungroup() %>%
  replicateData(nperm)

permResults <- vector(mode = "list", length = nperm)
pb <- txtProgressBar(min = 1, max = nperm, style = 3)
for (thisPerm in 1:nperm) {
  setTxtProgressBar(pb, thisPerm)
  permResults[[thisPerm]] <- permData %>%
    filter(Iteration == thisPerm) %>%
    shuffleData(control) %>%
    group_by(ElectrodeID, FrequencyBand, TrialNumber, TimePoint) %>%
    dcast(TrialNumber + ElectrodeID + FrequencyBand ~ TimePoint, value.var = "MeanPepisode") %>%
    group_by(ElectrodeID, FrequencyBand) %>%
    do(PreGtTele  = wilcox.test(.$Pre1, .$Tele, alternative = "greater", paired = TRUE),
       PreLtTele  = wilcox.test(.$Pre1, .$Tele, alternative = "less", paired = TRUE),
       TeleGtPost = wilcox.test(.$Tele, .$Post1, alternative = "greater", paired = TRUE),
       TeleLtPost = wilcox.test(.$Tele, .$Post1, alternative = "less", paired = TRUE)) %>%
    tidyWilcoxon() %>%
    mutate(Iteration = thisPerm)
}

permResults <- rbindlist(permResults)


# Add corrected P to wilcoxon results -------------------------------------

trueDataCorrected <- trueData %>%
  group_by(ElectrodeID, FrequencyBand, Contrast) %>%
  getCorrectedP(permResults, ., "true")
trueNSigElectrodes <- trueDataCorrected %>%
  ungroup() %>%
  group_by(FrequencyBand, Contrast) %>%
  filter(p.value < 0.05 & CorrP < 0.05) %>%
  summarise(Count = n())

# Determine # of significant electrodes by chance -------------------------

permAnalysisList <- permResults %>%
  select(ElectrodeID, FrequencyBand, Iteration) %>%
  group_by(ElectrodeID, FrequencyBand) %>%
  unique()

permMaxElectrodes <- vector(mode = "list", length = nperm)
for (thisPerm in 1:nperm) {
  permMaxElectrodes[[thisPerm]] <- permAnalysisList %>%
    sample_n(1) %>%
    mutate(ElectrodeIteration = thisPerm)
}

permMaxElectrodes <- rbindlist(permMaxElectrodes) %>%
  inner_join(permResults) %>%
  ungroup() %>%
  group_by(ElectrodeIteration) %>%
  do(getCorrectedP(permResults, ., "perm"))

permNSigElectrodes <- permMaxElectrodes %>%
  group_by(FrequencyBand, Contrast, ElectrodeIteration) %>%
  filter(CorrP < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(ElectrodeIteration) %>%
  summarise(MaxCount = max(Count))

# Get corrected P values for number of significant electrodes -------------
trueNSigElectrodesCorrected <- trueNSigElectrodes %>%
  group_by(FrequencyBand, Contrast, Count) %>%
  do(CorrP = getElectrodeCorrectedP(.$Count, permNSigElectrodes$MaxCount))



# Save data ---------------------------------------------------------------

save(file = 'Rda/allAnalyzedData.Rda', list = c('permResults', 'permMaxElectrodes', 'permNSigElectrodes', 'trueDataCorrected', 'trueNSigElectrodesCorrected'))