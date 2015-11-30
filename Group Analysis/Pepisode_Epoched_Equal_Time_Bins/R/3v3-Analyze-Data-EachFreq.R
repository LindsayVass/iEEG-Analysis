# Script Name:  3v3-Analyze-Data-EachFreq.R
# Author:       Lindsay Vass
# Date:         19 November 2015
# Purpose:      For each electrode and frequency, it will 
#               use a Wilcoxon signed-rank test to determine whether pepisode 
#               significantly differs between timepoints. Unlike other analyses,
#               this one does not group by frequency band.

library(dplyr)
library(reshape2)
library(data.table)
library(coin)
library(ggplot2)
library(ggthemes)

load('Rda/allCleanData_EachFreq.Rda')


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
    group_by(ElectrodeID, Frequency, Contrast) 
  
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
      inner_join(trueDataCorrectedGt, by = c("ElectrodeID", "Frequency", "statistic", "p.value", "Contrast", "Iteration")) %>%
      select(-Iteration) %>%
      inner_join(origData, by = c("ElectrodeID", "Frequency", "statistic", "p.value", "Contrast"))
    trueDataLt <- trueData %>%
      filter(Contrast == "Pre < Tele" | Contrast == "Tele < Post") %>%
      inner_join(trueDataCorrectedLt, by = c("ElectrodeID", "Frequency", "statistic", "p.value", "Contrast", "Iteration")) %>%
      select(-Iteration) %>%
      inner_join(origData, by = c("ElectrodeID", "Frequency", "statistic", "p.value", "Contrast"))
  } else {
    trueDataGt <- trueData %>%
      filter(Contrast == "Pre > Tele" | Contrast == "Tele > Post") %>%
      inner_join(trueDataCorrectedGt, by = c("ElectrodeID", "Frequency", "statistic", "p.value", "Contrast", "Iteration"))
    trueDataLt <- trueData %>%
      filter(Contrast == "Pre < Tele" | Contrast == "Tele < Post") %>%
      inner_join(trueDataCorrectedLt, by = c("ElectrodeID", "Frequency", "statistic", "p.value", "Contrast", "Iteration"))
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

# Wilcoxon signed-rank at each electrode and frequency -------------------------

# permutation testing will be automatically calculated using coin's version of the test
numPerm=10000

cleanData <- cleanData %>%
  group_by(ElectrodeID, Frequency, TrialNumber, TimePoint) %>%
  summarise(MeanPepisode = mean(Pepisode)) 

trueData <- cleanData %>%
  group_by(ElectrodeID, Frequency, TrialNumber, TimePoint) %>%
  dcast(TrialNumber + ElectrodeID + Frequency ~ TimePoint, value.var = "MeanPepisode") %>%
  group_by(ElectrodeID, Frequency) %>%
  do(TeleLtPre  = tryCatch(wilcoxsign_test(.$Tele ~ .$Pre1,  distribution = approximate(B=numPerm), alternative = "less"), error = function (e) NA),
     TeleGtPre  = tryCatch(wilcoxsign_test(.$Tele ~ .$Pre1,  distribution = approximate(B=numPerm), alternative = "greater"), error = function (e) NA),
     TeleGtPost = tryCatch(wilcoxsign_test(.$Tele ~ .$Post1, distribution = approximate(B=numPerm), alternative = "greater"), error = function (e) NA),
     TeleLtPost = tryCatch(wilcoxsign_test(.$Tele ~ .$Post1, distribution = approximate(B=numPerm), alternative = "less"), error = function (e) NA))

# can't melt stats, so separate the data frames manually
teleLtPre <- trueData %>%
  select(ElectrodeID, Frequency, TeleLtPre) 
goodInds <- which(is.na(teleLtPre$TeleLtPre) == FALSE)
badInds<- which(is.na(teleLtPre$TeleLtPre) == TRUE)
teleLtPreGood <- teleLtPre[goodInds, ] %>%
  rowwise() %>%
  do(statistic = statistic(.$TeleLtPre),
     pvalue = pvalue(.$TeleLtPre)) %>%
  cbind(teleLtPre[goodInds, 1:2])
teleLtPreBad <- teleLtPre[badInds, ] %>%
  mutate(statistic = 0,
         pvalue = 1) %>%
  select(-TeleLtPre)
teleLtPre <- rbind(teleLtPreGood, teleLtPreBad) %>%
  mutate(Contrast = 'TeleLtPre')

teleLtPost <- trueData %>%
  select(ElectrodeID, Frequency, TeleLtPost) 
goodInds <- which(is.na(teleLtPost$TeleLtPost) == FALSE)
badInds<- which(is.na(teleLtPost$TeleLtPost) == TRUE)
teleLtPostGood <- teleLtPost[goodInds, ] %>%
  rowwise() %>%
  do(statistic = statistic(.$TeleLtPost),
     pvalue = pvalue(.$TeleLtPost)) %>%
  cbind(teleLtPost[goodInds, 1:2])
teleLtPostBad <- teleLtPost[badInds, ] %>%
  mutate(statistic = 0,
         pvalue = 1) %>%
  select(-TeleLtPost)
teleLtPost <- rbind(teleLtPostGood, teleLtPostBad) %>%
  mutate(Contrast = 'TeleLtPost')

teleGtPre <- trueData %>%
  select(ElectrodeID, Frequency, TeleGtPre) 
goodInds <- which(is.na(teleGtPre$TeleGtPre) == FALSE)
badInds<- which(is.na(teleGtPre$TeleGtPre) == TRUE)
teleGtPreGood <- teleGtPre[goodInds, ] %>%
  rowwise() %>%
  do(statistic = statistic(.$TeleGtPre),
     pvalue = pvalue(.$TeleGtPre)) %>%
  cbind(teleGtPre[goodInds, 1:2])
teleGtPreBad <- teleGtPre[badInds, ] %>%
  mutate(statistic = 0,
         pvalue = 1) %>%
  select(-TeleGtPre)
teleGtPre <- rbind(teleGtPreGood, teleGtPreBad) %>%
  mutate(Contrast = 'TeleGtPre')

teleGtPost <- trueData %>%
  select(ElectrodeID, Frequency, TeleGtPost) 
goodInds <- which(is.na(teleGtPost$TeleGtPost) == FALSE)
badInds<- which(is.na(teleGtPost$TeleGtPost) == TRUE)
teleGtPostGood <- teleGtPost[goodInds, ] %>%
  rowwise() %>%
  do(statistic = statistic(.$TeleGtPost),
     pvalue = pvalue(.$TeleGtPost)) %>%
  cbind(teleGtPost[goodInds, 1:2])
teleGtPostBad <- teleGtPost[badInds, ] %>%
  mutate(statistic = 0,
         pvalue = 1) %>%
  select(-TeleGtPost)
teleGtPost <- rbind(teleGtPostGood, teleGtPostBad) %>%
  mutate(Contrast = 'TeleGtPost')

moltenTrueData <- rbind(teleGtPost, teleLtPost, teleGtPre, teleLtPre)
moltenTrueData$Contrast <- as.factor(moltenTrueData$Contrast)

# how many electrodes were significant for each contrast
trueNSigElectrodes <- moltenTrueData %>%
  group_by(Frequency, Contrast) %>%
  filter(pvalue < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  arrange(Contrast, Frequency)
trueZero <- moltenTrueData %>%
  select(Contrast, Frequency) %>%
  unique() %>%
  anti_join(trueNSigElectrodes) %>%
  mutate(Count = 0)
trueNSigElectrodes <- rbind(trueNSigElectrodes, trueZero)

trueNSigElectrodes$Contrast <- plyr::mapvalues(trueNSigElectrodes$Contrast, 
                                               from = c('TeleGtPost', 'TeleGtPre', 'TeleLtPost', 'TeleLtPre'), 
                                               to = c('Tele > Post', 'Tele > Pre', 'Tele < Post', 'Tele < Pre'))

# Determine number of electrodes significant by chance --------------------
elecPerm = 1000

# first get rid of rows containing NA
validData <- trueData
goodTeleLtPre <- which(is.na(validData$TeleLtPre) == FALSE)
validData <- validData[goodTeleLtPre, ]
goodTeleLtPost <- which(is.na(validData$TeleLtPost) == FALSE)
validData <- validData[goodTeleLtPost, ]
goodTeleGtPre <- which(is.na(validData$TeleGtPre) == FALSE)
validData <- validData[goodTeleGtPre, ]
goodTeleGtPost <- which(is.na(validData$TeleGtPost) == FALSE)
validData <- validData[goodTeleGtPost, ]

permData <- vector(mode = "list", length = elecPerm)
pb <- txtProgressBar(min = 1, max = elecPerm, style = 3)
for (i in 1:elecPerm) {
  setTxtProgressBar(pb, i)
  tmp <- validData %>%
    rowwise() %>%
    do(TeleLtPre_pvalue = pperm(.$TeleLtPre, rperm(.$TeleLtPre, 1)),
       TeleGtPre_pvalue = pperm(.$TeleGtPre, rperm(.$TeleGtPre, 1)),
       TeleLtPost_pvalue = pperm(.$TeleLtPost, rperm(.$TeleLtPost, 1)),
       TeleGtPost_pvalue = pperm(.$TeleGtPost, rperm(.$TeleGtPost, 1))) %>%
    mutate(Iteration = i) %>%
    cbind(validData) %>%
    select(-c(TeleLtPre, TeleLtPost, TeleGtPre, TeleGtPost))
  permData[[i]] <- tmp
}

permDataDf <- as.data.frame(rbindlist(permData)) %>%
  mutate_each(funs(as.numeric), TeleLtPre_pvalue, TeleGtPre_pvalue, TeleGtPost_pvalue, TeleLtPost_pvalue)

permNSigElectrodes <- permDataDf %>%
  melt(id.vars = c("ElectrodeID", "Frequency", "Iteration"), variable.name = "Contrast", value.name = "PValue") %>%
  group_by(Frequency, Iteration, Contrast) %>%
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

save(file = 'Rda/allAnalyzedData_EachFreq.Rda', list = c('cleanData', 'electrodePCorr', 'moltenTrueData', 'permData', 'permNSigElectrodes', 'trueData', 'validData', 'trueNSigElectrodes'))

## Plot Results
threshP <- electrodePCorr$Count[max(which(electrodePCorr$PValueElec > 0.05))] + 0.5

trueZero <- moltenTrueData %>%
  select(Contrast, Frequency) %>%
  unique()
trueZero$Contrast <- plyr::mapvalues(trueZero$Contrast, 
                                               from = c('TeleGtPost', 'TeleGtPre', 'TeleLtPost', 'TeleLtPre'), 
                                               to = c('Tele > Post', 'Tele > Pre', 'Tele < Post', 'Tele < Pre'))
trueZero <- anti_join(trueZero, trueNSigElectrodes) %>%
  mutate(Count = 0,
         PValueElec = 1)
trueNSigElectrodes <- rbind(trueNSigElectrodes, trueZero)

ggplot(trueNSigElectrodes, aes(x = Frequency, y = Count)) +
  geom_point(size = 4) +
  geom_hline(aes(yintercept = threshP), color = 'red') +
  scale_x_log10(breaks = c(2,4,8,16,32,64,128)) +
  scale_y_continuous(breaks = c(2,4,6,8,10,12,14)) +
  facet_wrap(~Contrast) +
  ylab('# Significant Electrodes') +
  xlab("Frequency (Hz)") +
  theme(text = element_text(size = 24))
ggsave('Figures/TimePointContrast_EachFrequency.pdf')