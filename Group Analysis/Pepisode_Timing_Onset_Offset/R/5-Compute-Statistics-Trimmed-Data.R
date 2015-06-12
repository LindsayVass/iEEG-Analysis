# Script Name:  5-Compute-Statistics-Trimmed-Data.R
# Author:       Lindsay Vass
# Date:         9 June 2015
# Purpose:      This script will calculate statistics for the group data plotted
#               in 4b-Plot-Trimmed-Data.R. For the episode timing, it will break
#               the time series into chunks and test whether pepisode varies
#               between the chunks using an ANOVA. Because we aim to prove the 
#               null hypothesis (no difference in pepisode between time points),
#               we will perform follow-up Bayesian statistics as described in 
#               Masson (2011) Behav Res "A tutorial on a practical Bayesian
#               alternaitve to null-hypothesis significance testing." For the
#               onsets and offsets, we will use the chi-square goodness of fit
#               test to determine whether the observed distribution varies from
#               a uniform distribution.

library(dplyr)
library(permute)
load('Rda/allTrimmedData.Rda')

# Functions ---------------------------------------------------------------

# return data from a specific frequency band and trial time type
filterRMData <- function(episodeData, timeType, frequencyBand) {
  filteredData <- episodeData %>%
    filter(TrialTimeType == timeType &
           FrequencyBand == frequencyBand) %>%
    ungroup() %>%
    group_by(ElectrodeID, TimeBin) %>%
    summarise(BinMean = mean(Mean))
}

# perform chi-square analysis
runChiSquaredTest <- function(inputData, observationType, timeType, frequencyBand, histBreaks, numTests) {
    
    # extract data
    theData  <- filterOnOffData(inputData, observationType, timeType, frequencyBand)
    
    # build histogram
    theHist  <- hist(theData$Time, breaks = histBreaks)
   
    # perform chi-square
    theChi  <- chisq.test(theHist$counts)
    
    # output for summary array
    output <- updateChiSummary(theChi, observationType, timeType, frequencyBand, numTests)
    
    return(list(summaryOutput = output, chiSquaredInput = theHist, chiSquaredOutput = theChi))
    
}

# perform post-hoc analysis on chi-squared data
runChiSquaredPostHoc <- function (histData, binNumber) {
  
  # use exact binomial test to determine whether one bin significantly differs
  # from other bins
  thisBinData <- histData$counts[binNumber]
  otherBinData <- sum(histData$counts[setdiff(seq(length(histData$counts)), binNumber)])
  
  binomResult <- binom.test(c(thisBinData, otherBinData), p = (1 / (length(histData$counts) - 1)))
  
  output <- data.frame(ActualP = binomResult$p.value, 
                       CorrectedP = ifelse(binomResult$p.value < 0.05, 
                              binomResult$p.value * length(histData$counts),
                              NA)) # Bonferroni correction
  
}

# return data from a specific frequency band and event type
filterOnOffData <- function(trimmedOnOffData, observationType, timeType, frequencyBand) {
  filteredData <- trimmedOnOffData %>%
    filter(FrequencyBand == frequencyBand &
           ObservationType == observationType &
           TrialTimeType == timeType) 
}

# update chi-squared summmary array
updateChiSummary <- function(chiOutput, observationType, timeType, frequencyBand, numTests) { 
  
  output <- data.frame(ObservationType = observationType, 
              TrialTimeType = timeType, 
              FrequencyBand = frequencyBand,
              ChiSquared = chiOutput$statistic,
              df = chiOutput$parameter, 
              RawPValue = chiOutput$p.value, 
              CorrectedPValue = ifelse(chiOutput$p.value < 0.05, chiOutput$p.value * numTests, NA),
              stringsAsFactors = FALSE,
              row.names = NULL)
  
}

# return the values from the one-way anova in a data frame
cleanAnova <- function(inputData) {
  aov.out <-aov(BinMean ~ TimeBin + Error(ElectrodeID/TimeBin), data = inputData)
  aov.selection <- summary(aov.out)$"Error: ElectrodeID:TimeBin"[[1]]
}

# perform bayesian null hypothesis testing (Masson 2011)
bayesNullHypTest <- function(anovaStats, numElectrodes, numConditions) {
  
  nIndObs  <- numElectrodes * (numConditions - 1)
  dfEffect <- anovaStats$Df[1]
  ssEffect <- anovaStats$"Sum Sq"[1]
  ssError  <- anovaStats$"Sum Sq"[2]
  
  sse1 <- ssError
  sse0 <- ssEffect + ssError
  
  deltaBIC <- (nIndObs * log(sse1 / sse0)) + (dfEffect * log(nIndObs))
  bf01 <- exp(deltaBIC / 2)
  pH0  <- bf01 / (1 + bf01)
  
  return(pH0)
}

# save anova results
saveAnovaResults <- function(anovaData, inputData, trialTimeType, thisFreqBand) {
  
  if (anovaData$"Pr(>F)"[1] < 0.05) {
    
    control <- how(within = Within(type = "free"),
                   blocks = inputData$ElectrodeID)
    trueP <- runAnovaPermutations(inputData, anovaData$"F value"[1], control, 1000)
    
    if (trueP < 0.05) {
      postHocResult <- runAnovaPostHoc(inputData)
    }
    
    anovaOutput <- data.frame(FrequencyBand = thisFreqBand,
                              TrialTimeType = trialTimeType,
                              anovaF = anovaData$"F value"[1],
                              anovaP = anovaData$"Pr(>F)"[1],
                              PermCorrectedP = trueP,
                              NullHypTestP = NA)
    postHocOutput <- data.frame(FrequencyBand = thisFreqBand,
                                TrialTimeType = trialTimeType,
                                Comparison = postHocResult$Comparison,
                                VStatistic = postHocResult$VStatistic,
                                ActualP = postHocResult$ActualP,
                                CorrectedP = postHocResult$CorrectedP)

  } else {
    
    pH0 <- bayesNullHypTest(anovaStats = anovaData, numElectrodes = nlevels(inputData$ElectrodeID), numConditions = nlevels(inputData$TimeBin))
    anovaOutput <- data.frame(FrequencyBand = thisFreqBand,
                              TrialTimeType = trialTimeType,
                              anovaF = anovaData$"F value"[1],
                              anovaP = anovaData$"Pr(>F)"[1],
                              PermCorrectedP = NA,
                              NullHypTestP = pH0)
    postHocOutput <- data.frame(Comparison = NA,
                                VStatistic = NA,
                                ActualP = NA,
                                CorrectedP = NA)
    
  }
  
  return(list(anovaOutput = anovaOutput, postHocOutput = postHocOutput))
  
}

# run posthoc analysis for anova
runAnovaPostHoc <- function(inputData) {
  
  output <- data.frame(Comparison = NA, VStatistic = NA, ActualP = NA, CorrectedP = NA)
  
  # make list of pairwise comparisons
  pairwise <- data.frame(Cond1 = c("Pre", "Pre", "Tele"),
                         Cond2 = c("Tele", "Post", "Post"))
  pairwise$Cond1 <- factor(pairwise$Cond1, levels = levels(inputData$TimeBin))
  pairwise$Cond2 <- factor(pairwise$Cond2, levels = levels(inputData$TimeBin))
  
  for (thisComparison in 1:nrow(pairwise)) {
    
    cond1 <- inputData %>% filter(TimeBin == pairwise$Cond1[thisComparison])
    cond2 <- inputData %>% filter(TimeBin == pairwise$Cond2[thisComparison])
    
    wilcoxResult <- wilcox.test(x = cond1$BinMean, y = cond2$BinMean, paired = TRUE)
    
    output[nrow(output) + 1, ] <- c(paste(pairwise$Cond1[thisComparison], "vs.", pairwise$Cond2[thisComparison]),
                                    wilcoxResult$statistic,
                                    wilcoxResult$p.value,
                                    ifelse(wilcoxResult$p.value < 0.05, 
                                           wilcoxResult$p.value * nrow(pairwise),
                                           NA))
    
  }
  output <- output %>%
    filter(is.na(Comparison) == FALSE)
  return(output)
  
}

# run permutation analysis for anova
runAnovaPermutations <- function(inputData, trueF, controlStructure, numPerms) {
  
  permF <- rep(0, numPerms)
  permOrders <- shuffleSet(nrow(inputData), numPerms, controlStructure)
  
  for (thisPerm in 1:numPerms) {
    
    tempData <- inputData
    tempData$BinMean <- tempData$BinMean[permOrders[thisPerm, ]]
    permAnova <- cleanAnova(tempData)
    permF[thisPerm] <-permAnova$"F value"[1]
    
  }
  
  permF[length(permF) + 1] <- trueF
  permF <- sort(permF, decreasing = TRUE)
  trueP <- which(permF == trueF) / length(permF)
  return(trueP)
}


# Perform repeated measures ANOVAs ----------------------------------------

# initialize output array
aovStatsOutput <- data.frame(FrequencyBand = NA, TrialTimeType = NA, anovaF = NA, anovaP = NA, PermCorrectedP = NA, NullHypTestP = NA)
aovPostHocOutput <- data.frame(FrequencyBand = NA, TrialTimeType = NA, Comparison = NA, VStatistic = NA, ActualP = NA, CorrectedP = NA)

for (thisFreqBand in 1:nlevels(binnedEpisodeData$FrequencyBand)) {
  
  print(paste0('Working on ', levels(binnedEpisodeData$FrequencyBand)[thisFreqBand]))
  print('Filtering data...')
  
  ntBinnedData <- filterRMData(binnedEpisodeData, "NT", levels(binnedEpisodeData$FrequencyBand)[thisFreqBand])
  ftBinnedData <- filterRMData(binnedEpisodeData, "FT", levels(binnedEpisodeData$FrequencyBand)[thisFreqBand])
  
  print('Running ANOVAs...')
  ntAov <- cleanAnova(ntBinnedData)
  ftAov <- cleanAnova(ftBinnedData)
  
  # If the ANOVA is significant, perform permutation analysis; otherwise,
  # perform Bayesian null hypothesis testing
  anovaResults <- saveAnovaResults(anovaData = ntAov,
                                   inputData = ntBinnedData,
                                   trialTimeType = "NT",
                                   thisFreqBand = levels(binnedEpisodeData$FrequencyBand)[thisFreqBand])
  aovStatsOutput <- rbind(aovStatsOutput, anovaResults$anovaOutput)
  if (is.na(anovaResults$postHocOutput$Comparison[1]) == FALSE) {
    aovPostHocOutput <- rbind(aovPostHocOutput, anovaResults$postHocOutput)
  }
  
  anovaResults <- saveAnovaResults(anovaData = ftAov,
                                   inputData = ftBinnedData,
                                   trialTimeType = "FT",
                                   thisFreqBand = levels(binnedEpisodeData$FrequencyBand)[thisFreqBand])
  aovStatsOutput <- rbind(aovStatsOutput, anovaResults$anovaOutput)
  if (is.na(anovaResults$postHocOutput$Comparison[1]) == FALSE) {
    aovPostHocOutput <- rbind(aovPostHocOutput, anovaResults$postHocOutput)
  }
  
   
}

# Perform chi-square on onsets and offsets --------------------------------

chiSquareOutput <- data.frame(Event = NA, TrialTimeType = NA, FrequencyBand = NA, ChiSquared = NA, df = NA, RawPValue = NA, CorrectedPValue = NA)
chiSquarePostHocOutput <- data.frame(Event = NA, TrialTimeType = NA, FrequencyBand = NA, Bin = NA, RawPValue = NA, CorrectedPValue = NA)

# hist pretty-ifys the breaks, which we don't want, so we'll set the breaks
# manually
ntBreaks <- seq(-1830, 3660, 305)
ftBreaks <- seq(-2830, 5660, 283)

# data.frame of analysis types
analysisDF <- data.frame(ObservationType = c("Onset", "Onset", "Offset", "Offset"),
                         TimeType = c("NT", "FT", "NT", "FT"),
                         HistBreaks = c("ntBreaks", "ftBreaks", "ntBreaks", "ftBreaks"))

for (thisFreqBand in 1:nlevels(trimmedOnOffData$FrequencyBand)) {
  
  for (thisAnalysis in 1:nrow(analysisDF)) {
    
    # perform chi-squared analysis for this frequency band and this analysis type
    chiOutput <- runChiSquaredTest(trimmedOnOffData, 
                                   as.character(analysisDF$ObservationType[thisAnalysis]),
                                   as.character(analysisDF$TimeType[thisAnalysis]),
                                   levels(trimmedOnOffData$FrequencyBand)[thisFreqBand],
                                   eval(as.name(as.character(analysisDF$HistBreaks[thisAnalysis]))),
                                   nlevels(trimmedOnOffData$FrequencyBand) * nrow(analysisDF))
    chiSquareOutput[nrow(chiSquareOutput) + 1, ] <- chiOutput$summaryOutput
    
    # posthoc tests if chi-squared is significant
    if (is.na(chiOutput$summaryOutput$CorrectedPValue) == FALSE &
          chiOutput$summaryOutput$CorrectedPValue < 0.05) {
      
      histData <- chiOutput$chiSquaredInput
      
      for (thisBin in 1:length(histData$counts)) {
        
        chiPostHoc <- runChiSquaredPostHoc(histData, thisBin)
        chiSquarePostHocOutput[nrow(chiSquarePostHocOutput) + 1, ] <- c(as.character(analysisDF$ObservationType[thisAnalysis]),
                                                      as.character(analysisDF$TimeType[thisAnalysis]),
                                                      levels(trimmedOnOffData$FrequencyBand)[thisFreqBand],
                                                      thisBin,
                                                      chiPostHoc$ActualP,
                                                      chiPostHoc$CorrectedP)
      }
      
    }
    
  }
  
}

aovStatsOutput         <- aovStatsOutput %>% filter(is.na(FrequencyBand) == FALSE)
aovPostHocOutput       <- aovPostHocOutput %>% filter(is.na(FrequencyBand) == FALSE)
chiSquareOutput        <- chiSquareOutput %>% filter(is.na(FrequencyBand) == FALSE)
chiSquarePostHocOutput <- chiSquarePostHocOutput %>% filter(is.na(FrequencyBand) == FALSE)

save(file = 'Rda/allStats.Rda', list = c("aovStatsOutput", "aovPostHocOutput", "chiSquareOutput", "chiSquarePostHocOutput"))