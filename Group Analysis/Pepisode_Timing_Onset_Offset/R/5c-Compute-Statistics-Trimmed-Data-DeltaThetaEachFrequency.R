# Script Name:  5c-Compute-Statistics-Trimmed-Data.R
# Author:       Lindsay Vass
# Date:         25 June 2015
# Purpose:      This script will calculate statistics for the group data plotted
#               in 4c-Plot-Trimmed-Data-DeltaThetaEachFrequency.R. It will
#               use the chi-square goodness of fit test to determine whether the
#               distributions of onsets and offsets differ from a uniform 
#               distribution. Posthoc binomial tests will identify specific bins
#               whose frequency significantly differs from the expected frequency.
#               Finally, it will compute the pairwise similarity of the 
#               histograms between frequencies using a chi-squared test. This
#               similarity will be visualized in matrix form.


library(dplyr)
library(permute)
load('Rda/allTrimmedData_DeltaThetaByFrequency.Rda')

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
runChiSquaredTest <- function(inputData, observationType, timeType, frequency, histBreaks, numTests) {
    
    # extract data
    theData  <- filterOnOffData(inputData, observationType, timeType, frequency)
    
    # build histogram
    theHist  <- hist(theData$Time, breaks = histBreaks)
   
    # perform chi-square
    theChi  <- chisq.test(theHist$counts)
    
    # output for summary array
    output <- updateChiSummary(theChi, observationType, timeType, frequency, numTests)
    
    return(list(summaryOutput = output, chiSquaredInput = theHist, chiSquaredOutput = theChi))
    
}

# perform chi-square analysis collapsed over frequency
runCollapsedChiSquaredTest <- function(inputData, observationType, timeType, histBreaks, numTests) {
  theData <- inputData %>%
    filter(ObservationType == observationType &
             TrialTimeType == timeType)
  theHist <- hist(theData$Time, breaks = histBreaks)
  theChi <- chisq.test(theHist$counts)
  output <- updateChiSummary(theChi, observationType, timeType, "Delta/Theta", numTests)
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

# return data from a specific frequency and event type
filterOnOffData <- function(trimmedOnOffData, observationType, timeType, frequency) {
  filteredData <- trimmedOnOffData %>%
    filter(Frequency == frequency &
           ObservationType == observationType &
           TrialTimeType == timeType) 
}

# update chi-squared summmary array
updateChiSummary <- function(chiOutput, observationType, timeType, frequency, numTests) { 
  
  output <- data.frame(ObservationType = observationType, 
              TrialTimeType = timeType, 
              Frequency = frequency,
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

# pairwise chi-squared test between frequencies
pairwiseChiSquared <- function(inputData, observationType, timeType, frequency1, frequency2, histBreaks) {
  
  # extract data
  data1 <- filterOnOffData(inputData, observationType, timeType, frequency1)
  data2 <- filterOnOffData(inputData, observationType, timeType, frequency2)
  
  # build histograms
  hist1 <- hist(data1$Time, breaks = histBreaks)
  hist2 <- hist(data2$Time, breaks = histBreaks)
  
  # chi squared test
  chiResult <- chisq.test(hist1$counts, hist2$counts)
  
  return(list(ChiSq = chiResult$statistic,
              PValue = chiResult$p.value))
  
}

# Perform chi-square on onsets and offsets --------------------------------

chiSquareOutput <- data.frame(Event = NA, TrialTimeType = NA, Frequency = NA, ChiSquared = NA, df = NA, RawPValue = NA, CorrectedPValue = NA)
chiSquarePostHocOutput <- data.frame(Event = NA, TrialTimeType = NA, Frequency = NA, Bin = NA, RawPValue = NA, CorrectedPValue = NA)

# hist pretty-ifys the breaks, which we don't want, so we'll set the breaks
# manually
ntBreaks <- seq(-1830, 3660, 305)
ftBreaks <- seq(-2830, 5660, 283)

# data.frame of analysis types
analysisDF <- data.frame(ObservationType = c("Onset", "Onset", "Offset", "Offset"),
                         TimeType = c("NT", "FT", "NT", "FT"),
                         HistBreaks = c("ntBreaks", "ftBreaks", "ntBreaks", "ftBreaks"))

for (thisFreq in 1:nlevels(trimmedOnOffData$Frequency)) {
  
  for (thisAnalysis in 1:nrow(analysisDF)) {
    
    # perform chi-squared analysis for this frequency band and this analysis type
    chiOutput <- runChiSquaredTest(trimmedOnOffData, 
                                   as.character(analysisDF$ObservationType[thisAnalysis]),
                                   as.character(analysisDF$TimeType[thisAnalysis]),
                                   levels(trimmedOnOffData$Frequency)[thisFreq],
                                   eval(as.name(as.character(analysisDF$HistBreaks[thisAnalysis]))),
                                   nlevels(trimmedOnOffData$Frequency) * nrow(analysisDF))
    chiSquareOutput[nrow(chiSquareOutput) + 1, ] <- chiOutput$summaryOutput
    
    # posthoc tests if chi-squared is significant
    if (is.na(chiOutput$summaryOutput$CorrectedPValue) == FALSE &
          chiOutput$summaryOutput$CorrectedPValue < 0.05) {
      
      histData <- chiOutput$chiSquaredInput
      
      for (thisBin in 1:length(histData$counts)) {
        
        chiPostHoc <- runChiSquaredPostHoc(histData, thisBin)
        chiSquarePostHocOutput[nrow(chiSquarePostHocOutput) + 1, ] <- c(as.character(analysisDF$ObservationType[thisAnalysis]),
                                                      as.character(analysisDF$TimeType[thisAnalysis]),
                                                      levels(trimmedOnOffData$Frequency)[thisFreq],
                                                      thisBin,
                                                      chiPostHoc$ActualP,
                                                      chiPostHoc$CorrectedP)
      }
      
    }
    
  }
  
}

chiSquareOutput        <- chiSquareOutput %>% filter(is.na(Frequency) == FALSE)
chiSquarePostHocOutput <- chiSquarePostHocOutput %>% filter(is.na(Frequency) == FALSE)


# Perform chi-square on onsets and offsets collapsed over frequency -------

chiSquareCollapsedOutput <- data.frame()
chiSquareCollapsedPostHocOutput <- data.frame()

for (thisAnalysis in 1:nrow(analysisDF)) {
  chiOutput <- runCollapsedChiSquaredTest(trimmedOnOffData,
                                          as.character(analysisDF$ObservationType[thisAnalysis]),
                                          as.character(analysisDF$TimeType[thisAnalysis]),
                                          eval(as.name(as.character(analysisDF$HistBreaks[thisAnalysis]))),
                                          nrow(analysisDF))
  chiSquareCollapsedOutput <- rbind(chiSquareCollapsedOutput, chiOutput$summaryOutput)
  
  if (is.na(chiOutput$summaryOutput$CorrectedPValue) == FALSE &
        chiOutput$summaryOutput$CorrectedPValue < 0.05) {
    
    histData <- chiOutput$chiSquaredInput
    
    for (thisBin in 1:length(histData$counts)) {
      
      chiPostHoc <- runChiSquaredPostHoc(histData, thisBin)
      chiSquareCollapsedPostHocOutput <- rbind(chiSquareCollapsedPostHocOutput, data.frame(ObservationType = as.character(analysisDF$ObservationType[thisAnalysis]),
                                                                                  TimeType = as.character(analysisDF$TimeType[thisAnalysis]),
                                                                                  Frequency = "Delta/Theta",
                                                                                  BinNumber = thisBin,
                                                                                  ActualP = chiPostHoc$ActualP,
                                                                                  CorrectedP = chiPostHoc$CorrectedP)) 
    }
    
  }
}



# Pairwise similarity between frequencies ---------------------------------

pairwiseData <- data.frame()

for (firstFreq in 1:(nlevels(trimmedOnOffData$Frequency) - 1)) {
  
  for (secondFreq in (firstFreq + 1):nlevels(trimmedOnOffData$Frequency)) {
    
    for (thisAnalysis in 1:nrow(analysisDF)) {
      
      thisChi <- pairwiseChiSquared(trimmedOnOffData,
                                    as.character(analysisDF$ObservationType[thisAnalysis]),
                                    as.character(analysisDF$TimeType[thisAnalysis]),
                                    levels(trimmedOnOffData$Frequency)[firstFreq],
                                    levels(trimmedOnOffData$Frequency)[secondFreq],
                                    eval(as.name(as.character(analysisDF$HistBreaks[thisAnalysis]))))
      
      thisOutput <- data.frame(Freq1 = c(substr(levels(trimmedOnOffData$Frequency)[firstFreq], 1, 4),
                                         substr(levels(trimmedOnOffData$Frequency)[secondFreq], 1, 4)),
                               Freq2 = c(substr(levels(trimmedOnOffData$Frequency)[secondFreq], 1, 4),
                                         substr(levels(trimmedOnOffData$Frequency)[firstFreq], 1, 4)),
                               ChiSq = thisChi$ChiSq,
                               PValue = thisChi$PValue,
                               ObservationType = analysisDF$ObservationType[thisAnalysis],
                               TimeType = analysisDF$TimeType[thisAnalysis])
      pairwiseData <- rbind(pairwiseData, thisOutput)
      
    }
  }
}


# Plot dissimilarity matrices ---------------------------------------------

dir.create('Figures/TimingSimilarityMatrices/')

goodHists <- c(substr(levels(trimmedOnOffData$Frequency)[3:nlevels(trimmedOnOffData$Frequency)], 1, 4))

plotData <- pairwiseData %>%
  filter(Freq1 %in% goodHists &
           Freq2 %in% goodHists)

for (thisPlot in 1:nrow(analysisDF)) {
  
  p <- plotData %>%
    filter(ObservationType == analysisDF$ObservationType[thisPlot] &
             TimeType == analysisDF$TimeType[thisPlot]) %>%
    ggplot(aes(x = Freq1, y = Freq2)) +
    geom_tile(aes(fill = ChiSq)) +
    scale_fill_gradient(low = "red", high = "blue") + 
    ggtitle(paste(analysisDF$TimeType[thisPlot], analysisDF$ObservationType[thisPlot], 'Timing Similarity'))
  ggsave(paste0('Figures/TimingSimilarityMatrices/', analysisDF$TimeType[thisPlot], '_', analysisDF$ObservationType[thisPlot], 'TimingSimilarity.png'),
         width = 12, height = 12)
  
}

save(file = 'Rda/allStats-DeltaThetaEachFrequency.Rda', list = c("chiSquareOutput", "chiSquarePostHocOutput", "pairwiseData", "chiSquareCollapsedOutput", "chiSquareCollapsedPostHocOutput"))