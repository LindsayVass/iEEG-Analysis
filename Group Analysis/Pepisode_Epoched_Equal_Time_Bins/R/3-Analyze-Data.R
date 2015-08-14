# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         3 June 2015
# Purpose:      This script will analyze the data frame produced by 
#               "2-Clean-Data.R". For each electrode and frequency band, it will 
#               use a Wilcoxon signed-rank test to determine whether pepisode 
#               significantly differs between timepoints. 

library(dplyr)
library(reshape2)

load('Rda/allCleanData.Rda')

# Wilcoxon signed-rank at each electrode and frequency band ---------------

wilcoxonResults <- data.frame(ElectrodeID = NA,
                              FrequencyBand = NA,
                              Contrast = NA,
                              P = NA)
thisRow <- 1

for (thisElectrode in 1:nlevels(cleanData$ElectrodeID)){
  
  # Get all data for this electrode
  electrodeData <- cleanData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElectrode]) %>%
    select(-ElectrodeID) 
  
  for (thisFreqBand in 1:nlevels(cleanData$FrequencyBand)) {
    
    # Filter data for this frequency band
    frequencyData <- electrodeData %>%
      filter(FrequencyBand == levels(FrequencyBand)[thisFreqBand]) %>%
      select(-FrequencyBand) %>%
      group_by(TrialNumber, TimePoint) %>%
      summarise(Pepisode = mean(Pepisode))
    
    # Cast to wide format
    frequencyDataWide <- frequencyData %>%
      dcast(TrialNumber ~ TimePoint, value.var = "Pepisode")
    
    # Wilcoxon signed-rank tests
    preGtTele  <- wilcox.test(frequencyDataWide$Pre1, frequencyDataWide$Tele, alternative = "greater", paired=TRUE)
    preLtTele  <- wilcox.test(frequencyDataWide$Pre1, frequencyDataWide$Tele, alternative = "less", paired=TRUE)
    teleGtPost <- wilcox.test(frequencyDataWide$Tele, frequencyDataWide$Post1, alternative = "greater", paired = TRUE)
    teleLtPost <- wilcox.test(frequencyDataWide$Tele, frequencyDataWide$Post1, alternative = "less", paired = TRUE)
    
    # Add to summary dataframe
    wilcoxonResults[thisRow, ]     <- c(levels(cleanData$ElectrodeID)[thisElectrode], 
                                        levels(cleanData$FrequencyBand)[thisFreqBand], 
                                        "Pre > Tele", 
                                        preGtTele$p.value)
    wilcoxonResults[thisRow + 1, ] <- c(levels(cleanData$ElectrodeID)[thisElectrode], 
                                        levels(cleanData$FrequencyBand)[thisFreqBand], 
                                        "Pre < Tele", 
                                        preLtTele$p.value)
    wilcoxonResults[thisRow + 2, ] <- c(levels(cleanData$ElectrodeID)[thisElectrode], 
                                        levels(cleanData$FrequencyBand)[thisFreqBand], 
                                        "Tele > Post", 
                                        teleGtPost$p.value)
    wilcoxonResults[thisRow + 3, ] <- c(levels(cleanData$ElectrodeID)[thisElectrode], 
                                        levels(cleanData$FrequencyBand)[thisFreqBand], 
                                        "Tele < Post", 
                                        teleLtPost$p.value)
    thisRow <- thisRow + 4
    
  } # end thisFreqBand
  
} # end thisElectrode


# Determine frequency of significant electrodes ---------------------------

# Reorder frequency band levels
freqBandOrder <- c("Delta-Theta", "Alpha", "Beta", "Gamma")
wilcoxonResults$FrequencyBand <- factor(wilcoxonResults$FrequencyBand, levels = freqBandOrder)

# Get the counts of the significant results
wilcoxonSigResults <- wilcoxonResults %>%
  group_by(FrequencyBand, Contrast) %>%
  filter(P < 0.05) %>%
  summarise(Count = n())

# Fill in missing frequency band/contrasts with zeros
allContrasts <- wilcoxonResults %>%
  group_by(FrequencyBand, Contrast) %>%
  summarise(Count = 0) %>%
  anti_join(wilcoxonSigResults, by = c("FrequencyBand", "Contrast"))
wilcoxonSigResults <- rbind(wilcoxonSigResults, allContrasts)


# Perform binomial test to determine whether counts are higher than expected by chance
binomP <- wilcoxonSigResults %>%
  rowwise() %>%
  do(BinomialTestP = binom.test(.$Count, 
                                nlevels(cleanData$ElectrodeID), 
                                p = 0.05, 
                                alternative = "greater")$p.value)
wilcoxonSigResults <- cbind(wilcoxonSigResults, binomP) 

# Make a dataframe for the significance markers
sigMarkers <- wilcoxonSigResults %>%
  filter(BinomialTestP < 0.05) %>%
  mutate(Count = Count + 0.2) 

save(file = 'Rda/allAnalyzedData.Rda', list = c('wilcoxonSigResults', 'wilcoxonResults', 'sigMarkers'))