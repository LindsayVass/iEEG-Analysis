# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         5 June 2015
# Purpose:      This script will analyze the clean data produced by 
#               "2-Clean-Data.R". It will first perform a within-electrode
#               analysis to identify when episodes occur during NT or FT trials
#               (NT and FT are analyzed separately since trial structure 
#               differs). It will then combine data across electrodes to look
#               at when episodes occur across the population.

library(dplyr)
library(reshape2)


# Within-electrode analysis -----------------------------------------------

frequencyNt <- cleanCharData %>%
  group_by(TrialTimeType) %>%
  inner_join(cleanNtData) %>%
  group_by(ObservationType) %>%
  melt(c("ElectrodeID", 
         "ObservationID", 
         "TrialNumber", 
         "TrialTimeType", 
         "FrequencyBand", 
         "ObservationType"),
       variable.name = "Time",
       value.name = "Value")

frequencyFt <- cleanCharData %>%
  group_by(TrialTimeType) %>%
  inner_join(cleanFtData) %>%
  group_by(ObservationType) %>%
  melt(c("ElectrodeID", 
         "ObservationID", 
         "TrialNumber", 
         "TrialTimeType", 
         "FrequencyBand", 
         "ObservationType"),
       variable.name = "Time",
       value.name = "Value")

# R added X to positive timepoints and X. to negative timepoints, so fix that
# now that they're no longer variable names
frequencyData <- rbind(frequencyNt, frequencyFt)
frequencyData$Time <- sub("X0", "0", frequencyData$Time)
frequencyData$Time <- sub("X.", "-", frequencyData$Time)
frequencyData$Time <- sub("X", "", frequencyData$Time)
frequencyData$Time <- as.numeric(frequencyData$Time)

# Get the mean and SEM for each time point
summaryData <- frequencyData %>%
  group_by(ElectrodeID, FrequencyBand, TrialTimeType, Time, ObservationType) %>%
  summarise(Mean = mean(Value), SEM = sd(Value) / sqrt(n()))

