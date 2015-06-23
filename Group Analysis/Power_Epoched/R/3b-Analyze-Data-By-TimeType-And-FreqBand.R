# Script Name:  3b-Analyze-Data-By-TimeType-And-FreqBand.R
# Author:       Lindsay Vass
# Date:         22 June 2015
# Purpose:      This script will test whether z-scored log(power) differs 
#               between delta and theta for each time type (NT/FT).


library(reshape2)
library(ggplot2)
library(dplyr)

analysisData <- cleanData %>%
  filter(TimePoint == "Tele") %>%
  filter(FrequencyBand == "Delta" | FrequencyBand == "Theta") %>%
  group_by(FrequencyBand, ElectrodeID, TrialNumber, TrialTimeType) %>%
  summarise(Power = mean(Power)) %>%
  select(ElectrodeID, TrialNumber, TrialTimeType, FrequencyBand, Power) 
ttestData <- analysisData %>%
  dcast(ElectrodeID + TrialNumber ~ FrequencyBand + TrialTimeType, value.var = "Power")

t.test(ttestData$Delta_NT, ttestData$Theta_NT, paired = TRUE)
t.test(ttestData$Delta_FT, ttestData$Theta_FT, paired = TRUE)

plotData <- analysisData %>%
  ungroup() %>%
  group_by(FrequencyBand, TrialTimeType) %>%
  summarise(Power = mean(Power)) 
qplot(data = plotData, x = FrequencyBand, y = Power, facets = ~TrialTimeType)