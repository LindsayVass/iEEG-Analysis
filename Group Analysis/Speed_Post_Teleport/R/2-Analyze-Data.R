# Script Name:  2-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         1 October 2015
# Purpose:      This script will get the mean speed for each electrode for pre 
#               and post-teleport epochs

library(dplyr)
library(reshape2)
library(broom)

load('Rda/allRawData.Rda')

speedData <- speedData %>%
  mutate(ElectrodeID = paste(Subject, Session, Depth, sep = "_")) %>%
  group_by(ElectrodeID, Interval) %>%
  filter(TrialNumber != 64)  # no speed data for post epoch of last trial

summarySpeedData <- speedData %>%
  summarise(MeanSpeed = mean(Speed),
            SEMSpeed = sd(Speed) / sqrt(n()))

groupSummarySpeedData <- speedData %>%
  ungroup() %>%
  group_by(Interval) %>%
  summarise(MeanSpeed = mean(Speed * 1000), # multiply by 1000 to convert to units/s
            SEMSpeed = sd(Speed * 1000) / sqrt(n()))

subjectSummarySpeedData <- speedData %>%
  ungroup() %>%
  group_by(Interval, Subject) %>%
  summarise(MeanSpeed = mean(Speed * 1000), # multiply by 1000 to convert to units/s
            SEMSpeed = sd(Speed * 1000) / sqrt(n()))

UCDMC1315Summary <- speedData %>%
  filter(Subject != "UCDMC14") %>%
  ungroup() %>%
  group_by(Interval) %>%
  summarise(MeanSpeed = mean(Speed * 1000), # multiply by 1000 to convert to units/s
            SEMSpeed = sd(Speed * 1000) / sqrt(n()))

summaryStats <- speedData %>%
  dcast(ElectrodeID + TrialNumber ~ Interval, value.var = "Speed") %>%
  group_by(ElectrodeID) %>%
  do(TTest = t.test(.$Post, .$Pre, paired = TRUE))
summaryStats <- tidy(summaryStats, TTest)

save(file = 'Rda/allAnalyzedData.Rda', list = c('speedData', 'summarySpeedData', 'summaryStats'))
