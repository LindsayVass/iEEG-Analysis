# Script Name:  3b-Analyze-Duration-Data.R
# Author:       Lindsay Vass
# Date:         30 November 2015
# Purpose:      This script will calculate summary statistics for the duration data.

library(dplyr)

load('Rda/allRawData_Duration.Rda')

durationStats <- latencyData %>%
  group_by(SubjectID, Session) %>%
  summarise(MeanDuration = mean(Duration),
            SEMDuration = sd(Duration) / sqrt(n()))