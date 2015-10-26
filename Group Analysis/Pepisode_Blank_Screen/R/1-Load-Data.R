# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         23 October 2015
# Purpose:      This script will load the pepisode data output by m/Run_Blank_Screen_Pepisode.m
#               as well as pepisode data measured during nav/tele previously.

library(dplyr)

# load blank screen data
csvPath <- 'csv/output_23-Oct-2015.csv'
blankData <- read.csv(csvPath, header = TRUE, sep = ",")

blankData <- blankData %>%
  mutate(ElectrodeID = paste(Subject, Teleporter, Electrode, sep = "_"))

# load previously measured nav/tele data
load('../Pepisode_Sustainedness/Rda/allCleanData.Rda')

save(file = 'Rda/allRawData.Rda', list = c('blankData', 'allPepisode'))