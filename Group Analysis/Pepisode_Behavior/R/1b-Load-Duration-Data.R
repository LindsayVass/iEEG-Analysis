# Script Name:  1b-Load-Duration-Data.R
# Author:       Lindsay Vass
# Date:         30 November 2015
# Purpose:      This script will load the csv files output by the Matlab
#               script "runCalculateTrialDuration.m", which contains
#               the duration for each trial.

library(plyr)
library(dplyr)
library(tidyr)

# Load latency csv files --------------------------------------------------

# latency data is stored in csv files
inputFiles  <- dir('csv/', pattern = "*Trial_Duration.csv")
latencyData <- data.frame()

for (thisFile in 1:length(inputFiles)) {
  
  tempData <- read.csv(paste0('csv/', inputFiles[thisFile]), header = TRUE)
  
  # add information about subject, session, and electrode to data frame based on
  # file name
  fileNameParts <- strsplit(inputFiles[thisFile], split = "_")
  
  latencyData        <- rbind(latencyData, tempData)
  
}

remove(tempData)

# Save data ---------------------------------------------------------------

save(file = 'Rda/allRawData_Duration.Rda', list = c('latencyData'))
