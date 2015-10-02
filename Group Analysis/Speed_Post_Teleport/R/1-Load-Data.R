# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         1 October 2015
# Purpose:      This script will load the csv files output by ../m/getMeanSpeedPTB.m
#               and getMeanSpeedUnityPulses.m

fileList <- dir('csv')

speedData <- vector(mode = "list", length = length(fileList))

for (i in 1:length(fileList)) {
  speedData[[i]] <- read.csv(paste0('csv/', fileList[i]))
}

speedData <- data.table::rbindlist(speedData)

save(file = 'Rda/allRawData.Rda', list = 'speedData')