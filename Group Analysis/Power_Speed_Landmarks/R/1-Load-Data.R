# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         7 August 2015
# Purpose:      This script will load the comma separated text files output by 
#               /m/PowerSpeedLandmarksAnalysis.m

txtDir <- "/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power_Speed_Landmarks/mat/"

fileList <- list.files(txtDir)
allRawData <- data.frame()
for (thisFile in 1:length(fileList)) {
  thisData <- read.delim(paste0(txtDir, fileList[thisFile]), header = TRUE, sep = ',')
  allRawData <- rbind(allRawData, thisData)
}

save(file = 'Rda/allRawData.Rda', list = 'allRawData')