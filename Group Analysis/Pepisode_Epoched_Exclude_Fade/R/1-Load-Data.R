# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         23 November 2015
# Purpose:      This script will load the csv files output by the Matlab
#               function "Run_Teleporter_Epoched_Exclude_Fade.m"


# Functions ---------------------------------------------------------------

makeSessionInfo <- function() {
  
  UCDMC13 <- data.frame(Subject = "UCDMC13", Session = "TeleporterA")
  UCDMC14 <- data.frame(Subject = "UCDMC14", Session = c("TeleporterA", "TeleporterB"))
  UCDMC15 <- data.frame(Subject = "UCDMC15", Session = c("TeleporterA", "TeleporterB"))
  
  sessionInfo <- rbind(UCDMC13, UCDMC14, UCDMC15)
  
  return(sessionInfo)
  
}

loadSubjectData <- function(subject, session, csvFileStem) {
  
  inputFile <- paste0('csv/', subject, '_', session, '_', csvFileStem)
  tempData  <- read.csv(inputFile, header = TRUE)
  return(tempData)
  
}


# Script ------------------------------------------------------------------

# Make a data frame that contains which subjects and sessions we have data for
sessionInfo <- makeSessionInfo()

# Load the data from each subject's session(s) and combine into one data frame
allData <- data.frame()
for (thisRow in 1:nrow(sessionInfo)) {
  
  subjectData <- loadSubjectData(sessionInfo$Subject[thisRow], 
                                 sessionInfo$Session[thisRow],
                                 'epoched_pepisode_equal_time_bins.csv')
    
  allData <- rbind(allData, subjectData)
}
dir.create('Rda')
save(file = 'Rda/allRawData.Rda', list = 'allData')
