# Script Name:  1b-Load-Data-ExcludeTrial64.R
# Author:       Lindsay Vass
# Date:         24 November 2015
# Purpose:      This script will load the csv files output by the Matlab
#               function "Run_Teleporter_Epoched_Pepisode_Equal_Time_Bins.m"

library(R.matlab)

# Functions ---------------------------------------------------------------

makeSessionInfo <- function() {
  
  UCDMC13 <- data.frame(Subject = "UCDMC13", 
                        Session = "TeleporterA", 
                        Electrode = c("LAD", "LHD"))
  UCDMC14 <- data.frame(Subject = "UCDMC14", 
                        Session = c(rep("TeleporterA", 4), rep("TeleporterB", 4)), 
                        Electrode = c("LAD", "LHD", "RAD", "RHD"))
  UCDMC15 <- data.frame(Subject = "UCDMC15", 
                        Session = c(rep("TeleporterA", 4), rep("TeleporterB", 4)), 
                        Electrode = c("LAD", "LHD", "RAD", "RHD"))
  
  sessionInfo <- rbind(UCDMC13, UCDMC14, UCDMC15)
  
  return(sessionInfo)
  
}

makeSessionInfoNoElectrode <- function() {
  
  UCDMC13 <- data.frame(Subject = "UCDMC13", 
                        Session = "TeleporterA")
  UCDMC14 <- data.frame(Subject = "UCDMC14", 
                        Session = c("TeleporterA", "TeleporterB"))
  UCDMC15 <- data.frame(Subject = "UCDMC15", 
                        Session = c("TeleporterA", "TeleporterB"))
  
  sessionInfo <- rbind(UCDMC13, UCDMC14, UCDMC15)
  
  return(sessionInfo)
  
}

loadSubjectData <- function(subject, session, csvFileStem) {
  
  inputFile <- paste0('csv/', subject, '_', session, '_', csvFileStem)
  tempData  <- read.csv(inputFile, header = TRUE)
  return(tempData)
  
}

loadTrialNumbers <- function(dir, sessionInfo, fileStem, variable, EDFoption) {
  
  realTrialNumbers <- data.frame()
  
  for (thisRow in 1:nrow(sessionInfo)) {
    
    subject    <- sessionInfo$Subject[thisRow]
    teleporter <- sessionInfo$Session[thisRow]
    electrode  <- sessionInfo$Electrode[thisRow]
    
    if (EDFoption == TRUE) {
      
      if (is.na(sessionInfo$EDF[thisRow]) == FALSE) {
        matPath <- paste0(dir,
                          subject,
                          '/Mat Files/',
                          subject,
                          '_',
                          teleporter,
                          '_',
                          electrode,
                          fileStem,
                          '_',
                          sessionInfo$EDF[thisRow],
                          '.mat')
      } else{
        matPath <- paste0(dir,
                          subject,
                          '/Mat Files/',
                          subject,
                          '_',
                          teleporter,
                          '_',
                          electrode,
                          fileStem,
                          '.mat')
      }
    } else {
      matPath <- paste0(dir,
                        subject,
                        '/Mat Files/',
                        subject,
                        '_',
                        teleporter,
                        '_',
                        electrode,
                        fileStem,
                        '.mat')
    }
    
    
    if (file.exists(matPath) == FALSE) {
      cat(paste0('File not found: \n', matPath, '\n'))
      next
    }
    
    realTrials <- readMat(matPath)[[variable]]
    
    thisData <- data.frame(Subject  = subject,
                           Session   = teleporter,
                           Electrode = electrode,
                           Trial     = realTrials)
    realTrialNumbers <- rbind(realTrialNumbers, thisData)
  }
  
  return(realTrialNumbers)
}

# Script ------------------------------------------------------------------

# Make a data frame that contains which subjects and sessions we have data for
sessionInfo <- makeSessionInfo()
sessionInfoNoElec <- makeSessionInfoNoElectrode()

# Load the data from each subject's session(s) and combine into one data frame
allData <- data.frame()
for (thisRow in 1:nrow(sessionInfoNoElec)) {
  
  subjectData <- loadSubjectData(sessionInfoNoElec$Subject[thisRow], 
                                 sessionInfoNoElec$Session[thisRow],
                                 'epoched_pepisode_equal_time_bins.csv')
    
  allData <- rbind(allData, subjectData)
}

realTeleporterTrialNumbers <- loadTrialNumbers(dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/',
                                               sessionInfo = sessionInfo,
                                               fileStem = '_noSpikes_noWaves_goodEpochs',
                                               variable = 'goodEpochs',
                                               EDFoption = FALSE)

save(file = 'Rda/allRawData_ExcludeTrial64.Rda', list = c('allData', 'realTeleporterTrialNumbers', 'sessionInfo'))
