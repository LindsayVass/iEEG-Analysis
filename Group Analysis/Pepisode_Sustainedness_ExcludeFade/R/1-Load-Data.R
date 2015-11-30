# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         13 July 2015
# Purpose:      This script will load the csv files output by the Matlab
#               function "Run_Teleporter_Epoched_Pepisode_Equal_Time_Bins.m". 
#               This data is from when subjects are actively navigating before
#               entering the teleporter. It will also load the matching csv
#               files from ../Pepisode_Epoched_ExcludeFade/, which is data
#               from the teleportation epoch. These files contain the pepisode
#               values for each time bin.
#               

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

makeSessionInfoWithEDF <- function() {
  
  UCDMC13 <- data.frame(Subject = "UCDMC13", 
                        Session = "TeleporterA", 
                        Electrode = c("LAD", "LHD"), 
                        EDF = NA)
  UCDMC14 <- data.frame(Subject = "UCDMC14", 
                        Session = c(rep("TeleporterA", 4), rep("TeleporterB", 4)), 
                        Electrode = c("LAD", "LHD", "RAD", "RHD"), 
                        EDF = NA)
  UCDMC15 <- data.frame(Subject = "UCDMC15", 
                        Session = c(rep("TeleporterA", 8), rep("TeleporterB", 8)), 
                        Electrode = rep(c("LAD", "LHD", "RAD", "RHD"), 4), 
                        EDF = c(rep("EDF1", 4), rep("EDF2", 4)))
  
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

loadSubjectData <- function(subject, session, path, csvFileStem) {
  
  inputFile <- paste0(path, '/csv/', subject, '_', session, '_', csvFileStem)
  tempData  <- read.csv(inputFile, header = TRUE)
  return(tempData)
  
}

loadAllSubjectData <- function(dir, sessionInfo, fileStem) {
  output <- data.frame()
  for (thisRow in 1:nrow(sessionInfo)) {
    subjectData <- loadSubjectData(sessionInfo$Subject[thisRow],
                                   sessionInfo$Session[thisRow],
                                   dir,
                                   fileStem)
    output <- rbind(output, subjectData)
  }
  return(output)
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
sessionInfo    <- makeSessionInfo()
sessionInfoEdf <- makeSessionInfoWithEDF()
sessionInfoNoElec <- makeSessionInfoNoElectrode()

# Load the data from each subject's session(s) and combine into one data frame
allNavigationPepisodeData <- loadAllSubjectData(getwd(), 
                                                sessionInfoNoElec,
                                                'epoched_navigation_pepisode_equal_time_bins.csv')
allTeleporterPepisodeData <- loadAllSubjectData('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Epoched_Exclude_Fade',
                                                sessionInfoNoElec,
                                                'epoched_pepisode_equal_time_bins.csv')



# Load the real trial numbers
realNavigationTrialNumbers <- loadTrialNumbers(dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/',
                                               sessionInfo = sessionInfoEdf,
                                               fileStem = '_navigation_goodEpochs',
                                               variable = 'goodNavEpochs',
                                               EDFoption = TRUE)
realTeleporterTrialNumbers <- loadTrialNumbers(dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/',
                                               sessionInfo = sessionInfo,
                                               fileStem = '_noSpikes_noWaves_goodEpochs',
                                               variable = 'goodEpochs',
                                               EDFoption = FALSE)

save(file = 'Rda/allRawData.Rda', list = c('realNavigationTrialNumbers', 'realTeleporterTrialNumbers', 'allNavigationPepisodeData', 'allTeleporterPepisodeData', 'sessionInfo'))

