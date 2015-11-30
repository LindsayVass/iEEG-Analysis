# Load matlab .mat files containing the intervals in seconds between the end of 
# the navigation epoch and the start of the teleportation epoch

library(R.matlab)
library(dplyr)

patients <- c('UCDMC13', 'UCDMC14')
patientDir <- '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/'

intervalData <- data.frame(Patient = NULL,
                           Session = NULL,
                           Electrode = NULL,
                           Interval = NULL)

# load UCDMC13 & 14
for (thisPatient in 1:length(patients)) {
  matPath <- paste0(patientDir, patients[thisPatient], '/Mat Files/')
  fileList <- dir(matPath)
  fileInds <- grep('_intervals.mat', fileList)
  for (thisFile in 1:length(fileInds)) {
    # load data
    thisFile <- fileList[fileInds[thisFile]]
    thisPath <- paste0(matPath, thisFile)
    data <- readMat(thisPath)
    data <- data$navTeleIntervalSecs
    
    # put into df
    fileParts <- strsplit(thisFile, '_')
    thisData <- data.frame(Patient = fileParts[[1]][1],
                           Session = fileParts[[1]][2],
                           Electrode = fileParts[[1]][3],
                           Interval = t(data)) %>%
      filter(Interval != "NaN")
    intervalData <- rbind(intervalData, thisData)
  }
}

# load UCDMC15
matPath <- paste0(patientDir, 'UCDMC15/Mat Files/')
fileList <- dir(matPath)
fileInds <- grep('_intervals_EDF1.mat', fileList)
for (thisFile in 1:length(fileInds)) {
  # load data
  thisFile <- fileList[fileInds[thisFile]]
  thisPath <- paste0(matPath, thisFile)
  data <- readMat(thisPath)
  data <- data$navTeleIntervalSecs
  
  # put into df
  fileParts <- strsplit(thisFile, '_')
  thisData <- data.frame(Patient = fileParts[[1]][1],
                         Session = fileParts[[1]][2],
                         Electrode = fileParts[[1]][3],
                         Interval = t(data)) %>%
    filter(Interval != "NaN")
  intervalData <- rbind(intervalData, thisData)
}

fileInds <- grep('_intervals_EDF2.mat', fileList)
for (thisFile in 1:length(fileInds)) {
  # load data
  thisFile <- fileList[fileInds[thisFile]]
  thisPath <- paste0(matPath, thisFile)
  data <- readMat(thisPath)
  data <- data$navTeleIntervalSecs
  
  # put into df
  fileParts <- strsplit(thisFile, '_')
  thisData <- data.frame(Patient = fileParts[[1]][1],
                         Session = fileParts[[1]][2],
                         Electrode = fileParts[[1]][3],
                         Interval = t(data)) %>%
    filter(Interval != "NaN")
  intervalData <- rbind(intervalData, thisData)
}

# summary stats
intervalSummary <- intervalData %>%
  summarise(Mean = mean(Interval),
            SEM = sd(Interval) / sqrt(n()))