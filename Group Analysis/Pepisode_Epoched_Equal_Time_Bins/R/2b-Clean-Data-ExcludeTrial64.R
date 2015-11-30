# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         3 June 2015
# Purpose:      This script will clean the data frame produced by "1-Load-Data.R"

library(dplyr)

load('Rda/allRawData_ExcludeTrial64.Rda')


# Functions ---------------------------------------------------------------

makeTrialNumberDf <- function(dataFrame, subjectID, teleporter, depthElectrode) {
  navTrials <- filterData(dataFrame,
                          subjectID,
                          teleporter,
                          depthElectrode)
  badTrials <- data.frame(TrialNumber = seq(1, length(unique(navTrials$RealTrialNumber)), 1))
  navTrials <- cbind(navTrials, badTrials)
}

filterData <- function(dataFrame, subjectID, teleporter, depthElectrode) {
  dataFrame <- dataFrame %>%
   filter(SubjectID == subjectID,
           Teleporter == teleporter,
           DepthElectrode == depthElectrode) %>%
    select(-DepthElectrode)
}


# Clean up trial numbers --------------------------------------------------

cleanTeleTrialNumbers <- realTeleporterTrialNumbers %>%
  rename(SubjectID       = Subject,
         Teleporter      = Session,
         DepthElectrode  = Electrode,
         RealTrialNumber = Trial)
remove(realTeleporterTrialNumbers)

allData <- allData %>%
  mutate(DepthElectrode = substr(Electrode, 1, 3))
  
cleanData <- data.frame()
for (thisRow in 1:nrow(sessionInfo)) {
  
  subjectID  <- sessionInfo$Subject[thisRow]
  teleporter <- sessionInfo$Session[thisRow]
  electrode  <- sessionInfo$Electrode[thisRow]
  
  # get real trial numbers for this session
  teleTrials <- makeTrialNumberDf(cleanTeleTrialNumbers, 
                                  subjectID,
                                  teleporter,
                                  electrode)
  
  # get data for this session
  thisTelePepisode <- filterData(allData,
                                 subjectID,
                                 teleporter,
                                 electrode)
  
  # fix the data
  fixTelePepisode <- inner_join(teleTrials, thisTelePepisode) 
  
  cleanData <- rbind(cleanData, fixTelePepisode)
}


# Put our conditions in order
timePointOrder <- c('Pre1', 'Tele', 'Post1')
cleanData$TimePoint <- factor(cleanData$TimePoint, levels = timePointOrder)

trialTypeOrder <- c('NSNT','NSFT','FSNT','FSFT')
cleanData$TrialType <- factor(cleanData$TrialType, levels = trialTypeOrder)

spaceOrder <- c('NS','FS')
cleanData$TrialSpaceType <- factor(cleanData$TrialSpaceType, levels = spaceOrder)

timeOrder <- c('NT','FT')
cleanData$TrialTimeType <- factor(cleanData$TrialTimeType, levels = timeOrder)

# Make a new variable for electrode ID (subjectID + Teleporter + Electrode)
cleanData <- cleanData %>%
  mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = "_")) %>%
  mutate(ElectrodeID = factor(ElectrodeID))

# Cut up the frequencies into bands
frequencies    <- unique(cleanData$Frequency)
freqBandBreaks <- c(0, 8, 12, 30, 182)
freqBandNames  <- c("Delta-Theta","Alpha","Beta","Gamma")

cleanData <- cleanData %>%
  mutate(FrequencyBand = cut(Frequency, freqBandBreaks, labels = freqBandNames)) %>%
  select(ElectrodeID, RealTrialNumber, TrialSpaceType, TrialTimeType, TimePoint, FrequencyBand, Pepisode) %>%
  rename(TrialNumber = RealTrialNumber)

# Exclude Post for trial 64
noPost <- data.frame(TrialNumber = 64,
                     TimePoint = 'Post1')
cleanData <- anti_join(cleanData, noPost)

save(file = 'Rda/allCleanData_ExcludeTrial64.Rda', list = 'cleanData')