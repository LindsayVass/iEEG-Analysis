# Script Name:  2b-Clean-Power-Data.R
# Author:       Lindsay Vass
# Date:         25 September 2015
# Purpose:      This script will clean the data frames produced by "1b-Load-Power-Data.R"

library(dplyr)
load('Rda/allRawPowerData.Rda')


# Functions ---------------------------------------------------------------

factorVar <- function(dataFrame, variable, levelOrder) {
  
  varNum <- which(names(dataFrame) == variable)
  dataFrame[[varNum]] <- factor(dataFrame[[varNum]], levels = levelOrder)
  return(dataFrame)
  
}

orderPowerVar <- function(dataFrame, timePointOrder, trialTypeOrder, spaceOrder, timeOrder) {
  
  dataFrame <- factorVar(dataFrame, 'TimePoint', timePointOrder)
  dataFrame <- factorVar(dataFrame, 'TrialType', trialTypeOrder)
  dataFrame <- factorVar(dataFrame, 'TrialSpaceType', spaceOrder)
  dataFrame <- factorVar(dataFrame, 'TrialTimeType', timeOrder)
  
}

addVariables <- function(dataFrame) {
  
  # Make a new variable for electrode ID (subjectID + Teleporter + Electrode)
  dataFrame <- dataFrame %>%
    mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = "_")) %>%
    mutate(ElectrodeID = factor(ElectrodeID))
  
  # Cut up the frequencies into bands
  frequencies    <- unique(dataFrame$Frequency)
  freqBandBreaks <- c(0, 8, 12, 30, 182)
  freqBandNames  <- c("Delta-Theta","Alpha","Beta","Gamma")
  
  dataFrame <- dataFrame %>%
    mutate(FrequencyBand = cut(Frequency, freqBandBreaks, labels = freqBandNames),
           DepthElectrode = substring(Electrode, 1, 3))
  
  return(dataFrame)
  
}

filterData <- function(dataFrame, subjectID, teleporter, depthElectrode) {
  dataFrame <- dataFrame %>%
    filter(SubjectID == subjectID,
           Teleporter == teleporter,
           DepthElectrode == depthElectrode)
}

makeTrialNumberDf <- function(dataFrame, subjectID, teleporter, depthElectrode) {
  navTrials <- filterData(dataFrame,
                          subjectID,
                          teleporter,
                          depthElectrode)
  badTrials <- data.frame(TrialNumber = seq(1, length(unique(navTrials$RealTrialNumber)), 1))
  navTrials <- cbind(navTrials, badTrials)
}


# Clean data --------------------------------------------------------------

timePointOrder <- c('Pre1', 'Tele', 'Post1')
trialTypeOrder <- c('NSNT','NSFT','FSNT','FSFT')
spaceOrder     <- c('NS','FS')
timeOrder      <- c('NT','FT')

cleanNavPower  <- orderPowerVar(allNavigationPowerData, timePointOrder, trialTypeOrder, spaceOrder, timeOrder) %>%
  select(-TrialSpaceType)
cleanTelePower <- orderPowerVar(allTeleporterPowerData, timePointOrder, trialTypeOrder, spaceOrder, timeOrder)

cleanNavPower  <- addVariables(cleanNavPower)
cleanTelePower <- addVariables(cleanTelePower)

cleanNavTrialNumbers <- realNavigationTrialNumbers %>%
  rename(SubjectID       = Subject,
         Teleporter      = Session,
         DepthElectrode  = Electrode,
         RealTrialNumber = Trial)
cleanTeleTrialNumbers <- realTeleporterTrialNumbers %>%
  rename(SubjectID       = Subject,
         Teleporter      = Session,
         DepthElectrode  = Electrode,
         RealTrialNumber = Trial)

# Fix trial numbers -------------------------------------------------------

navPower  <- data.frame()
telePower <- data.frame()

for (thisRow in 1:nrow(sessionInfo)) {
  
  subjectID  <- sessionInfo$Subject[thisRow]
  teleporter <- sessionInfo$Session[thisRow]
  electrode  <- sessionInfo$Electrode[thisRow]
  
  # get real trial numbers for this session
  navTrials <- makeTrialNumberDf(cleanNavTrialNumbers, 
                          subjectID,
                          teleporter,
                          electrode)
  teleTrials <- makeTrialNumberDf(cleanTeleTrialNumbers, 
                                 subjectID,
                                 teleporter,
                                 electrode)
  
  # get data for this session
  thisnavPower  <- filterData(cleanNavPower,
                                subjectID,
                                teleporter,
                                electrode)
  thistelePower <- filterData(cleanTelePower,
                                 subjectID,
                                 teleporter,
                                 electrode)

  
  # fix the data
  fixnavPower  <- inner_join(navTrials, thisnavPower)
  fixtelePower <- inner_join(teleTrials, thistelePower)
  
  navPower  <- rbind(navPower, fixnavPower)
  telePower <- rbind(telePower, fixtelePower)
}

# Keep trials with both navigation and teleporter data --------------------

navPower <- navPower %>%
  select(ElectrodeID, RealTrialNumber, TimePoint, Frequency, FrequencyBand, Power) %>%
  rename(navPower = Power)
telePower <- telePower %>%
  select(ElectrodeID, RealTrialNumber, TrialSpaceType, TrialTimeType, TrialType, TimePoint, Frequency, FrequencyBand, Power) %>%
  rename(telePower = Power)
allPower <- inner_join(telePower, navPower) %>%
  group_by(ElectrodeID, RealTrialNumber, TimePoint, FrequencyBand)

# Save data ---------------------------------------------------------------

save(file = 'Rda/allCleanPowerData.Rda', list = c('allPower', 'sessionInfo'))




