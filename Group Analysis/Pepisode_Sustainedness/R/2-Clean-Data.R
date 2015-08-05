# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         22 July 2015
# Purpose:      This script will clean the data frames produced by "1-Load-Data.R"

library(dplyr)
load('Rda/allRawData.Rda')


# Functions ---------------------------------------------------------------

factorVar <- function(dataFrame, variable, levelOrder) {
  
  varNum <- which(names(dataFrame) == variable)
  dataFrame[[varNum]] <- factor(dataFrame[[varNum]], levels = levelOrder)
  return(dataFrame)
  
}

orderPepisodeVar <- function(dataFrame, timePointOrder, trialTypeOrder, spaceOrder, timeOrder) {
  
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

cleanNavPepisode  <- orderPepisodeVar(allNavigationPepisodeData, timePointOrder, trialTypeOrder, spaceOrder, timeOrder) %>%
  select(-TrialSpaceType)
cleanTelePepisode <- orderPepisodeVar(allTeleporterPepisodeData, timePointOrder, trialTypeOrder, spaceOrder, timeOrder)

cleanNavSustain  <- allNavigationSustainData %>%
  select(-(c(TimePoint, SpaceType))) %>%
  factorVar('TimeType', timeOrder) %>%
  rename(TrialNumber   = Trial, 
         TrialTimeType = TimeType)
cleanTeleSustain <- allTeleporterSustainData %>%
  select(-c(TimePoint)) %>%
  factorVar('SpaceType', spaceOrder) %>%
  factorVar('TimeType', timeOrder) %>%
  rename(TrialNumber    = Trial,
         TrialSpaceType = SpaceType,
         TrialTimeType  = TimeType)

cleanNavPepisode  <- addVariables(cleanNavPepisode)
cleanTelePepisode <- addVariables(cleanTelePepisode)
cleanNavSustain   <- addVariables(cleanNavSustain)
cleanTeleSustain  <- addVariables(cleanTeleSustain)

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

navPepisode  <- data.frame()
telePepisode <- data.frame()
navSustain   <- data.frame()
teleSustain  <- data.frame()

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
  thisNavPepisode  <- filterData(cleanNavPepisode,
                                subjectID,
                                teleporter,
                                electrode)
  thisTelePepisode <- filterData(cleanTelePepisode,
                                 subjectID,
                                 teleporter,
                                 electrode)
  thisNavSustain   <- filterData(cleanNavSustain,
                                 subjectID,
                                 teleporter,
                                 electrode)
  thisTeleSustain  <- filterData(cleanTeleSustain,
                                 subjectID,
                                 teleporter,
                                 electrode)
  
  # fix the data
  fixNavPepisode  <- inner_join(navTrials, thisNavPepisode)
  fixTelePepisode <- inner_join(teleTrials, thisTelePepisode)
  fixNavSustain   <- inner_join(navTrials, thisNavSustain)
  fixTeleSustain  <- inner_join(teleTrials, thisTeleSustain)
  
  navPepisode  <- rbind(navPepisode, fixNavPepisode)
  telePepisode <- rbind(telePepisode, fixTelePepisode)
  navSustain   <- rbind(navSustain, fixNavSustain)
  teleSustain  <- rbind(teleSustain, fixTeleSustain)
}

# Keep trials with both navigation and teleporter data --------------------

navPepisode <- navPepisode %>%
  select(ElectrodeID, RealTrialNumber, TimePoint, Frequency, FrequencyBand, Pepisode) %>%
  rename(NavPepisode = Pepisode)
telePepisode <- telePepisode %>%
  select(ElectrodeID, RealTrialNumber, TrialSpaceType, TrialTimeType, TrialType, TimePoint, Frequency, FrequencyBand, Pepisode) %>%
  rename(TelePepisode = Pepisode)
allPepisode <- inner_join(telePepisode, navPepisode) %>%
  group_by(ElectrodeID, RealTrialNumber, TimePoint, FrequencyBand)

# For the episode duration analysis, I don't think it really matters that the 
# trials are matched, but if so, uncomment the lines below.

# goodTrials <- allPepisode %>%
#   ungroup() %>%
#   select(ElectrodeID, RealTrialNumber) %>%
#   unique()
# 
# navSustain  <- inner_join(navSustain, goodTrials)
# teleSustain <- inner_join(teleSustain, goodTrials)

# Save data ---------------------------------------------------------------

save(file = 'Rda/allCleanData.Rda', list = c('allPepisode', 'navSustain', 'teleSustain', 'sessionInfo'))




