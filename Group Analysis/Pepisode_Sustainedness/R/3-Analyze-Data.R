# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         23 July 2015
# Purpose:      This script will analyze the data from 2-Clean-Data.R. It will
#               perform the following analyses:
#               1. Compare pepisode between teleporter and navigation epochs. In
#                  this analysis, we will simply compare pepisode at each time
#                  bin (pre/tele/post) for the teleporter epoch and the preceding
#                  active navigation epoch. 
#               2. Compare the duration of oscillatory episodes for teleporter
#                  and navigation epochs. This analysis will identify episodes
#                  that cross a boundary (either teleporter entry or exit) and
#                  will test whether the mean duration of these episodes differs
#                  between teleporter epochs and their matched navigation epochs.
#                  NB There is no true "boundary" event for the navigation epochs;
#                  we are simply using the same temporal event structure as the
#                  teleporter epochs.
#               3. Similar to #2, we will compare how long the episode lasts
#                  after the boundary event rather than the full duration.


library(permute)
library(reshape2)
library(dplyr)
load('Rda/allCleanData.Rda')

# Functions ---------------------------------------------------------------

testPepisode <- function(dataFrame) {
  
  trueWilcox  <- wilcox.test(dataFrame$MeanTelePepisode, dataFrame$MeanNavPepisode, paired = TRUE)
  trueV       <- trueWilcox$statistic
  truePUncorr <- trueWilcox$p.value
  output      <- data.frame(V = trueV, UncorrP = truePUncorr)
}

runPerms <- function(dataFrame, nperm = 1000) {
  # set up permutations
  tidyData <- melt(dataFrame,
                   id.vars = c("ElectrodeID", "RealTrialNumber", "TimePoint", "FrequencyBand"),
                   measure.vars = c("MeanTelePepisode", "MeanNavPepisode"), 
                   variable.name = "Condition",
                   value.name = "Pepisode")
  
  # run permutations
  permV <- data.frame()
  for (i in 1:nperm) {
    permData <- tidyData
    thisPerm <- shuffle(nrow(tidyData))
    permData$Pepisode <- permData$Pepisode[thisPerm]
    permData <- dcast(permData, ElectrodeID + RealTrialNumber + TimePoint + FrequencyBand ~ Condition,
                      value.var = "Pepisode")
    thisV <- data.frame(V = wilcox.test(permData$MeanTelePepisode, permData$MeanNavPepisode, paired = TRUE, alternative = "less")$statistic, row.names = NULL)
    permV <- rbind(permV, thisV)
  }
  return(permV)
}

sortPerms <- function(permData) {
  rowPos <- data.frame(Row = seq(1, nrow(permData)))
  permData <- permData %>%
    arrange(desc(V))
  vVals <- data.frame(V = unique(permData$V))
  pVals <- data.frame()
  for (thisV in 1:nrow(vVals)) {
    vInd  <- min(which(permData == vVals$V[thisV]))
    pVal  <- data.frame(P = vInd / nrow(permData))
    pVals <- rbind(pVals, pVal)
  }
  vpVals <- cbind(vVals, pVals)
  output <- inner_join(permData, vpVals, by = "V")
}

filterData <- function(dataFrame, electrode, freqband, timepoint) {
  dataFrame <- dataFrame %>%
    filter(ElectrodeID == electrode,
           FrequencyBand == freqband,
           TimePoint == timepoint)
}

meanOscDuration <- function(dataFrame, eventTime, trialTimeType) {
  output <- dataFrame %>%
    filter(Onset < eventTime & Offset > eventTime & TrialTimeType == trialTimeType) %>%
    group_by(ElectrodeID, FrequencyBand, RealTrialNumber) %>%
    summarise(MeanDuration = mean(Duration))
}

testOscDuration <- function(navDf, teleDf, electrode, freqBand, timeType, boundaryType, varName, minTrials = 5) {
  
  navVarInd  <- which(names(navDf) == varName)
  teleVarInd <- which(names(teleDf) == varName)
  navData <- navDf %>%
    filter(ElectrodeID == electrode & FrequencyBand == freqBand)
  teleData <- teleDf %>%
    filter(ElectrodeID == electrode & FrequencyBand == freqBand)
  if (nrow(navData) < minTrials | nrow(teleData) < minTrials){
    output <- data.frame(ElectrodeID = electrode,
                         FrequencyBand = freqBand,
                         TimeType = timeType,
                         BoundaryType = boundaryType,
                         W = NA,
                         UncorrP = NA)
    inputData <- NA
  } else{
    theTest <- wilcox.test(navData[[navVarInd]], teleData[[teleVarInd]])
    output <- data.frame(ElectrodeID = electrode,
                         FrequencyBand = freqBand,
                         TimeType = timeType,
                         BoundaryType = boundaryType,
                         W = theTest$statistic,
                         UncorrP = theTest$p.value)
    
    navDataLab  <- cbind(navData, data.frame(Condition = "Navigation"))
    teleDataLab <- cbind(teleData, data.frame(Condition = "Teleporter"))
    inputData   <- rbind(navDataLab, teleDataLab)
  }
  
  return(list(trueResult = output, inputData = inputData))
  
}

runOscPerms <- function(inputData, trueW, nperm = 1000) {
  output <- data.frame()
  for (i in 1:nperm) {
    permData <- inputData
    thisPerm <- shuffle(nrow(permData))
    permData$MeanDuration <- permData$MeanDuration[thisPerm]
    
    navData <- permData %>%
      filter(Condition == "Navigation")
    teleData <- permData %>%
      filter(Condition == "Teleporter")
    
    thisW <- data.frame(V = wilcox.test(teleData$MeanDuration, navData$MeanDuration, alternative = "less")$statistic, row.names = NULL)
    output <- rbind(output, thisW)
  }
  sortedPerms <- sortPerms(output) %>%
    rename(W = V)
  sortedPermsAndTrue <- rbind(output, trueW) %>%
    arrange(desc(V))
  trueP <- data.frame(CorrP = min(which(sortedPermsAndTrue == trueW)) / nrow(sortedPermsAndTrue))
  return(list(sortedPerms = sortedPerms, trueP = trueP))
}

meanPostEventOscDuration <- function(dataFrame, eventTime, trialTimeType) {
  output <- dataFrame %>%
    filter(Onset < eventTime & Offset > eventTime & TrialTimeType == trialTimeType) %>%
    mutate(PostEventDuration = Offset - eventTime) %>%
    group_by(ElectrodeID, FrequencyBand, RealTrialNumber) %>%
    summarise(MeanPostEventDuration = mean(PostEventDuration))
}

calcMaxSigElectrodes <- function(dataFrame, nPerm = 1000) {
  maxElec <- data.frame()
  for (i in 1:nPerm) {
    theSample <- dataFrame %>%
      sample_n(1) 
    thisMax <-  data.frame(V = max(theSample$NSigElectrodes))
    maxElec <- rbind(maxElec, thisMax)
  }
  return(maxElec)
}

# Pepisode Analysis -------------------------------------------------------

meanPepisode <- allPepisode %>%
  summarise(MeanTelePepisode = mean(TelePepisode), MeanNavPepisode = mean(NavPepisode)) %>%
  ungroup() %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint)

pepisodeResults <- data.frame()
sigElectrodeResults <- data.frame()
minTrials <- 5

for (thisFreqBand in 1:nlevels(meanPepisode$FrequencyBand)) {
  
  cat(paste('\nFrequency Band', thisFreqBand, 'of', nlevels(meanPepisode$FrequencyBand)))
  
  for (thisTimePoint in 1:nlevels(meanPepisode$TimePoint)) {
    
    cat(paste('\nTime Point', thisTimePoint, 'of', nlevels(meanPepisode$TimePoint)))
    
    allPermData <- data.frame()
    
    for (thisElectrode in 1:nlevels(meanPepisode$ElectrodeID)) {
      
      cat(paste('\nElectrode', thisElectrode, 'of', nlevels(meanPepisode$ElectrodeID)))
      
      # extract data
      electrode <- levels(meanPepisode$ElectrodeID)[thisElectrode]
      freqband  <- levels(meanPepisode$FrequencyBand)[thisFreqBand]
      timepoint <- levels(meanPepisode$TimePoint)[thisTimePoint]
      thisData  <- filterData(meanPepisode, electrode, freqband, timepoint)
      
      if (nrow(thisData) < minTrials) {
        
        trueResult <- data.frame(V = NA,
                                 UncorrP = NA,
                                 ElectrodeID = electrode,
                                 FrequencyBand = freqband,
                                 TimePoint = timepoint,
                                 CorrP = NA)
        next
      }
      
      # run wilcoxon
      trueResult <- testPepisode(thisData) %>%
        mutate(ElectrodeID   = electrode,
               FrequencyBand = freqband,
               TimePoint     = timepoint)
      
      # run permutations
      permData <- runPerms(thisData)
      permDataPlusTrue   <- rbind(permData, trueResult$V)
      sortedPermPlusTrue <- sortPerms(permDataPlusTrue)
      trueP <- min(which(sortedPermPlusTrue == trueResult$V)) / nrow(sortedPermPlusTrue)
      trueResult <- trueResult %>%
        mutate(CorrP = trueP)
      pepisodeResults <- rbind(pepisodeResults, trueResult)
      
      # save perm data for later perms on # of significant electrodes      
      sortedPermData <- sortPerms(permData) %>%
        mutate(ElectrodeID = electrode)
      allPermData    <- rbind(allPermData, sortedPermData)
      
    }
    
    # randomly select one permutation from each electrode to build a distribution
    # of how many significant electrodes you get by chance
    nperm <- 1000
    allPermData <- allPermData %>%
      group_by(ElectrodeID)
    elecPerm <- data.frame()
    for (thisPerm in 1:nperm) {
      sampleData <- allPermData %>%
        sample_n(1) %>%
        filter(P < 0.05)
      thisPerm <- data.frame(V = nrow(sampleData))
      elecPerm <- rbind(elecPerm, thisPerm)
    }
    sortedElecPerm <- sortPerms(elecPerm) %>%
      mutate(FrequencyBand = freqband,
             TimePoint = timepoint)
    sigElectrodeResults <- rbind(sigElectrodeResults, sortedElecPerm)
  }
}

# determine whether more electrodes than expected by chance
sigElectrodeResults <- sigElectrodeResults %>%
  group_by(FrequencyBand, TimePoint) %>%
  rename(NSigElectrodes = V)

maxSigElectrodePerms <- calcMaxSigElectrodes(sigElectrodeResults) %>%
  sortPerms() %>%
  rename(NSigElectrodes = V) %>%
  unique()

pepisodeSigElectrodes <- pepisodeResults %>%
  filter(CorrP < 0.05) %>%
  group_by(FrequencyBand, TimePoint) %>%
  summarise(NSigElectrodes = n())
pepisodeSigElectrodes <- inner_join(pepisodeSigElectrodes, maxSigElectrodePerms)

# Oscillatory episode duration analysis -----------------------------------

navNtEntryDur <- meanOscDuration(navSustain, 0, "NT")
navFtEntryDur <- meanOscDuration(navSustain, 0, "FT")
navNtExitDur  <- meanOscDuration(navSustain, 1830, "NT")
navFtExitDur  <- meanOscDuration(navSustain, 2830, "FT")

teleNtEntryDur <- meanOscDuration(teleSustain, 0, "NT")
teleFtEntryDur <- meanOscDuration(teleSustain, 0, "FT")
teleNtExitDur  <- meanOscDuration(teleSustain, 1830, "NT")
teleFtExitDur  <- meanOscDuration(teleSustain, 2830, "FT")

analysisInfo <- list(navInput     = list(navNtEntryDur, navNtExitDur, navFtEntryDur, navFtExitDur),
                     teleInput    = list(teleNtEntryDur, teleNtExitDur, teleFtEntryDur, teleFtExitDur),
                     timeType     = list("NT", "NT", "FT", "FT"),
                     boundaryType = list("Entry", "Exit", "Entry", "Exit"))

meanEpisodeDuration    <- data.frame()
allEpisodePermData     <- data.frame()
sigElectrodeEpisodeDur <- data.frame()
for (thisFreqBand in 1:nlevels(navNtEntryDur$FrequencyBand)) {
  cat('\n\nFrequency Band', thisFreqBand, 'of', nlevels(navNtEntryDur$FrequencyBand))
  for (thisAnalysis in 1:length(analysisInfo)) {
    cat('\nAnalysis', thisAnalysis, 'of', length(analysisInfo), '\nElectrode ')
    for (thisElectrode in 1:nlevels(navNtEntryDur$ElectrodeID)) {
      cat(thisElectrode, ' ')
      # run wilcoxon
      electrode <- levels(navNtEntryDur$ElectrodeID)[thisElectrode]
      freqBand  <- levels(navNtEntryDur$FrequencyBand)[thisFreqBand]
      
      thisResult <- testOscDuration(analysisInfo$navInput[[thisAnalysis]],
                                    analysisInfo$teleInput[[thisAnalysis]],
                                    electrode,
                                    freqBand,
                                    analysisInfo$timeType[[thisAnalysis]],
                                    analysisInfo$boundaryType[[thisAnalysis]],
                                    "MeanDuration")
      
      # run permutations
      if (is.na(thisResult$inputData) == FALSE) {
        permData <- thisResult$inputData
        thePerms <- runOscPerms(permData, thisResult$trueResult$W)
        thePerms$sortedPerms <- cbind(thePerms$sortedPerms, data.frame(ElectrodeID = electrode))
        allEpisodePermData   <- rbind(allEpisodePermData, thePerms$sortedPerms)
        
        thisResult$trueResult <- cbind(thisResult$trueResult, thePerms$trueP)
        meanEpisodeDuration   <- rbind(meanEpisodeDuration, thisResult$trueResult)
        
      } else {
        thisResult$trueResult <- cbind(thisResult$trueResult, data.frame(CorrP = NA))
        meanEpisodeDuration   <- rbind(meanEpisodeDuration, thisResult$trueResult)
      }
    }
    
    # randomly select one permutation from each electrode to build a distribution
    # of how many significant electrodes you get by chance
    allEpisodePermData <- allEpisodePermData %>%
      group_by(ElectrodeID)
    elecPerm <- data.frame()
    for (thisPerm in 1:nperm) {
      sampleData <- allEpisodePermData %>%
        sample_n(1) %>%
        filter(P < 0.05)
      thisPerm <- data.frame(V = nrow(sampleData))
      elecPerm <- rbind(elecPerm, thisPerm)
    }
    sortedElecPerm <- sortPerms(elecPerm) %>%
      mutate(FrequencyBand = freqBand,
             TimeType = analysisInfo$timeType[[thisAnalysis]],
             BoundaryType = analysisInfo$boundaryType[[thisAnalysis]])
    sigElectrodeEpisodeDur <- rbind(sigElectrodeEpisodeDur, sortedElecPerm)
  }
}

# determine whether more electrodes than expected by chance
sigElectrodeEpisodeDur <- sigElectrodeEpisodeDur %>%
  group_by(FrequencyBand, TimeType, BoundaryType) %>%
  rename(NSigElectrodes = V)

maxSigElectrodePerms <- calcMaxSigElectrodes(sigElectrodeEpisodeDur) %>%
  sortPerms() %>%
  rename(NSigElectrodes = V) %>%
  unique()

episodeDurSigElectrodes <- meanEpisodeDuration %>%
  filter(CorrP < 0.05) %>%
  group_by(FrequencyBand, TimeType, BoundaryType) %>%
  summarise(NSigElectrodes = n())
episodeDurSigElectrodes <- inner_join(episodeDurSigElectrodes, maxSigElectrodePerms)


# Oscillatory episode post-event duration analysis ------------------------

navNtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "NT")
navFtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "FT")
navNtPostExitDur  <- meanPostEventOscDuration(navSustain, 1830, "NT")
navFtPostExitDur  <- meanPostEventOscDuration(navSustain, 2830, "FT")

teleNtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "NT")
teleFtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "FT")
teleNtPostExitDur  <- meanPostEventOscDuration(teleSustain, 1830, "NT")
teleFtPostExitDur  <- meanPostEventOscDuration(teleSustain, 2830, "FT")

analysisInfo <- list(navInput     = list(navNtPostEntryDur, navNtPostExitDur, navFtPostEntryDur, navFtPostExitDur),
                     teleInput    = list(teleNtPostEntryDur, teleNtPostExitDur, teleFtPostEntryDur, teleFtPostExitDur),
                     timeType     = list("NT", "NT", "FT", "FT"),
                     boundaryType = list("Entry", "Exit", "Entry", "Exit"))

meanPostEpisodeDuration    <- data.frame()
allPostEpisodePermData     <- data.frame()
sigElectrodePostEpisodeDur <- data.frame()

for (thisFreqBand in 1:nlevels(navNtPostEntryDur$FrequencyBand)) {
  cat('\n\nFrequency Band', thisFreqBand, 'of', nlevels(navNtPostEntryDur$FrequencyBand))
  for (thisAnalysis in 1:length(analysisInfo)) {
    cat('\nAnalysis', thisAnalysis, 'of', length(analysisInfo), '\nElectrode ')
    for (thisElectrode in 1:nlevels(navNtPostEntryDur$ElectrodeID)) {
      cat(thisElectrode, ' ')
      # run wilcoxon
      electrode <- levels(navNtPostEntryDur$ElectrodeID)[thisElectrode]
      freqBand  <- levels(navNtPostEntryDur$FrequencyBand)[thisFreqBand]
      
      thisResult <- testOscDuration(analysisInfo$navInput[[thisAnalysis]],
                                    analysisInfo$teleInput[[thisAnalysis]],
                                    electrode,
                                    freqBand,
                                    analysisInfo$timeType[[thisAnalysis]],
                                    analysisInfo$boundaryType[[thisAnalysis]],
                                    "MeanPostEventDuration")
      
      # run permutations
      if (is.na(thisResult$inputData) == FALSE) {
        permData <- thisResult$inputData %>%
          rename(MeanDuration = MeanPostEventDuration)
        thePerms <- runOscPerms(permData, thisResult$trueResult$W)
        thePerms$sortedPerms   <- cbind(thePerms$sortedPerms, data.frame(ElectrodeID = electrode))
        allPostEpisodePermData <- rbind(allPostEpisodePermData, thePerms$sortedPerms)
        
        thisResult$trueResult   <- cbind(thisResult$trueResult, thePerms$trueP)
        meanPostEpisodeDuration <- rbind(meanPostEpisodeDuration, thisResult$trueResult)
        
      } else {
        thisResult$trueResult <- cbind(thisResult$trueResult, data.frame(CorrP = NA))
        meanPostEpisodeDuration   <- rbind(meanPostEpisodeDuration, thisResult$trueResult)
      }
    }
    
    # randomly select one permutation from each electrode to build a distribution
    # of how many significant electrodes you get by chance
    allPostEpisodePermData <- allPostEpisodePermData %>%
      group_by(ElectrodeID)
    elecPerm <- data.frame()
    for (thisPerm in 1:nperm) {
      sampleData <- allPostEpisodePermData %>%
        sample_n(1) %>%
        inner_join(allPostEpisodePermData, by = c("W", "P", "ElectrodeID")) %>% # get corrected P value
        unique()  %>%
        filter(P < 0.05)
      thisPerm <- data.frame(V = nrow(sampleData))
      elecPerm <- rbind(elecPerm, thisPerm)
    }
    sortedElecPerm <- sortPerms(elecPerm) %>%
      mutate(FrequencyBand = freqBand,
             TimeType = analysisInfo$timeType[[thisAnalysis]],
             BoundaryType = analysisInfo$boundaryType[[thisAnalysis]])
    sigElectrodePostEpisodeDur <- rbind(sigElectrodePostEpisodeDur, sortedElecPerm)
  }
}

# determine whether more electrodes than expected by chance
postEpisodeDurSigElectrodes <- meanPostEpisodeDuration %>%
  filter(CorrP < 0.05) %>%
  group_by(FrequencyBand, TimeType, BoundaryType) %>%
  summarise(NSigElectrodes = n()) %>%
  group_by(FrequencyBand, TimeType, BoundaryType)

sigElectrodePostEpisodeDur <- sigElectrodePostEpisodeDur %>%
  group_by(FrequencyBand, TimeType, BoundaryType) %>%
  rename(NSigElectrodes = V)

permList <- vector(mode = "list", length = nperm)
for (i in 1:nperm) {
  theSample <- sigElectrodePostEpisodeDur %>%
    sample_n(1)
  permList[[i]] <- max(theSample$NSigElectrodes)
}
permVec <- unlist(permList) %>%
  sort(decreasing = TRUE)
maxElectrodes <- data.frame()
for (i in min(permVec):max(permVec)) {
  thisResult <- data.frame(NSigElectrodes = i,
                           ElecP = min(which(permVec == i)) / length(permVec))
  maxElectrodes <- rbind(maxElectrodes, thisResult)
}
postEpisodeDurSigElectrodes <- inner_join(postEpisodeDurSigElectrodes, maxElectrodes)


# Save data ---------------------------------------------------------------

save(file = 'Rda/allAnalyzedData.Rda', list = c('pepisodeSigElectrodes', 'pepisodeResults', 'episodeDurSigElectrodes', 'meanEpisodeDuration', 'postEpisodeDurSigElectrodes', 'meanPostEpisodeDuration', 'maxSigElectrodePerms'))
