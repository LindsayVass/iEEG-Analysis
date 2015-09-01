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
#                  that cross a boundary (teleporter entry) and
#                  will test whether the mean duration of these episodes AFTER 
#                  the boundary differs between teleporter epochs and their 
#                  matched navigation epochs.
#                  NB There is no true "boundary" event for the navigation epochs;
#                  we are simply using the same temporal event structure as the
#                  teleporter epochs.



library(permute)
library(reshape2)
library(dplyr)
library(broom)
library(data.table)
load('Rda/allCleanData.Rda')

# Functions ---------------------------------------------------------------

tidyWilcoxon <- function(inputData) {
  NavGtTele <- tidy(inputData, NavGtTele) 
}

replicateData <- function(inputData, nrep) {
  repData <- inputData[rep(1:nrow(inputData), times = nrep), ] %>%
    mutate(Iteration = gl(nrep, nrow(inputData)))
}

shuffleData <- function(thisData, control, variable) {
  ind <- which(names(thisData) == variable)
  shuffleOrder <- shuffle(nrow(thisData), control = control)
  thisData[[ind]] <- thisData[[ind]][shuffleOrder]
  return(thisData)
}

getCorrectedP <- function(permData, trueData, trueOrPerm = "true") {
  
  if (tolower(trueOrPerm) == "perm") {
    origData <- trueData
    trueData <- trueData %>%
      ungroup() %>%
      select(-c(Iteration, ElectrodeIteration))
  }
  
  trueData <- trueData %>%
    mutate(Iteration = 0)
  trueData$Iteration <- factor(trueData$Iteration)
  
  trueDataCorrected <- rbind(permData, trueData)
  trueDataCorrected <- trueDataCorrected %>%
    ungroup() %>%
    group_by(ElectrodeID, FrequencyBand, TimePoint) 
  
  trueDataCorrected <- trueDataCorrected %>%
    arrange(desc(statistic)) %>%
    mutate(CorrP = row_number(desc(statistic)) / n()) %>%
    filter(Iteration == 0)
  
  if (tolower(trueOrPerm) == "perm") {
    trueDataCorrected <- trueDataCorrected %>%
      select(-Iteration) %>%
      inner_join(origData, by = c("ElectrodeID", "FrequencyBand", "statistic", "p.value", "TimePoint"))
  }
  
  return(trueDataCorrected)
}

getElectrodeCorrectedP <- function(trueNSigElectrodes, permNSigElectrodes) {
  permNSigElectrodes <- permNSigElectrodes %>%
    c(trueNSigElectrodes) %>%
    sort(decreasing = TRUE)
  CorrP <- min(which(permNSigElectrodes == trueNSigElectrodes)) / length(permNSigElectrodes)
}





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



# Pepisode Analysis (1) ---------------------------------------------------
nperm <- 1000

meanPepisode <- allPepisode %>%
  summarise(MeanTelePepisode = mean(TelePepisode), MeanNavPepisode = mean(NavPepisode)) %>%
  ungroup() %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint)

# get wilcoxon results for true data
pepisodeTrueData <- meanPepisode %>%
  do(NavGtTele = wilcox.test(.$MeanNavPepisode, .$MeanTelePepisode, alternative = "greater", paired = TRUE)) %>%
  tidy(NavGtTele)

# get wilcoxon results for permuted data
pepisodeData <- meanPepisode %>%
  melt(c('ElectrodeID', 'RealTrialNumber', 'TimePoint', 'FrequencyBand')) %>%
  mutate(Observation = paste(ElectrodeID, FrequencyBand, RealTrialNumber, TimePoint, sep = "_"))

control <- how(within = Within(type = "free"), blocks = pepisodeData$Observation)
permData <- pepisodeData %>%
  ungroup() %>%
  replicateData(nperm)

pepisodePermResults <- vector(mode = "list", length = nperm)
pb <- txtProgressBar (min = 1, max = nperm, style = 3)
for (thisPerm in 1:nperm) {
  setTxtProgressBar(pb, thisPerm)
  pepisodePermResults[[thisPerm]] <- permData %>%
    filter(Iteration == thisPerm) %>%
    shuffleData(control, "value") %>%
    group_by(ElectrodeID, FrequencyBand, RealTrialNumber, TimePoint) %>%
    dcast(ElectrodeID + FrequencyBand + RealTrialNumber + TimePoint ~ variable, value.var = "value") %>%
    group_by(ElectrodeID, FrequencyBand, TimePoint) %>%
    do(NavGtTele = wilcox.test(.$MeanNavPepisode, .$MeanTelePepisode, alternative = "greater", paired = TRUE)) %>%
    tidy(NavGtTele) %>%
    mutate(Iteration = thisPerm)
}
pepisodePermResults <- rbindlist(pepisodePermResults)

# add corrected P values
pepisodeTrueDataCorrected <- pepisodeTrueData %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint) %>%
  getCorrectedP(pepisodePermResults, ., "true")
pepisodeTrueNSigElectrodes <- pepisodeTrueDataCorrected %>%
  ungroup() %>%
  group_by(FrequencyBand, TimePoint) %>%
  filter(CorrP < 0.05) %>%
  summarise(Count = n())

# determine number of significant electrodes expected by chance

pepisodePermAnalysisList <- pepisodePermResults %>%
  select(ElectrodeID, FrequencyBand, TimePoint, Iteration) %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint) %>%
  unique()
pepisodePermMaxElectrodes <- vector(mode = "list", length = nperm)
for (thisPerm in 1:nperm) {
  pepisodePermMaxElectrodes[[thisPerm]] <- pepisodePermAnalysisList %>%
    sample_n(1) %>%
    mutate(ElectrodeIteration = thisPerm)
}
pepisodePermMaxElectrodes <- rbindlist(pepisodePermMaxElectrodes) %>%
  inner_join(pepisodePermResults) %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint) %>%
  mutate(CorrP = row_number(desc(statistic)) / n())
pepisodePermNSigElectrodes <- pepisodePermMaxElectrodes %>%
  group_by(FrequencyBand, TimePoint, ElectrodeIteration) %>%
  filter(CorrP < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(ElectrodeIteration) %>%
  summarise(MaxCount = max(Count))
  
pepisodeTrueNSigElectrodesCorrected <- pepisodeTrueNSigElectrodes %>%
  group_by(FrequencyBand, TimePoint, Count) %>%
  do(CorrP = getElectrodeCorrectedP(.$Count, pepisodePermNSigElectrodes$MaxCount))


# Oscillatory episode post-event duration analysis ------------------------

navNtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Navigation")
navFtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Navigation")

teleNtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Teleportation")
teleFtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Teleportation")

allData <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType, Condition) %>%
  filter(n() > 5) %>%
  summarise(MeanDuration = mean(MeanPostEventDuration),
            SEM = sd(MeanPostEventDuration) / sqrt(n()))
navData <- allData %>%
  filter(Condition == "Navigation")
teleData <- allData %>%
  filter(Condition == "Teleportation")

validData <- inner_join(navData, teleData, by = c('ElectrodeID', 'FrequencyBand', 'TimeType')) %>%
  mutate(Difference = (MeanDuration.y - MeanDuration.x)) %>%
  select(ElectrodeID, TimeType, Difference, FrequencyBand) %>%
  inner_join(allData)

# get wilcoxon results for true data
durationObs <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType, Condition) %>%
  select(-RealTrialNumber) %>%
  dcast(ElectrodeID + FrequencyBand + TimeType ~ Condition, fun.aggregate = length, value.var = "MeanPostEventDuration") %>%
  filter(Navigation >= 5 & Teleportation >=5) %>%
  select(-c(Navigation, Teleportation))
durationTrueData <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  select(-RealTrialNumber) %>%
  inner_join(durationObs) %>%
  dcast(ElectrodeID + FrequencyBand + TimeType ~ Condition, fun.aggregate = list, value.var = "MeanPostEventDuration") %>%
  group_by(ElectrodeID, FrequencyBand, TimeType) %>%
  do(NavGtTele = wilcox.test(unlist(.$Navigation), unlist(.$Teleportation), alternative = "greater")) %>%
  tidy(NavGtTele)

# get wilcoxon results for permuted data
durationData <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  select(-RealTrialNumber) %>%
  inner_join(durationObs) %>%
  mutate(Observation = paste(ElectrodeID, FrequencyBand, TimeType, sep = "_"))

control <- how(within = Within(type = "free"), blocks = durationData$Observation)

permDurationData <- durationData %>%
  ungroup() %>%
  replicateData(nperm)

durationPermResults <- vector(mode = "list", length = nperm)
pb <- txtProgressBar(min = 1, max = nperm, style = 3)
for (thisPerm in 1:nperm) {
  setTxtProgressBar(pb, thisPerm)
  durationPermResults[[thisPerm]] <- permDurationData %>%
    filter(Iteration == thisPerm) %>%
    shuffleData(control, "MeanPostEventDuration") %>%
    dcast(ElectrodeID + FrequencyBand + TimeType ~ Condition, fun.aggregate = list, value.var = "MeanPostEventDuration") %>%
    group_by(ElectrodeID, FrequencyBand, TimeType) %>%
    do(NavGtTele = wilcox.test(unlist(.$Navigation), unlist(.$Teleportation), alternative = "greater")) %>%
    tidy(NavGtTele)
}

durationPermResults <- rbindlist(durationPermResults)

# add corrected p values
durationTrueDataCorrected <- durationTrueData %>%
  group_by(ElectrodeID, FrequencyBand, TimeType) %>%
  getCorrectedP(durationPermResults, ., "true")
durationTrueNSigElectrodes <- durationTrueDataCorrected %>%
  ungroup() %>%
  group_by(FrequencyBand, TimeType) %>%
  filter(CorrP < 0.05) %>%
  summarise(Count = n())

# determine number of significant electrodes expected by chance
durationPermAnalysisList <- durationPermResults %>%
  select(ElectrodeID, FrequencyBand, TimeType, Iteration) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType) %>%
  unique()
durationPermMaxElectrodes <- vector(mode = "list", length = nperm)
for (thisPerm in 1:nperm) {
  durationPermMaxElectrodes[[thisPerm]] <- durationPermAnalysisList %>%
    sample_n(1) %>%
    mutate(ElectrodeIteration = thisPerm)
}
durationPermMaxElectrodes <- rbindlist(durationPermMaxElectrodes) %>%
  inner_join(durationPermResults) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType) %>%
  mutate(CorrP = row_number(desc(statistic)) / n())
durationPermNSigElectrodes <- durationPermMaxElectrodes %>%
  group_by(FrequencyBand, TimeType, ElectrodeIteration) %>%
  filter(CorrP < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(ElectrodeIteration) %>%
  summarise(MaxCount = max(Count))
durationTrueNSigElectrodesCorrected <- durationTrueNSigElectrodes %>%
  group_by(FrequencyBand, TimeType, Count) %>%
  do(CorrP = getElectrodeCorrectedP(.$Count, durationPermNSigElectrodes$MaxCount))



















# Old analysis ------------------------------------------------------------

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
