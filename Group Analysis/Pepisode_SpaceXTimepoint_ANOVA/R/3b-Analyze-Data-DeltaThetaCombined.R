# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         16 June 2015
# Purpose:      This script will analyze the data output by "2-Clean-Data.R". It
#               will perform 2-way repeated measure ANOVAs with permutation
#               analyses to identify whether pepisode varies as a function of
#               Space (NS/FS) and Timepoint (Pre/Tele/Post). It will group 
#               together all frequencies within Delta and Theta (< 8 Hz).

library(dplyr)
library(reshape2)
library(permute)
library(car)

load('Rda/allCleanData.Rda')


# Functions ---------------------------------------------------------------

# return data from a specific trial time type
filterRMData <- function(episodeData, timeType) {
  filteredData <- episodeData %>%
    filter(TrialTimeType == timeType) %>%
    ungroup() %>%
    group_by(ElectrodeID, TimeBin)
}

# run two-way anova
runTwoWayAnova <- function(theData, iData) {
  
  # put into wide format for first-level ANOVA
  wideData <- dcast(theData, 
                    ElectrodeID ~ TrialSpaceType + TimeBin, 
                    value.var = 'Pepisode')
  
  # Convert to matrix form and remove ElectrodeID
  wideMatrix <- data.matrix(wideData[2:dim(wideData)[2]])
  
  # Use lm() to generate a linear model on all columns
  firstLM <- lm(wideMatrix ~ 1)
  
  # run second-level ANOVA
  twoWayAnova <- Anova(firstLM, 
                       idata = iData, 
                       idesign = ~spaceType * timePoint)
  twoWayAnova <- summary(twoWayAnova)$univariate.tests
  
}

# run permutation analysis
runPermAnalysis <- function(theData, anovaOutput, control, iData, numPerm = 10000) {
  
  # extract F statistics from anovaOutput so we can compare them to the
  # distribution; skip first row because it's the intercept
  trueF <- data.frame(Contrast = row.names(anovaOutput)[2:nrow(anovaOutput)],
                      FStatistic = anovaOutput[2:nrow(anovaOutput), 5],
                      row.names = NULL)
  
  # initialize array to hold permuted F values
  permF <- data.frame(Contrast = NA, FStatistic = NA)
  
  # create progress bar
  pB <- txtProgressBar(min = 1, max = numPerm, initial = 1, style = 3)
  
  # run two-way anova for each permutation
  for (thisPerm in 1:numPerm) {
    
    # update progress bar
    setTxtProgressBar(pB, thisPerm)
    
    # make a shuffle order for this permutation
    permOrder <- shuffle(nrow(theData), control)
    
    # apply the shuffle to the data
    tempData <- theData
    tempData$Pepisode <- tempData$Pepisode[permOrder]
    
    # run the two-way anova
    # suppress message that first model only has an intercept and warning that
    # first model only has an intercept
    anovaOutput <- suppressWarnings(suppressMessages(runTwoWayAnova(tempData, iData)))
    
    # extract F statistics and add to data holder
    thisPermF <- data.frame(Contrast = row.names(anovaOutput)[2:nrow(anovaOutput)],
                            FStatistic = anovaOutput[2:nrow(anovaOutput), 5],
                            row.names = NULL)
    permF <- rbind(permF, thisPermF)
    
  }
  
  # add the true F data to the perm data
  permF <- rbind(permF, trueF) %>%
    filter(is.na(Contrast) == FALSE) %>%
    group_by(Contrast) %>%
    arrange(desc(FStatistic))
  
  trueP <- data.frame(Contrast = NA, Pcorr = NA)
  for (thisContrast in 1:nrow(trueF)) {
    
    thisF <- permF %>%
      filter(Contrast == trueF$Contrast[thisContrast])
    thisP <- which(thisF$FStatistic == trueF$FStatistic[thisContrast]) / nrow(thisF)
    
    trueP <- rbind(trueP, c(Contrast = as.character(trueF$Contrast[thisContrast]),
                            Pcorr = thisP)) 
    
  }
  
  trueP <- trueP %>%
    filter(is.na(Contrast) == FALSE) 
  trueP$Contrast <- factor(trueP$Contrast)
  
  output <- suppressMessages(inner_join(trueF, trueP))
  
}

# Get mean across values within frequency band ----------------------------

cleanData <- cleanData %>%
  filter(FrequencyBand == "Delta-Theta") %>%
  group_by(ElectrodeID, TrialSpaceType, TrialTimeType, TimeBin) %>%
  summarise(Pepisode = mean(Pepisode))

# Set up data frames for ANOVA ---------------------------------------

# iData for two-way ANOVA
spaceType  <- c(rep("NS", nlevels(cleanData$TimeBin)), 
                rep("FS", nlevels(cleanData$TimeBin)))
timePoint  <- rep(levels(cleanData$TimeBin), 2)
iData <- data.frame(spaceType, timePoint)

# Initialize data frame to save the ANOVA results
anovaResults <- data.frame(TimeType = NA,
                           Contrast = NA,
                           Fstat = NA,
                           Pcorr = NA)

# permutation analysis parameters
numPerm <- 10000

# Run two-way ANOVAs ------------------------------------------------------

for (thisTimeType in 1:nlevels(cleanData$TrialTimeType)) {
  
  cat(paste0('\n',
             levels(cleanData$TrialTimeType)[thisTimeType],
             '\n'))
  
  # extract data from cleanData
  theData <- filterRMData(cleanData, 
                          levels(cleanData$TrialTimeType)[thisTimeType])
  
  # run two-way anova
  anovaOutput <- suppressWarnings(suppressMessages(runTwoWayAnova(theData, iData)))
  
  # set up control for permutation analysis
  control  <- how(within = Within(type = "free"), # freely shuffle timepoints
                  plots = Plots(strata = theData$TrialSpaceType, type = "free"), # freely shuffle space type
                  blocks = theData$ElectrodeID,# don't shuffle between electrodes
                  nperm = numPerm) 
  
  # run permutation analysis 
  permOutput <- runPermAnalysis(theData, anovaOutput, control, iData, numPerm)  
  
  # output data to summary array
  tempResults <- data.frame(TimeType = as.character(levels(cleanData$TrialTimeType)[thisTimeType]),
                            Contrast = as.character(permOutput$Contrast),
                            Fstat = as.numeric(permOutput$FStatistic),
                            Pcorr = as.numeric(permOutput$Pcorr))
  anovaResults <- rbind(anovaResults, tempResults)
  
}


anovaResults <- anovaResults %>% filter(is.na(Contrast) == FALSE)

# save data
save(file = 'Rda/allAnalyzedData_DeltaThetaCombined.Rda', 'anovaResults')
