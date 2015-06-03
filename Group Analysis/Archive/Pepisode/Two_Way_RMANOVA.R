# Script Name:  Two_Wavy_RMANOVA.R
# Author:       Lindsay Vass
# Date:         26 May 2015
# Purpose:      This script will perform two sets of 2-way repeated measures 
#               ANOVAs. The first will test whether pepisode varies as a
#               function of TrialSpaceType x TimePoint. The second will test
#               whether it varies as a function of TrialTimeType x TimePoint.
#               After performing the ANOVAs, it will then run permutation tests
#               to determine the true P values of the effects. Each of these 
#               ANOVAs will be performed separately for each frequency band
#               (Delta, Theta, Alpha, Beta, Gamma).

# Set up workspace and load data ------------------------------------------

library(dplyr)
library(reshape2)
library(permute)
library(car)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/')
load("All_Subjects_Pepisode_CLEAN_3sBuffer2015-05-22.Rda")

# Set up idata structures for ANOVA ---------------------------------------

spaceType  <- c(rep("NS", nlevels(allPepisodeData$TimePoint)), 
               rep("FS", nlevels(allPepisodeData$TimePoint)))
timeType   <- c(rep("NT", nlevels(allPepisodeData$TimePoint)), 
               rep("FT", nlevels(allPepisodeData$TimePoint)))
timePoint  <- rep(levels(allPepisodeData$TimePoint), 2)
spaceIdata <- data.frame(spaceType, timePoint)
timeIdata  <- data.frame(timeType, timePoint)

# Set up parameters for permutation tests ---------------------------------

numPerm <- 10000L # how many permutation iterations to perform

# Perform ANOVAs and permutation tests for each frequency band ------------

# Initialize data frame to save the ANOVA results
anovaResults <- data.frame(FreqBand = NA,
                           ANOVA = NA,
                           Contrast = NA,
                           Fstat = NA,
                           Pcorr = NA)
thisRow <- 1

for (thisFreqBand in 1:nlevels(allPepisodeData$FrequencyBand)){
  
  # Summarize the data for this frequency band
  summaryPepisode <- allPepisodeData %>%
    select(ElectrodeID, 
           TrialSpaceType, 
           TrialTimeType, 
           FrequencyBand, 
           TimePoint, 
           Pepisode) %>%
    filter(FrequencyBand == levels(allPepisodeData$FrequencyBand)[thisFreqBand]) 
  
  # Make separate data frames for the time analysis and space analysis
  spacePepisode <- summaryPepisode %>%
    group_by(ElectrodeID,
             TrialSpaceType,
             TimePoint) %>%
    summarise(Pepisode = mean(Pepisode))
  
  timePepisode <- summaryPepisode %>%
    group_by(ElectrodeID,
             TrialTimeType,
             TimePoint) %>%
    summarise(Pepisode = mean(Pepisode))
  
  # Cast to wide format for ANOVA
  spacePepisodeWide <- dcast(spacePepisode, 
                             ElectrodeID ~ TrialSpaceType + TimePoint, 
                             value.var = "Pepisode")
  timePepisodeWide  <- dcast(timePepisode, 
                             ElectrodeID ~ TrialTimeType + TimePoint, 
                             value.var = "Pepisode")
  
  # Convert to matrix form and remove ElectrodeID
  spaceMatrix <- data.matrix(spacePepisodeWide[2:dim(spacePepisodeWide)[2]])
  timeMatrix  <- data.matrix(timePepisodeWide[2:dim(timePepisodeWide)[2]])
  
  # Use lm() to generate a linear model on all columns
  spaceLM <- lm(spaceMatrix ~ 1)
  timeLM  <- lm(timeMatrix ~ 1)
  
  # Run the ANOVA
  spaceANOVA <- Anova(spaceLM, 
                      idata = spaceIdata, 
                      idesign = ~spaceType * timePoint)
  timeANOVA  <- Anova(timeLM,
                      idata = timeIdata,
                      idesign = ~timeType * timePoint)
  
  # Extract the F statistics
  spaceAnovaResults <- summary(spaceANOVA)$univariate.tests
  timeAnovaResults  <- summary(timeANOVA)$univariate.tests
  
  trueFSpaceSpaceType   <- spaceAnovaResults[2, 5]
  trueFSpaceTimePoint   <- spaceAnovaResults[3, 5]
  trueFSpaceInteraction <- spaceAnovaResults[4, 5]
  
  trueFTimeTimeType     <- timeAnovaResults[2, 5]
  trueFTimeTimePoint    <- timeAnovaResults[3, 5]
  trueFTimeInteraction  <- timeAnovaResults[4, 5]
  
  # Setup control structures for permutation analyses
  controlSpace <- how(within = Within(type = "free"), # freely shuffle timepoints
                      plots = Plots(strata = spacePepisode$TrialSpaceType, type = "free"), # freely shuffle space type
                      blocks = spacePepisode$ElectrodeID,# don't shuffle between electrodes
                      nperm = numPerm) 
  controlTime  <- how(within = Within(type = "free"), # freely shuffle timepoints
                      plots = Plots(strata = timePepisode$TrialTimeType, type = "free"), # freely shuffle space type
                      blocks = timePepisode$ElectrodeID,# don't shuffle between electrodes
                      nperm = numPerm) 
  
  # Initialize vectors to store F statistics from each iteration
  fSpaceSpaceType   <- rep(0, numPerm)
  fSpaceTimePoint   <- rep(0, numPerm)
  fSpaceInteraction <- rep(0, numPerm)
  fTimeTimeType     <- rep(0, numPerm)
  fTimeTimePoint    <- rep(0, numPerm)
  fTimeInteraction  <- rep(0, numPerm)
  
  # Loop through iterations
  for (i in 1:numPerm) {
    
    # Shuffle the pepisode data
    spacePerm <- shuffle(nrow(spacePepisode), control = controlSpace)
    timePerm  <- shuffle(nrow(timePepisode), control = controlTime)
    
    permSpacePepisode <- spacePepisode
    permSpacePepisode$Pepisode <- permSpacePepisode$Pepisode[spacePerm]
    
    permTimePepisode <- timePepisode
    permTimePepisode$Pepisode <- permTimePepisode$Pepisode[timePerm]
    
    # Cast to wide format
    permSpacePepisodeWide <- dcast(permSpacePepisode, 
                                   ElectrodeID ~ TrialSpaceType + TimePoint,
                                   value.var = "Pepisode")
    permTimePepisodeWide  <- dcast(permTimePepisode,
                                   ElectrodeID ~ TrialTimeType + TimePoint,
                                   value.var = "Pepisode")
    
    # Convert to matrix and remove ElectrodeID
    permSpaceMatrix <- data.matrix(permSpacePepisodeWide[2:dim(permSpacePepisodeWide)[2]])
    permTimeMatrix  <- data.matrix(permTimePepisodeWide[2:dim(permTimePepisodeWide)[2]])
    
    # Use lm() to generate a linear model on all columns
    permSpaceLM <- lm(permSpaceMatrix ~ 1)
    permTimeLM  <- lm(permTimeMatrix ~ 1)
    
    # Run the ANOVA
    permSpaceAnova <- Anova(permSpaceLM,
                            idata = spaceIdata,
                            idesign = ~spaceType * timePoint)
    permTimeAnova  <- Anova(permTimeLM,
                            idata = timeIdata,
                            idesign = ~timeType * timePoint)
    
    # Extract the F statistics
    permSpaceStats <- summary(permSpaceAnova)$univariate.tests
    permTimeStats  <- summary(permTimeAnova)$univariate.tests
    
    # Put the F statistics into the distribution
    fSpaceSpaceType[i]   <- permSpaceStats[2, 5]
    fSpaceTimePoint[i]   <- permSpaceStats[3, 5]
    fSpaceInteraction[i] <- permSpaceStats[4, 5]
    
    fTimeTimeType[i]     <- permTimeStats[2, 5]
    fTimeTimePoint[i]    <- permTimeStats[3, 5]
    fTimeInteraction[i]  <- permTimeStats[4, 5]
  } # end loop over iterations
  
  # Add the true F statistics to the distribution
  fSpaceSpaceType[numPerm + 1]   <- trueFSpaceSpaceType
  fSpaceTimePoint[numPerm + 1]   <- trueFSpaceTimePoint
  fSpaceInteraction[numPerm + 1] <- trueFSpaceInteraction
  
  fTimeTimeType[numPerm + 1]     <- trueFTimeTimeType
  fTimeTimePoint[numPerm + 1]    <- trueFTimeTimePoint
  fTimeInteraction[numPerm + 1]  <- trueFTimeInteraction
  
  # Sort the F distribution
  fSpaceSpaceTypeSorted   <- sort(fSpaceSpaceType, decreasing = TRUE)
  fSpaceTimePointSorted   <- sort(fSpaceTimePoint, decreasing = TRUE)
  fSpaceInteractionSorted <- sort(fSpaceInteraction, decreasing = TRUE)
  
  fTimeTimeTypeSorted     <- sort(fTimeTimeType, decreasing = TRUE)
  fTimeTimePointSorted    <- sort(fTimeTimePoint, decreasing = TRUE)
  fTimeInteractionSorted  <- sort(fTimeInteraction, decreasing = TRUE)
  
  # Determine the true P value for our F statistics
  truePSpaceSpaceType <- which(fSpaceSpaceTypeSorted == trueFSpaceSpaceType) / 
    length(fSpaceSpaceTypeSorted)
  truePSpaceTimePoint <- which(fSpaceTimePointSorted == trueFSpaceTimePoint) /
    length(fSpaceTimePointSorted)
  truePSpaceInteraction <- which(fSpaceInteractionSorted == trueFSpaceInteraction) / 
    length(fSpaceInteractionSorted)
  
  truePTimeTimeType <- which(fTimeTimeTypeSorted == trueFTimeTimeType) /
    length(fTimeTimeTypeSorted)
  truePTimeTimePoint <- which(fTimeTimePointSorted == trueFTimeTimePoint) / 
    length(fTimeTimePointSorted)
  truePTimeInteraction <- which(fTimeInteractionSorted == trueFTimeInteraction) /
    length(fTimeInteractionSorted)
  
  # Add the data to the summary data frame
  anovaResults[thisRow, ]     <- c(levels(allPepisodeData$FrequencyBand)[thisFreqBand],
                                   "SpaceXTimepoint",
                                   "SpaceType",
                                   trueFSpaceSpaceType,
                                   truePSpaceSpaceType)
  anovaResults[thisRow + 1, ] <- c(levels(allPepisodeData$FrequencyBand)[thisFreqBand],
                                  "SpaceXTimepoint",
                                  "TimePoint",
                                  trueFSpaceTimePoint,
                                  truePSpaceTimePoint)
  anovaResults[thisRow + 2, ] <- c(levels(allPepisodeData$FrequencyBand)[thisFreqBand],
                                   "SpaceXTimepoint",
                                   "SpaceTypeXTimePoint",
                                   trueFSpaceInteraction,
                                   truePSpaceInteraction)
  anovaResults[thisRow + 3, ] <- c(levels(allPepisodeData$FrequencyBand)[thisFreqBand],
                                   "TimeXTimepoint",
                                   "TimeType",
                                   trueFTimeTimeType,
                                   truePTimeTimeType)
  anovaResults[thisRow + 4, ] <- c(levels(allPepisodeData$FrequencyBand)[thisFreqBand],
                                   "TimeXTimepoint",
                                   "TimePoint",
                                   trueFTimeTimePoint,
                                   truePTimeTimePoint)
  anovaResults[thisRow + 5, ] <- c(levels(allPepisodeData$FrequencyBand)[thisFreqBand],
                                   "TimeXTimepoint",
                                   "TimeTypeXTimePoint",
                                   trueFTimeInteraction,
                                   truePTimeInteraction)
  
  thisRow <- thisRow + 6
                    
  
} # end loop over frequency bands

# Save ANOVA results
today <- Sys.Date()
saveFile <- paste0("All_Subjects_Pepisode_CLEAN_3sBuffer_Two_Way_RMANOVA_", today, ".Rda")
save(anovaResults, file = saveFile)



