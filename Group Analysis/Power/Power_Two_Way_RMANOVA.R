# Script Name:  Power_Two_Way_RMANOVA.R
# Author:       Lindsay Vass
# Date:         27 May 2015
# Purpose:      This script will perform two sets of 2-way repeated measures 
#               ANOVAs. The first will test whether power varies as a
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

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power/')
load("All_Subjects_Power_CLEAN_3sBuffer2015-05-27.Rda")

# Set up idata structures for ANOVA ---------------------------------------

spaceType  <- c(rep("NS", nlevels(allPowerData$TimePoint)), 
               rep("FS", nlevels(allPowerData$TimePoint)))
timeType   <- c(rep("NT", nlevels(allPowerData$TimePoint)), 
               rep("FT", nlevels(allPowerData$TimePoint)))
timePoint  <- rep(levels(allPowerData$TimePoint), 2)
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

for (thisFreqBand in 1:nlevels(allPowerData$FrequencyBand)){
  
  # Summarize the data for this frequency band
  summaryPower <- allPowerData %>%
    select(ElectrodeID, 
           TrialSpaceType, 
           TrialTimeType, 
           FrequencyBand, 
           TimePoint, 
           Power) %>%
    filter(FrequencyBand == levels(allPowerData$FrequencyBand)[thisFreqBand]) 
  
  # Make separate data frames for the time analysis and space analysis
  spacePower <- summaryPower %>%
    group_by(ElectrodeID,
             TrialSpaceType,
             TimePoint) %>%
    summarise(Power = mean(Power))
  
  timePower <- summaryPower %>%
    group_by(ElectrodeID,
             TrialTimeType,
             TimePoint) %>%
    summarise(Power = mean(Power))
  
  # Cast to wide format for ANOVA
  spacePowerWide <- dcast(spacePower, 
                             ElectrodeID ~ TrialSpaceType + TimePoint, 
                             value.var = "Power")
  timePowerWide  <- dcast(timePower, 
                             ElectrodeID ~ TrialTimeType + TimePoint, 
                             value.var = "Power")
  
  # Convert to matrix form and remove ElectrodeID
  spaceMatrix <- data.matrix(spacePowerWide[2:dim(spacePowerWide)[2]])
  timeMatrix  <- data.matrix(timePowerWide[2:dim(timePowerWide)[2]])
  
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
                      plots = Plots(strata = spacePower$TrialSpaceType, type = "free"), # freely shuffle space type
                      blocks = spacePower$ElectrodeID,# don't shuffle between electrodes
                      nperm = numPerm) 
  controlTime  <- how(within = Within(type = "free"), # freely shuffle timepoints
                      plots = Plots(strata = timePower$TrialTimeType, type = "free"), # freely shuffle space type
                      blocks = timePower$ElectrodeID,# don't shuffle between electrodes
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
    
    # Shuffle the Power data
    spacePerm <- shuffle(nrow(spacePower), control = controlSpace)
    timePerm  <- shuffle(nrow(timePower), control = controlTime)
    
    permSpacePower <- spacePower
    permSpacePower$Power <- permSpacePower$Power[spacePerm]
    
    permTimePower <- timePower
    permTimePower$Power <- permTimePower$Power[timePerm]
    
    # Cast to wide format
    permSpacePowerWide <- dcast(permSpacePower, 
                                   ElectrodeID ~ TrialSpaceType + TimePoint,
                                   value.var = "Power")
    permTimePowerWide  <- dcast(permTimePower,
                                   ElectrodeID ~ TrialTimeType + TimePoint,
                                   value.var = "Power")
    
    # Convert to matrix and remove ElectrodeID
    permSpaceMatrix <- data.matrix(permSpacePowerWide[2:dim(permSpacePowerWide)[2]])
    permTimeMatrix  <- data.matrix(permTimePowerWide[2:dim(permTimePowerWide)[2]])
    
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
  anovaResults[thisRow, ]     <- c(levels(allPowerData$FrequencyBand)[thisFreqBand],
                                   "SpaceXTimepoint",
                                   "SpaceType",
                                   trueFSpaceSpaceType,
                                   truePSpaceSpaceType)
  anovaResults[thisRow + 1, ] <- c(levels(allPowerData$FrequencyBand)[thisFreqBand],
                                  "SpaceXTimepoint",
                                  "TimePoint",
                                  trueFSpaceTimePoint,
                                  truePSpaceTimePoint)
  anovaResults[thisRow + 2, ] <- c(levels(allPowerData$FrequencyBand)[thisFreqBand],
                                   "SpaceXTimepoint",
                                   "SpaceTypeXTimePoint",
                                   trueFSpaceInteraction,
                                   truePSpaceInteraction)
  anovaResults[thisRow + 3, ] <- c(levels(allPowerData$FrequencyBand)[thisFreqBand],
                                   "TimeXTimepoint",
                                   "TimeType",
                                   trueFTimeTimeType,
                                   truePTimeTimeType)
  anovaResults[thisRow + 4, ] <- c(levels(allPowerData$FrequencyBand)[thisFreqBand],
                                   "TimeXTimepoint",
                                   "TimePoint",
                                   trueFTimeTimePoint,
                                   truePTimeTimePoint)
  anovaResults[thisRow + 5, ] <- c(levels(allPowerData$FrequencyBand)[thisFreqBand],
                                   "TimeXTimepoint",
                                   "TimeTypeXTimePoint",
                                   trueFTimeInteraction,
                                   truePTimeInteraction)
  
  thisRow <- thisRow + 6
                    
  
} # end loop over frequency bands

# Save ANOVA results
today <- Sys.Date()
saveFile <- paste0("All_Subjects_Power_CLEAN_3sBuffer_Two_Way_RMANOVA_", today, ".Rda")
save(anovaResults, file = saveFile)



