library(dplyr)


# Load data ---------------------------------------------------------------

# turn direction data is stored in csv files
inputFiles  <- dir('csv/', pattern = "*Turn_Direction.csv")

# load turn direction data
allData <- data.frame()
for (thisFile in 1:length(inputFiles)) {
  tempData <- read.csv(paste0('csv/', inputFiles[thisFile]), header = TRUE)
  
  # add information about subject, session, and electrode to data frame based on
  # file name
  fileNameParts <- strsplit(inputFiles[thisFile], split = "_")
  
  tempData$Subject <- fileNameParts[[1]][1]
  tempData$Session <- fileNameParts[[1]][2]
  allData          <- rbind(allData, tempData)
}
remove(tempData)


# Get trial counts --------------------------------------------------------

allData <- allData %>%
  group_by(Subject, Session, Outcome) 

allOutcomesSummary <- allData %>%
  summarise(Count = n())

turnSummary <- allData %>%
  filter(Outcome == 'WrongDir' | Outcome == 'RightDir') %>%
  ungroup() %>%
  group_by(Subject, Session) %>%
  mutate(Total = n()) %>%
  group_by(Subject, Session, Outcome, Total) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / Total)
