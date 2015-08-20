# Script Name:  2-Run-Regressions.R
# Author:       Lindsay Vass
# Date:         7 August 2015
# Purpose:      This script will build and run the regressions using the data 
#               loaded by 1-Load-Data.R.

library(dplyr)
library(broom)
library(permute)

#load('Rda/allRawData.Rda')

# Functions ---------------------------------------------------------------

replicateData <- function(inputData, nrep) {
  repData <- inputData[rep(1:nrow(inputData), times = nrep), ] %>%
    mutate(Iteration = gl(nrep, nrow(inputData)))
}

shuffleData <- function(thisData, control, nrep) {
  permData <- replicateData(thisData, nrep)
  shuffleOrder <- shuffleSet(nrow(thisData), nrep, control = control) %>%
    t() %>%
    as.vector()
  permData$Power <- permData$Power[shuffleOrder]
  return(permData)
}

runRegression <- function(inputData) {
  stats <- inputData %>%
    do(Fit = lm(Power ~ Speed + Landmarks, data = .)) %>%
    tidy(Fit) %>%
    filter(term != "(Intercept)")
  return(stats)
}

getCorrectedP <- function(permData, trueData, trueOrPerm = "true") {
  
  if (tolower(trueOrPerm) == "perm") {
    origData <- trueData
    trueData <- trueData %>%
    select(-c(Iteration, ElectrodeIteration))
  }
  
  trueData <- trueData %>%
    mutate(Iteration = 0)
  trueData$Iteration <- factor(trueData$Iteration)
  
  trueDataCorrected <- rbind(permData, trueData)
  trueDataCorrected$term <- factor(trueDataCorrected$term)
  trueDataCorrected <- trueDataCorrected %>%
    ungroup() %>%
    group_by(ElectrodeID, Frequency, term) %>%
    arrange(desc(statistic)) %>%
    mutate(CorrP = row_number(desc(statistic)) / n()) %>%
    filter(Iteration == 0)
  if (tolower(trueOrPerm) == "perm") {
    trueDataCorrected <- trueDataCorrected %>%
      select(-Iteration) %>%
      inner_join(origData)
  }
  return(trueDataCorrected)
}

getElectrodeCorrectedP <- function(trueNSigElectrodes, permNSigElectrodes) {
  permNSigElectrodes <- permNSigElectrodes %>%
    c(trueNSigElectrodes) %>%
    sort(decreasing = TRUE)
  CorrP <- min(which(permNSigElectrodes == trueNSigElectrodes)) / length(permNSigElectrodes)
}

termLabeller <- function(variable, value) {
  return(termNames[value])
}

# Run regressions on true data --------------------------------------------

allRawData <- allRawData %>%
  mutate(ElectrodeID = paste(Subject, Session, Electrode, sep = '_'))

noCentralData <- allRawData %>%
  filter(Landmarks != "CENTRAL") %>%
  group_by(ElectrodeID, Frequency)
noCentralStats <- runRegression(noCentralData)
noCentralStats$ElectrodeID <- factor(noCentralStats$ElectrodeID)
noCentralStats$Frequency   <- factor(noCentralStats$Frequency)


# Run regressions on permuted data ----------------------------------------

# I think for computational reasons (i.e., not enough memory), permutations need
# to be done separately in a loop
nperm <- 1000
control <- how(within = Within(type = "free"))
permList <- vector("list", nlevels(noCentralStats$ElectrodeID) * nlevels(noCentralStats$Frequency))
i <- 1
for (thisElectrode in 1:nlevels(noCentralStats$ElectrodeID)) {
  for (thisFrequency in 1:nlevels(noCentralStats$Frequency)) {

    cat('\n\nAnalysis', i, 'of', length(permList))
    thisData <- allRawData %>%
      filter(ElectrodeID == levels(noCentralStats$ElectrodeID)[thisElectrode],
             Frequency   == levels(noCentralStats$Frequency)[thisFrequency],
             Landmarks != "CENTRAL") %>%
      shuffleData(control, nperm) %>%
      group_by(ElectrodeID, Frequency, Iteration)
#       do(NoCentralFit = lm(Power ~ Speed + Landmarks, data = .))
#     thisStats <- tidy(thisData, NoCentralFit) %>%
#       filter(term != "(Intercept)")
#     permList[[i]] <- thisStats
    permList[[i]] <- runRegression(thisData)
    i <- i + 1
    
  }
}

save(file = 'Rda/allRegressionData.Rda', list = c('noCentralStats', 'permList'))

# Get corrected P value for true regressions ------------------------------

permDf <- do.call(rbind.data.frame, permList) %>%
  group_by(ElectrodeID, Frequency, Iteration, term)

noCentralStatsCorrected <- getCorrectedP(permDf, noCentralStats)

# Get # of significant electrodes -----------------------------------------

noCentralNSigElectrodes <- noCentralStatsCorrected %>%
  ungroup() %>%
  group_by(Frequency, term) %>%
  filter(CorrP < 0.05) %>%
  summarise(Count = n())
noCentralNSigElectrodes$Frequency <- as.numeric(noCentralNSigElectrodes$Frequency)
noCentralNSigElectrodes <- noCentralNSigElectrodes %>%
  arrange(Frequency)


# Build distribution of # of significant electrodes -----------------------

permAnalysisList <- permDf %>%
  select(ElectrodeID, Frequency, Iteration) %>%
  group_by(ElectrodeID, Frequency) %>%
  unique()

noCentralMaxElectrodes <- vector(mode = "list", length = nperm)
for (thisPerm in 1:nperm) {
  noCentralMaxElectrodes[[thisPerm]] <- permAnalysisList %>%
    sample_n(1) %>%
    mutate(ElectrodeIteration = thisPerm)
}
noCentralMaxElectrodes <- do.call(rbind.data.frame, noCentralMaxElectrodes) %>%
  inner_join(permDf) %>%
  group_by(ElectrodeID, Frequency, term)

noCentralMaxElectrodesCorrected <- noCentralMaxElectrodes %>%
  group_by(ElectrodeIteration) %>%
  do(getCorrectedP(permDf, ., "perm"))

noCentralPermNSigElectrodes <- noCentralMaxElectrodesCorrected %>%
  group_by(Frequency, term, ElectrodeIteration) %>%
  filter(CorrP < 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(ElectrodeIteration) %>%
  summarise(MaxCount = max(Count))


# Get corrected P values for number of significant electrodes -------------
noCentralNSigElectrodesCorrected <- noCentralNSigElectrodes %>%
  group_by(Frequency, term, Count) %>%
  do(CorrP = getElectrodeCorrectedP(.$Count, noCentralPermNSigElectrodes$MaxCount))


# Plot histograms ---------------------------------------------------------

termNames <- c("Landmarks (Rich > Poor)", "Speed")

p <- ggplot(noCentralNSigElectrodesCorrected, aes(x = Frequency, y = Count)) +
  geom_bar(stat = "identity") +
  facet_grid(term ~ ., labeller = termLabeller) +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = "red") +
  ylab("# of Significant Electrodes") +
  theme(text = element_text(size = 24))
p
ggsave(filename = 'Figures/Histogram_NoCentralPlaza.png')
save(file = 'Rda/allRegressionData.Rda', list = c('noCentralMaxElectrodesCorrected', 'noCentralNSigElectrodesCorrected', 'noCentralPermNSigElectrodes', 'noCentralStatsCorrected'))