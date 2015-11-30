library(dplyr)
library(ggplot2)

load('Rda/UCDMC14_TeleporterB_LAD1.Rda')

powerDf <- as.data.frame(powerDf) %>%
  group_by(Frequency)

sampRate <- 512
minuteSamps <- sampRate * 60
minutes <- c(1, 2.5, 5, 10, 20)
sampsList <- minutes * minuteSamps

numIterations <- 100
allSubsampPowerData <- vector(mode = 'list', length = length(minutes))
for (thisLength in 1:length(minutes)) {
  subsampPowerData <- vector(mode = 'list', length = numIterations)
  for (i in 1:numIterations) {
    subsampPowerData[[i]] <- sample_n(powerDf, sampsList[thisLength]) %>%
      summarise(Mean = mean(Power),
                SEM = sd(Power) / sqrt(n())) %>%
      mutate(Iteration = i)
  }
  subsampPowerData <- data.table::rbindlist(subsampPowerData) %>%
    mutate(Minutes = minutes[thisLength])
  allSubsampPowerData[[thisLength]] <- subsampPowerData
}
allPowerData <- data.table::rbindlist(allSubsampPowerData)

summaryPowerData <- allPowerData %>%
  ungroup() %>%
  group_by(Minutes, Frequency) %>%
  summarise(MeanPow = mean(Mean),
            SEMPow = sd(Mean) / sqrt(n()))

## plot data
p <- ggplot(summaryPowerData, aes(x = Frequency, y = MeanPow, ymin = MeanPow - SEMPow, ymax = MeanPow + SEMPow)) +
  geom_line() +
  geom_linerange() +
  scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128)) +
  scale_y_log10() + 
  facet_wrap(~ Minutes)