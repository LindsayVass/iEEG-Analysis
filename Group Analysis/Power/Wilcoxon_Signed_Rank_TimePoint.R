library(dplyr)
library(ggplot2)
library(reshape2)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power/')
filePath <- 'Tall_Power_by_condition_06-May-2015.csv'
powerData <- read.csv(filePath, header = TRUE)

# put our conditions in order
timePointOrder <- c('Pre','Tele','Post')
conditionOrder <- c('NSNT','NSFT','FSNT','FSFT')
freqOrder <- c('deltaBand','thetaBand','alphaBand','betaBand','gammaBand')
powerData$Condition <- factor(powerData$TrialType, levels = conditionOrder)
powerData$TimePoint <- factor(powerData$TimePoint, levels = timePointOrder)
powerData$FreqBand <- factor(powerData$FreqBand, levels = freqOrder)

# separate data by frequency band
theta <- powerData %>%
  select(-TrialType) %>%
  filter(FreqBand == "thetaBand") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

delta <- powerData %>%
  select(-TrialType) %>%
  filter(FreqBand == "deltaBand") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

alpha <- powerData %>%
  select(-TrialType) %>%
  filter(FreqBand == "alphaBand") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

beta <- powerData %>%
  select(-TrialType) %>%
  filter(FreqBand == "betaBand") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

gamma <- powerData %>%
  select(-TrialType) %>%
  filter(FreqBand == "gammaBand") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))


# to wide format
theta_wide <- dcast(theta, Electrode ~ TimePoint, value.var = "Power")
delta_wide <- dcast(delta, Electrode ~ TimePoint, value.var = "Power")
alpha_wide <- dcast(alpha, Electrode ~ TimePoint, value.var = "Power")
beta_wide <- dcast(beta, Electrode ~ TimePoint, value.var = "Power")
gamma_wide <- dcast(gamma, Electrode ~ TimePoint, value.var = "Power")

# Wilcoxon signed-rank test
theta_pre_tele <- wilcox.test(theta_wide$Pre, theta_wide$Tele, paired=TRUE)
theta_tele_post <- wilcox.test(theta_wide$Tele, theta_wide$Post, paired=TRUE)

delta_pre_tele <- wilcox.test(delta_wide$Pre, delta_wide$Tele, paired=TRUE)
delta_tele_post <- wilcox.test(delta_wide$Tele, delta_wide$Post, paired=TRUE)

alpha_pre_tele <- wilcox.test(alpha_wide$Pre, alpha_wide$Tele, paired=TRUE)
alpha_tele_post <- wilcox.test(alpha_wide$Tele, alpha_wide$Post, paired=TRUE)

beta_pre_tele <- wilcox.test(beta_wide$Pre, beta_wide$Tele, paired=TRUE)
beta_tele_post <- wilcox.test(beta_wide$Tele, beta_wide$Post, paired=TRUE)

gamma_pre_tele <- wilcox.test(gamma_wide$Pre, gamma_wide$Tele, paired=TRUE)
gamma_tele_post <- wilcox.test(gamma_wide$Tele, gamma_wide$Post, paired=TRUE)