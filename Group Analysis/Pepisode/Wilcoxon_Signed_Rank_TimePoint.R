library(dplyr)
library(ggplot2)
library(reshape2)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/')
filePath <- 'Tall_Pepisode_by_condition_30-Apr-2015.csv'
pepisode <- read.csv(filePath, header = FALSE)
names(pepisode) <- c('Electrode','Condition','FreqBand','TimePoint','Power')

# put our conditions in order
timePointOrder <- c('Pre','Tele','Post')
conditionOrder <- c('NSNT','NSFT','FSNT','FSFT')
freqOrder <- c('Delta','Theta','Alpha','Beta','Gamma')
pepisode$Condition <- factor(pepisode$Condition, levels = conditionOrder)
pepisode$TimePoint <- factor(pepisode$TimePoint, levels = timePointOrder)
pepisode$FreqBand <- factor(pepisode$FreqBand, levels = freqOrder)

# separate data by frequency band
theta <- pepisode %>%
  select(-Condition) %>%
  filter(FreqBand == "Theta") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

delta <- pepisode %>%
  select(-Condition) %>%
  filter(FreqBand == "Delta") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

alpha <- pepisode %>%
  select(-Condition) %>%
  filter(FreqBand == "Alpha") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

beta <- pepisode %>%
  select(-Condition) %>%
  filter(FreqBand == "Beta") %>%
  group_by(Electrode, TimePoint) %>%
  summarise(Power = mean(Power))

gamma <- pepisode %>%
  select(-Condition) %>%
  filter(FreqBand == "Gamma") %>%
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
