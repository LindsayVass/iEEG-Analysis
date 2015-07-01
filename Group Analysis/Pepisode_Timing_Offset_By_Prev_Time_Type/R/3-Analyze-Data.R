# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         30 June 2015
# Purpose:      This script will analyze the data output from 2-Clean-Data.R. It
#               will run a chi-square test to determine whether the distribution
#               of offsets differs depending on the previous trial's time type.


# Load data ---------------------------------------------------------------

load('Rda/allCleanData.Rda')


# Chi-squared test --------------------------------------------------------

# variables are named as PrevType_CurrentType; for example, "NT_FT" indicates 
# that the previous trial was NT and the current trial is FT
NT_FT <- offsetData %>%
  filter(TrialTimeType == "FT" & 
           PrevTimeType == "NT" &
           Time >= -2830 &
           Time <= 5660)
FT_FT <- offsetData %>%
  filter(TrialTimeType == "FT" & 
           PrevTimeType == "FT" &
           Time >= -2830 &
           Time <= 5660)

# histograms
# hist pretty-ifys the breaks, which we don't want, so we'll set the breaks
# manually
ntBreaks <- seq(-1830, 3660, 305)
ftBreaks <- seq(-2830, 5660, 283)

histNT_FT <- hist(NT_FT$Time, breaks = ftBreaks)
histFT_FT <- hist(FT_FT$Time, breaks = ftBreaks)

chiResult <- chisq.test(x = histNT_FT$counts, y = histFT_FT$counts)
