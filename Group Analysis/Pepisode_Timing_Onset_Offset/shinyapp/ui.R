library(shiny)
library(dplyr)
library(ggplot2)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Pepisode Timing"),
  
  sidebarPanel(
    selectInput("freqBand", "Frequency Band:",
                list("Delta",
                     "Theta",
                     "Alpha",
                     "Beta",
                     "Gamma")),
    selectInput("variable", "Variable:",
                list("Episodes" = "Episode",
                     "Onsets" = "Onset",
                     "Offsets" = "Offset")),
    selectInput("trialType", "Trial Type:",
                list("NT", "FT"))
    ),
  
  mainPanel(
    plotOutput("plot")
    )
))