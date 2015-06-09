library(shiny)
library(dplyr)
library(ggplot2)

# Load data
load('~/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Timing_Onset_Offset/Rda/allAnalyzedData.Rda')

# Define server logic required to plot various variables
shinyServer(function(input, output) {
  
  # Select data for this plot
  plotData <- reactive({
    if (input$variable == "Episode") {
      episodeData %>%
        filter(FrequencyBand == input$freqBand &
                 TrialTimeType == input$trialType)
    } else {
      onOffData %>%
        ungroup() %>%
        filter(ObservationType == input$variable & 
                 TrialTimeType == input$trialType & 
                 FrequencyBand == input$freqBand)
    }
    
  })
  
  plotTitle <- reactive({
    paste(input$freqBand, input$trialType, input$variable)
  })
  
  plotLines <- reactive({
    if (input$trialType == "NT") {
      c(0, 1830)
    } else{
      c(0, 2830)
    }
  })
  
  # Generate plot
  output$plot <- renderPlot({
    if (input$variable == "Episode") {
      p <- ggplot(plotData(), aes(x = Time, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM)) + 
        geom_ribbon(color = "lightskyblue", fill = "lightskyblue") +
        geom_line(color = "steelblue4") +
        geom_vline(xintercept = plotLines()[1]) +
        geom_vline(xintercept = plotLines()[2]) +
        facet_wrap(~ElectrodeID) +
        theme_few() +
        labs(y = "Mean Pepisode", title = plotTitle()) +
        theme(text = element_text(size = 24),
              axis.text = element_text(size = 18),
              axis.title.x = element_text(vjust = -0.5),
              axis.title.y = element_text(vjust = 1),
              panel.margin = unit(1, "lines"))
      print(p)
    } else {
      p <- ggplot(plotData(), aes(x = Time)) +
        geom_histogram(colour = "steelblue4", fill = "lightskyblue") +
        facet_wrap(~ElectrodeID, scales = "free_y") +
        geom_vline(xintercept = plotLines()[1]) +
        geom_vline(xintercept = plotLines()[2]) +
        theme_few() +
        theme(text = element_text(size = 24),
              panel.margin = unit(1, "lines")) +
        labs(y = "# of Episodes", title = plotTitle())
      print(p)
    }
  })
  
})