library(shiny)
library(shinythemes)

shinyUI(fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("Project 4.0 (beta)"),
  sidebarLayout(
    sidebarPanel(
      
      # Upload -> SDS Export File
      fileInput(inputId = "datafile",
                label = "Select SDS Export File",
                accept = ".txt"),
      
      conditionalPanel(
        condition = "output.fileUploaded",
        
        # Horizontal Bar Separator
        tags$hr(),
        
        # Header for Optional Section
        tags$b("Options"),
        
        # Upload -> qPCR Template File
        checkboxInput(inputId = "submitTemplate",
                      label = "Submit qPCR Plate Template"),
        conditionalPanel(
          condition = "input.submitTemplate == true",
          fileInput(inputId = "template",
                    label = NULL,
                    accept = ".csv")
        ),
        
        # Checkbox -> Remove Samples Flagged by SDS
        checkboxInput(inputId = "rmFlag",
                      label = "Remove Samples Flagged by SDS",
                      value = T),
        
        # Numeric Input -> Grubb's Outlier Score Cutoff
        checkboxInput(inputId = "rmOutliers",
                      label = "Remove Outliers",
                      value = T),
        conditionalPanel(
          condition = "input.rmOutliers == true",
          numericInput(inputId = "outCutoff",
                       label = "Grubb's Outlier Score Cutoff",
                       value = 1.15)
        ),
        
        # Horizontal Bar Separator
        tags$hr(),
        
        # Radio Buttons -> Quantification Method
        radioButtons(inputId = "method",
                     label = "Select Quantification Method",
                     choices = c("Absolute (Standard Curve)" = "absolute",
                                 "Relative (ΔΔCt)" = "relative")),
        
        # Checkbox Group Input -> Select Housekeeping Gene
        checkboxGroupInput(inputId = "housekeepingGenes",
                           label = "Select Housekeeping Gene(s)",
                           choices = "CLPTM",
                           selected = "CLPTM"),
        
        # Checkbox Group Input -> Select Control Condition
        conditionalPanel(
          condition = "input.method == 'relative'",
          checkboxGroupInput(inputId = "control",
                             label = "Select Control Condition(s)",
                             choices = "Control")
        ),
        
        # Horizontal Bar Separator
        tags$hr(),
        
        # Select Menu -> Desired Output
        selectInput(input = "outfile",
                    label = "Select Data Output",
                    choices = c("Errors",
                                "Fold Change",
                                "Normalized",
                                "Raw Quantities"),
                    selected = "Raw Quantities"),
        
        # Download Button -> Download Processed Data
        downloadButton(outputId = "downloadData",
                       label = "Download"),
        
        # Download All -> Download All Processed Data
        downloadButton(outputId = "downloadAll",
                       label = "Download All")
      )
      
    ),
    mainPanel(
      conditionalPanel(
        condition = "output.fileUploaded",
        tableOutput(outputId = "table")
      )
    )
  )
))
