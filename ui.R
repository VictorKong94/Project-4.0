library(shiny)
library(shinythemes)

shinyUI(fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("Project 4.0"),
  sidebarLayout(
    sidebarPanel(
      
      # Upload -> SDS Output File
      fileInput(inputId = "datafile",
                label = "Select SDS Output File",
                accept = ".txt"),
      
      conditionalPanel(
        condition = "output.fileUploaded",
        
        tags$hr(),
        
        # Header for Optional Section
        tags$b("Optional (beta)"),
        
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
        
        # Upload -> qPCR Template File
        checkboxInput(inputId = "submitTemplate",
                      label = "Submit qPCR Plate Template"),
        conditionalPanel(
          condition = "input.submitTemplate == true",
          fileInput(inputId = "template",
                    label = NULL,
                    accept = ".csv")
        ),
        
        # Text Input -> Sort by Replicates
        checkboxInput(inputId = "sortByReplicates",
                      label = "Enter String to Sort by Replicates"),
        conditionalPanel(
          condition = "input.sortByReplicates == true",
          textInput(inputId = "repIndicator",
                    label = NULL),
          checkboxInput(inputId = "poolControls",
                        label = "Pool ΔCt (Control)")
        ),
        
        # Horizontal Bar Separator
        tags$hr(),
        
        # Select Input -> Select Housekeeping Gene
        selectInput(inputId = "housekeepingGene",
                    label = "Select Housekeeping Gene",
                    choices = "CLPTM"),
        
        # Radio Buttons -> Quantification Method
        radioButtons(inputId = "method",
                     label = "Select Quantification Method",
                     choices = c("Absolute (Standard Curve)" = "absolute",
                                 "Relative (ΔΔCt)" = "relative")),
        
        # Select Input -> Select Control Condition
        conditionalPanel(
          condition = "input.method == 'relative'",
          selectInput(inputId = "control",
                      label = "Select Control Condition",
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
