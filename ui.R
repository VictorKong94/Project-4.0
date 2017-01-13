library(shiny)
library(shinythemes)

fluidPage(
  theme = shinytheme("cerulean"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "stylesheet.css")
  ),
  titlePanel("Project 5.0"),
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
        
        # Text Input -> String to Subset Data
        checkboxInput(inputId = "subsetData",
                      label = "Enter String to Subset Data"),
        conditionalPanel(
          condition = "input.subsetData == true",
          textInput(inputId = "subsetString",
                    label = NULL)
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
        
        # Checkbox Group Input -> Select Housekeeping Gene(s)
        checkboxGroupInput(inputId = "housekeepingGenes",
                           label = "Select Housekeeping Gene(s)",
                           choices = "CLPTM",
                           selected = "CLPTM"),
        
        # Selectize Input -> Select Control Condition(s)
        selectizeInput(inputId = "control",
                         label = "Select Control Condition(s)",
                         choices = NULL,
                         multiple = TRUE),

        # Numeric Input -> Number of Treatment Condition(s)
        numericInput(inputId = "nTreatments",
                     label = "Number of Treatment Condition(s)",
                     value = 1,
                     min = 1),

        # Selectize Input(s) -> Select Treatment Condition(s)
        uiOutput(outputId = "selectTreatments"),

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
                       label = "Data",
                       class = "top-left"),
        
        # Download Button -> Download Plot
        downloadButton(outputId = "downloadPlot",
                       label = "Plot",
                       class = "top-right"),
        
        # Line Break
        tags$br(),
        
        # Download All -> Download All Processed Data
        downloadButton(outputId = "downloadAll",
                       label = "All Output",
                       class = "bottom")
      )
      
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Data", tableOutput(outputId = "tableDisplay")),
                  tabPanel("Plot", plotOutput(outputId = "plotDisplay")))
    )
  )
)
