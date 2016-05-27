library(shiny)

shinyUI(fluidPage(
  titlePanel("Project 4.0"),
  sidebarLayout(
    sidebarPanel(
      paste("Bandwidth is limited. Please close this",
            "tab when you're finished. Thanks!"),
      tags$hr(),
      
      # Upload -> SDS Output File
      fileInput("datafile", "Select SDS Output File", accept = ".txt"),
      
      conditionalPanel(
        condition = "output.fileUploaded",
        
        # Radio Buttons -> Select Housekeeping Gene
        radioButtons("housekeepingGene", "Select Housekeeping Gene",
                    choices = "CLPTM"),
        
        # Horizontal Bar Separator
        tags$hr(),
        
        # Header for Optional Section
        tags$b("Optional (beta)"),
        
        # Upload -> qPCR Template File
        checkboxInput("submitTemplate", "Submit qPCR Template File"),
        conditionalPanel(
          condition = "input.submitTemplate == true",
          fileInput("template", NULL, accept = ".csv")
        ),
        
        # Text Input -> Sort by Replicates
        checkboxInput("sortByReplicates", "Enter String to Sort by Replicates"),
        conditionalPanel(
          condition = "input.sortByReplicates == true",
          textInput("repIndicator", NULL)
        ),
        
        # Horizontal Bar Separator
        tags$hr(),
        
        # Select Menu -> Desired Output
        selectInput("outfile", "Select Data Output",
                    choices = c("Fold Change", "Normalized",
                                "Outliers", "Raw Quantities"),
                    selected = "Raw Quantities"),
        
        # Download Button -> Download Processed Data
        downloadButton("downloadData", "Download")
      )
      
    ),
    mainPanel(
      conditionalPanel(
        condition = "output.fileUploaded",
        tableOutput("table")
      )
    )
  )
))
