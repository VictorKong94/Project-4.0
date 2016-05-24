library(shiny)

shinyUI(fluidPage(
  titlePanel("Project 4.0"),
  sidebarLayout(
    sidebarPanel(
      paste("Bandwidth is limited. Please close this",
            "tab when you're finished. Thanks!"),
      tags$hr(),
      
      # Upload SDS Output File
      fileInput("datafile", "Select SDS Output File", accept = ".txt"),
      conditionalPanel(
        condition = "output.fileUploaded",
        radioButtons("housekeepingGene", "Select Housekeeping Gene",
                    choices = "CLPTM")
      ),
      
      # Upload Optional qPCR Template File
      conditionalPanel(
        condition = "output.fileUploaded",
        fileInput("template", "Select qPCR Template File (Optional)",
                  accept = ".csv")
      ),
    
      conditionalPanel(
        condition = "output.fileUploaded",
        tags$hr(),
        
        # Select Desired Output
        selectInput("outfile", "Select Data Output",
                    choices = c("Fold Change", "Normalized", "Outliers",
                                "Raw Quantities"),
                    selected = "Raw Quantities"),
        
        # Button to Download Processed Data
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
