library(shiny)

shinyServer(function(input, output, session) {
  
  step1 = reactive({
    
    # Do not proceed if no data file uploaded
    if (is.null(input$datafile)) {
      return(NULL)
    }
    
    # Import plaintext output from SDS software
    infile = input$datafile$datapath
    rawData = readLines(infile)
    xRows = grep("^Position", rawData) - 1
    oCols = grep("Qty stddev", unlist(strsplit(rawData[xRows + 1], "\t")))
    xCols = count.fields(infile, sep = "\t", skip = xRows)[1] - oCols
    df = read.delim(infile, row.names = 1, stringsAsFactors = F, skip = xRows,
                    colClasses = c(rep(NA, oCols), rep("NULL", xCols)))
    
    # Replace sample names with data from pertinent qPCR template file
    if (!is.null(input$template)) {
      sampleNames = read.csv(input$template$datapath, header = F,
                             stringsAsFactors = F)
      df$Sample = rep(as.vector(t(sampleNames)), each = 3)
    }
    
    # Clean imported data
    df = df[df$Flag == "Passed" & !(df$Task %in% c("NTC", "Standard")),
            names(df) %in% setdiff(names(df), c("Flag", "Task"))]
    df$Detector = factor(df$Detector)
    df$Sample = factor(df$Sample)
    
    # Identify and remove outliers. Observations with a Grubbs' Test Statistic
    # greater than or equal to 1.15 are determined to be outliers
    outlierScore = abs(df$Quantity - df$Qty.mean) / df$Qty.stddev
    outliers = df[outlierScore >= 1.15,]
    df = df[setdiff(row.names(df), row.names(outliers)),]
    outliers = data.frame("Sample" = outliers$Sample,
                          "Detector" = outliers$Detector,
                          "Outlier Score" = outlierScore[outlierScore >= 1.15])
    
    # Compute means of remaining observations
    means = data.frame("Target" = factor(levels = levels(df$Detector)),
                       "Sample" = factor(levels = levels(df$Sample)),
                       "Mean" = numeric())
    for (target in levels(df$Detector)) {
      for (sample in levels(df$Sample)) {
        means[nrow(means) + 1,] = c(target, sample,
                                    mean(df$Quantity[df$Detector == target &
                                                       df$Sample == sample]))
      }
    }
    
    # Organize processed data into matrix
    quantity = matrix(as.numeric(means$Mean), nrow = nlevels(means$Sample),
                      dimnames = list(levels(df$Sample), levels(df$Detector)))
    
    # Extract names of target genes
    genes = list()
    for (gene in colnames(quantity)) {
      genes[[gene]] = gene
    }
    
    # Define list of useful processed data sets
    return(list("genes" = genes, "outliers" = outliers, "quantity" = quantity))
    
  })
  
  observe({
    updateRadioButtons(session, "housekeepingGene", choices = step1()$genes)
  })
  
  step2 = reactive({
    
    housekeepingGene = input$housekeepingGene
    foldChange = step1()$quantity / step1()$quantity[, housekeepingGene]
    normalized = step1()$quantity * foldChange
    
    # Define list of useful processed data sets
    return(list("foldChange" = foldChange, "normalized" = normalized))
    
  })  
  
  output$fileUploaded = reactive({
    return(!is.null(input$datafile))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = F)
  
  output$table = renderTable({
    switch(input$outfile,
           "Fold Change" = step2()$foldChange,
           "Normalized" = step2()$normalized,
           "Outliers" = step1()$outliers,
           "Raw Quantities" = step1()$quantity)
  })
  
  output$downloadData = downloadHandler(
    filename = function() {
      name = strsplit(input$datafile$name, ".txt")[[1]]
      paste0(name, " (", input$outfile, ").csv")
    },
    content = function(con) {
      write.csv(switch(input$outfile,
                       "Fold Change" = step2()$foldChange,
                       "Normalized" = step2()$normalized,
                       "Outliers" = step1()$outliers,
                       "Raw Quantities" = step1()$quantity),
                con)
    }
  )
  
})
