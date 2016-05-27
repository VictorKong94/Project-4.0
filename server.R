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
    skipRow = grep("^Position", rawData) - 1
    colClass = rep("NULL", count.fields(infile, sep = "\t", skip = skipRow)[1])
    colClass[unlist(strsplit(rawData[skipRow + 1], "\t")) %in%
      c("Position", "Sample", "Detector", "Task", "Quantity")] = NA
    df = read.delim(infile, row.names = 1, stringsAsFactors = F,
                    skip = skipRow, colClasses = colClass)
    
    
    # Replace sample names with data from pertinent qPCR template file
    if (input$submitTemplate == T & !is.null(input$template)) {
      sampleNames = matrix(NA, ncol = 8, nrow = 16)
      template = read.csv(input$template$datapath, header = F, na.strings = "",
                          stringsAsFactors = F)
      I = as.matrix(expand.grid(as.numeric(rownames(template)),
                                as.numeric(gsub("V", "", colnames(template)))))
      sampleNames[I] = template[I]
      sampleNames = rep(as.vector(t(sampleNames)), each = 3)
      names(sampleNames) = paste0(rep(LETTERS[1:16], each = 24), 1:24)
      df$Sample = sampleNames[rownames(df)]
    }
    
    # Clean imported data
    df = df[!is.na(df$Quantity) & !(df$Task %in% c("NTC", "Standard")),
            names(df) %in% setdiff(names(df), c("Flag", "Task"))]
    df$Detector = factor(df$Detector)
    df$Sample = factor(df$Sample)
    
    # Identify and remove outliers. Observations with a Grubbs' Test Statistic
    # greater than or equal to 1.15 are determined to be outliers
    outlierScore = abs(df$Quantity - df$Qty.mean) / df$Qty.stddev
    names(outlierScore) = rownames(df)
    outliers = df[outlierScore >= 1.15 & !is.na(outlierScore),]
    df = df[setdiff(row.names(df), row.names(outliers)),]
    outliers = data.frame("Sample" = outliers$Sample,
                          "Detector" = outliers$Detector,
                          "Outlier Score" = outlierScore[rownames(outliers)])
    
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
    
    # Extract names of target genes
    genes = lapply(levels(df$Detector), sprintf)
    
    # Organize processed data into matrix
    if (input$sortByReplicates == T & input$repIndicator != "") {
      means$Sample = sapply(strsplit(levels(means$Sample),
                                     trimws(casefold(input$repIndicator))),
                            function(x) trimws(x[1]))
      Txs = sapply(unique(means$Target),
                   function(x) table(means$Sample[means$Target == x]))
      reps = max(Txs)
      qty = matrix(0, nrow = nrow(Txs), ncol = ncol(Txs) * reps,
                   dimnames = list(rownames(Txs),
                                   sapply(levels(means$Target), function(x)
                                     c(x, rep("", reps - 1)))))
      for (target in levels(means$Target)) {
        for (Tx in rownames(Txs)) {
          cols = seq(which(target == colnames(qty)),
                     which(target == colnames(qty)) + reps - 1)
          qty[Tx, cols] = as.numeric(means$Mean[means$Target == target &
                                                  means$Sample == Tx])
        }
      }
      qty[qty == 0] = NA
    } else {
      qty = matrix(as.numeric(means$Mean), nrow = nlevels(means$Sample),
                   dimnames = list(levels(df$Sample), levels(df$Detector)))
    }
    
    # Define list of useful processed data sets
    return(list("genes" = genes, "outliers" = outliers, "qty" = qty))
    
  })
  
  observe({
    updateRadioButtons(session, "housekeepingGene", choices = step1()$genes)
  })
  
  step2 = reactive({
    
    hkGene = input$housekeepingGene
    foldChange = step1()$qty / step1()$qty[, hkGene]
    normalFactor = step1()$qty[, hkGene] / mean(step1()$qty[, hkGene])
    normalized = step1()$qty * normalFactor
    
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
           "Raw Quantities" = step1()$qty)
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
                       "Raw Quantities" = step1()$qty),
                con)
    }
  )
  
})
