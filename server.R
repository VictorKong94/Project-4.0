library(shiny)

shinyServer(function(input, output, session) {
  
  step1 = reactive({
    
    # Do not proceed if no data file uploaded
    if (is.null(input$datafile)) {
      return(NULL)
    }
    
    # Import plaintext output from SDS software
    options(stringsAsFactors = F)
    infile = input$datafile$datapath
    rawData = readLines(infile)
    skip = grep("^Position", rawData) - 1
    colClass = rep("NULL", count.fields(infile, sep = "\t", skip = skip)[1])
    colClass[unlist(strsplit(rawData[skip + 1], "\t")) %in%
      c("Position", "Flag", "Sample", "Detector", "Task", "Ct", "Ct.median",
        "Quantity", "Qty mean", "Qty stddev")] = NA
    df = read.delim(infile, row.names = 1, skip = skip, colClasses = colClass)
    df$X = NULL
    
    # Replace sample names with data from pertinent qPCR template file
    if (input$submitTemplate == T & !is.null(input$template)) {
      sampleNames = matrix(NA, ncol = 8, nrow = 16)
      template = read.csv(input$template$datapath, header = F, na.strings = "")
      I = as.matrix(expand.grid(as.numeric(rownames(template)),
                                as.numeric(gsub("V", "", colnames(template)))))
      sampleNames[I] = template[I]
      sampleNames = rep(as.vector(t(sampleNames)), each = 3)
      names(sampleNames) = paste0(rep(LETTERS[1:16], each = 24), 1:24)
      df$Sample = sampleNames[rownames(df)]
    }
    
    # Clean imported data
    df = df[!(df$Task %in% c("", "NTC", "Standard")),
            names(df) %in% setdiff(names(df), "Task")]
    df$Detector = factor(df$Detector)
    df$Sample = factor(df$Sample)
    
    # Identify errors, including:
    # - samples with Grubb's Outlier Score >= 1.15,
    # - missing data,
    # - samples flagged by SDS software
    df$Outlier.score = abs(df$Quantity - df$Qty.mean) / df$Qty.stddev
    errors = df[df$Flag != "Passed" | df$Outlier.score >= 1.15 |
                  is.na(df$Quantity),]
    df = df[setdiff(rownames(df), rownames(errors)),]
    errors$Reason = paste("Grubb's Outlier Score =", errors$Outlier.score)
    errors$Reason[errors$Flag == "Flagged"] = "Flagged by SDS"
    errors$Reason[is.na(errors$Quantity)] = "Missing data"
    errors$Flag = NULL
    errors$Outlier.Score = NULL
    
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
                                     trimws(input$repIndicator)),
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
    return(list("genes" = genes, "errors" = errors, "qty" = qty))
    
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
           "Errors" = step1()$errors,
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
                       "Errors" = step1()$errors,
                       "Raw Quantities" = step1()$qty),
                con)
    }
  )
  
})

