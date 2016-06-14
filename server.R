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
    colClass[unlist(strsplit(rawData[skip + 1], "\t")) %in% c(
      "Position", "Flag", "Sample", "Detector", "Task", "Ct", "Quantity")] = NA
    df = read.delim(infile, row.names = 1, skip = skip, colClasses = colClass,
                    na.strings = c(NA, "", "Undetermined"))
    df$X = NULL
    
    # Replace sample names with data from pertinent qPCR template file
    if (input$submitTemplate & !is.null(input$template)) {
      sampleNames = matrix(NA, ncol = 8, nrow = 16)
      template = read.csv(input$template$datapath, header = F, na.strings = "")
      I = as.matrix(expand.grid(as.numeric(rownames(template)),
                                as.numeric(gsub("V", "", colnames(template)))))
      sampleNames[I] = template[I]
      sampleNames = rep(as.vector(t(sampleNames)), each = 3)
      names(sampleNames) = paste0(rep(LETTERS[1:16], each = 24), 1:24)
      sampleNames[is.na(sampleNames)] = names(sampleNames[is.na(sampleNames)])
      df$Sample = sampleNames[rownames(df)]
    }
    
    # Clean imported data
    df = df[!(df$Task %in% c(NA, "NTC", "Standard")),
            names(df) %in% setdiff(names(df), "Task")]
    df$Detector = factor(df$Detector)
    df$Sample = factor(df$Sample)
    
    # Define variable of interest according to quantification method
    if (input$method == "absolute") df$X = df$Quantity else df$X = df$Ct
    
    # Identify errors, including:
    # - samples with Grubb's Outlier Score >= 1.15,
    # - missing data,
    # - samples flagged by SDS software
    df$Outlier.score = 0
    for (target in levels(df$Detector)) {
      for (sample in levels(df$Sample)) {
        i = df$Detector == target & df$Sample == sample
        df$Outlier.score[i] = abs((df$X[i] - mean(df$X[i])) / sd(df$X[i]))
      }
    }
    if (input$rmFlag) rmFlag = "Passed" else rmFlag = c("Passed", "Flagged")
    if (input$rmOutliers) {
      outCutoff = input$outCutoff
    } else {
      outCutoff = max(df$Outlier.score, na.rm = T) + 1
    }
    errors = df[!(df$Flag %in% rmFlag) |
                  !(df$Outlier.score < outCutoff) |
                  is.na(df$X),]
    df = df[setdiff(rownames(df), rownames(errors)),]
    if (nrow(errors) != 0) {
      errors$Reason = paste("Grubb's Outlier Score =", errors$Outlier.score)
    }
    if (input$rmFlag) errors$Reason[errors$Flag == "Flagged"] = "Flagged by SDS"
    errors$Reason[is.na(errors$X)] = "Missing data"
    errors = errors[, setdiff(colnames(errors), c("X", "Outlier.score"))]
    
    # Compute means 
    means = data.frame("Target" = factor(levels = levels(df$Detector)),
                       "Sample" = factor(levels = levels(df$Sample)),
                       "Mean" = numeric())
    for (target in levels(df$Detector)) {
      for (sample in levels(df$Sample)) {
        i = df$Detector == target & df$Sample == sample
        means[nrow(means) + 1,] = c(target, sample, mean(df$X[i]))
      }
    }
    
    # Organize processed data into matrix
    if (input$sortByReplicates & input$repIndicator != "") {
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
    
    # Extract names of treatment conditions
    conditions = lapply(rownames(qty), sprintf)
    
    # Extract names of target genes
    genes = lapply(levels(df$Detector), sprintf)
    
    # Define set of raw quantity data
    return(list("conditions" = conditions,
                "errors" = errors,
                "genes" = genes,
                "qty" = qty))
    
  })
  
  observe({
    
    # Update radio buttons for selection of housekeeping gene
    updateSelectInput(session,
                      inputId = "housekeepingGene",
                      choices = step1()$genes)
    
    # Update radio buttons for selection of control condition
    updateSelectInput(session,
                      inputId = "control",
                      choices = step1()$conditions)
  
  })
  
  output$fileUploaded = reactive({
    return(!is.null(input$datafile))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = F)
  
  step2 = reactive({
    
    qty = step1()$qty
    cntlCond = input$control
    hkGene = input$housekeepingGene
    if (input$sortByReplicates) {
      hkGene = qty[, seq(which(colnames(colnames(qty)) == hkGene), by = 1,
                         length.out = nrow(colnames(qty)))]
      hkGeneSet = matrix(rep(hkGene, ncol(colnames(qty))),
                         nrow = nrow(qty))
      if (input$method == "absolute") {
        foldChange = qty / hkGeneSet
        normalFactor = hkGeneSet / mean(hkGene)
        normalized = qty / normalFactor
      } else {
        normalized = qty - hkGeneSet
        normalized = apply(normalized, 2, function(x) x - x[cntlCond])
        foldChange = 2^(-normalized)
      }
    } else {
      if (input$method == "absolute") {
        foldChange = qty / qty[, hkGene]
        normalFactor = qty[, hkGene] / mean(qty[, hkGene])
        normalized = qty / normalFactor
      } else {
        normalized = t(apply(qty, 1, function(x) x - x[hkGene]))
        normalized = apply(normalized, 2, function(x) x - x[cntlCond])
        foldChange = 2^(-normalized)
      }
    }
    
    # Define list of useful processed data sets
    return(list("foldChange" = foldChange,
                "normalized" = normalized))
    
  })  
  
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
                       "Errors" = step1()$errors,
                       "Fold Change" = step2()$foldChange,
                       "Normalized" = step2()$normalized,
                       "Raw Quantities" = step1()$qty),
                con)
    }
  )
  
  output$downloadAll = downloadHandler(
    filename = function() {
      name = paste0(strsplit(input$datafile$name, ".txt")[[1]], ".zip")
    },
    content = function(con) {
      name = strsplit(input$datafile$name, ".txt")[[1]]
      files = paste(name, c("(Errors).csv",
                            "(Fold Change).csv",
                            "(Normalized).csv",
                            "(Raw Quantities).csv"))
      tmpdir = tempdir()
      setwd(tempdir())
      write.csv(step1()$errors, files[1])
      write.csv(step2()$foldChange, files[2])
      write.csv(step2()$normalized, files[3])
      write.csv(step1()$qty, files[4])
      zip(zipfile = con, files = files)
    }
  )
  
})

