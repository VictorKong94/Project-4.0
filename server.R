library(shiny)

source("helpers.R")

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
    # if (input$sortByReplicates & input$repIndicator != "") {
    #   means$Sample = sapply(strsplit(levels(means$Sample),
    #                                  trimws(input$repIndicator)),
    #                         function(x) trimws(x[1]))
    #   Txs = sapply(unique(means$Target),
    #                function(x) table(means$Sample[means$Target == x]))
    #   reps = max(Txs)
    #   qty = matrix(0, nrow = nrow(Txs), ncol = ncol(Txs) * reps,
    #                dimnames = list(rownames(Txs),
    #                                sapply(levels(means$Target), function(x)
    #                                  c(x, rep("", reps - 1)))))
    #   for (target in levels(means$Target)) {
    #     for (Tx in rownames(Txs)) {
    #       cols = seq(which(target == colnames(qty)),
    #                  which(target == colnames(qty)) + reps - 1)
    #       qty[Tx, cols] = as.numeric(means$Mean[means$Target == target &
    #                                               means$Sample == Tx])
    #     }
    #   }
    #   qty[qty == 0] = NA
    # } else {
      qty = matrix(as.numeric(means$Mean), nrow = nlevels(means$Sample),
                   dimnames = list(levels(df$Sample), levels(df$Detector)))
    # }
    
    # Set up choices for output shown
    if (input$method == "absolute") {
      choices = c("Errors", "Normalized", "Raw Quantities")
    } else if (input$method == "relative") {
      choices = c("Errors", "Fold Change", "Normalized", "Raw Quantities")
    }
    
    # Extract names of treatment conditions
    conditions = rownames(qty)
    
    # Extract names of target genes
    genes = levels(df$Detector)
    
    # Define set of raw quantity data
    return(list("choices" = choices,
                "conditions" = conditions,
                "errors" = errors,
                "genes" = genes,
                "qty" = qty))
    
  })
  
  observe({
    
    # Update checkbox group input for selection of housekeeping gene
    updateCheckboxGroupInput(session,
                             inputId = "housekeepingGenes",
                             choices = step1()$genes,
                             selected = step1()$genes[1])
    
    # Update checkbox group input for selection of control condition
    updateCheckboxGroupInput(session,
                             inputId = "control",
                             choices = step1()$conditions)
    
    # Update select input for selection of desired output
    updateSelectInput(session,
                      inputId = "outfile",
                      choices = step1()$choices,
                      selected = "Raw Quantities")
  
  })
  
  output$fileUploaded = reactive({
    return(!is.null(input$datafile))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = F)
  
  step2 = reactive({
    
    genes = step1()$genes
    qty = step1()$qty
    cntlCond = input$control
    hkGenes = input$housekeepingGenes
    if (length(hkGenes) == 0 | length(hkGenes) == length(genes)) {
      normalized = matrix(paste("Please select at least 1 and at most",
                                length(genes) - 1, "housekeeping gene(s)."))
      foldChange = normalized
    } else if (!input$sortByReplicates & input$method == "absolute") {
      hkCols = which(colnames(qty) %in% hkGenes)
      nf = apply(as.matrix(qty[, hkCols]), 1, geometricMean)
      nf = nf / mean(nf)
      nonHkCols = setdiff(seq(1, ncol(qty)), hkCols)
      normalized = qty
      normalized[, nonHkCols] = qty[, nonHkCols] / nf
      foldChange = normalized
    # } else if (input$sortByReplicates & input$method == "absolute") {
    #   hkCols = (which(colnames(qty) %in% hkGenes) - 1) * nrow(qty) + 1
    #   nf = matrix(nrow = nrow(qty), ncol = nrow(colnames(qty)))
    #   for (col in 1:ncol(nf)) {
    #     nf[, col] = apply(qty[, hkCols + col - 1], 1, geometricMean)
    #   }
    #   nf = nf / mean(nf)
    #   nonHkCols = setdiff(
    #     seq(1, by = ncol(nf), length.out = ncol(colnames(qty))), hkCols)
    #   normalized = qty
    #   for (col in nonHkCols) {
    #     normalized[, col:(col + ncol(nf) - 1)] =
    #       qty[, col:(col + ncol(nf) - 1)] / nf
    #   }
    #   foldChange = normalized
    } else if (!input$sortByReplicates & input$method == "relative") {
      hkCols = which(colnames(qty) %in% hkGenes)
      gMeanHkGenes = apply(as.matrix(qty[, hkCols]), 1, geometricMean)
      qty = cbind(qty, gMeanHkGenes)
      normalized = t(apply(qty, 1, function(x) x - x["gMeanHkGenes"]))
      normalized = rbind(normalized, "controlCondition" =
                           apply(normalized, 2, function(x) mean(x[cntlCond])))
      normalized = apply(normalized, 2, function(x) x - x["controlCondition"])
      foldChange = 2^(-normalized)
    # } else if (input$sortByReplicates & input$method == "relative") {
    #   hkCols = which(colnames(qty) %in% hkGenes)
    #   reps = nrow(colnames(qty))
    #   hkCt = matrix(nrow = nrow(qty), ncol = reps)
    #   for (col in 1:ncol(hkCt)) {
    #     hkCt[, col] = apply(as.matrix(qty[, hkCols + col - 1]),
    #                         1, geometricMean)
    #   }
    #   qty = cbind(qty, hkCt)
    #   colnames(qty)[length(genes) * reps + 1] = "gMeanHkGenes"
    #   normalized = rbind(qty, "controlCondition" = 
    #                      apply(normalized, 2, function(x) mean(x[cntlCond])))
    #   for (gene in 1:(length(genes) + 1)) {
    #     set = seq((gene - 1) * ncol(hkCt) + 1, by = 1, length.out = ncol(hkCt))
    #     normalized[, set] = normalized[, set] -
    #       normalized[, (ncol(normalized) - reps + 1):ncol(normalized)]
    #     normalized[, set] = normalized[, set] -
    #       mean(normalized["controlCondition", set])
    #   }
    #   foldChange = 2^(-normalized)
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
      if (input$method == "absolute") {
        name = strsplit(input$datafile$name, ".txt")[[1]]
        files = paste(name, c("(Errors).csv",
                              "(Normalized).csv",
                              "(Raw Quantities).csv"))
        tmpdir = tempdir()
        setwd(tempdir())
        write.csv(step1()$errors, files[1])
        write.csv(step2()$normalized, files[2])
        write.csv(step1()$qty, files[3])
        zip(zipfile = con, files = files)
      } else if (input$method == "relative") {
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
    }
  )
  
})

