library(ggplot2)
library(shiny)

source("helpers.R")

function(input, output, session) {
  
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
    
    # Extract a subset of the data according to user preferences
    if (input$subsetData & !is.null(input$subsetString)) {
      subsetData = grep(input$subsetString, df$Sample)
      if (length(subsetData) > 0) df = df[subsetData,]
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
    means = t(mapply(function(x, y)
      c(x, y, mean(df$X[df$Detector == x & df$Sample == y], na.rm = T)),
      rep(levels(df$Detector), each = nlevels(df$Sample)),
      rep(levels(df$Sample), times = nlevels(df$Detector))
    ))
    means = data.frame("Target" = as.factor(means[, 1]),
                       "Sample" = as.factor(means[, 2]),
                       "Mean" = as.numeric(means[, 3]))
    
    qty = matrix(as.numeric(means$Mean), nrow = nlevels(means$Sample),
                 dimnames = list(levels(df$Sample), levels(df$Detector)))
    
    # Set up choices for output shown
    if (input$method == "absolute") {
      choices = c("Errors",
                  "Normalized Quantities",
                  "Raw Quantities")
    } else if (input$method == "relative") {
      choices = c("ΔΔCt Values",
                  "Errors",
                  "Fold Changes",
                  "Raw Ct Values")
    }
    
    # Extract names of treatment conditions
    conditions = rownames(qty)
    
    # Extract names of target genes
    genes = levels(df$Detector)
    
    # Initialize name to save output file as
    if (grepl("export", input$datafile$name)) splt = "export" else splt = ".txt"
    if (input$subsetData & !is.null(input$subsetString)) {
      name = paste0(trimws(strsplit(input$datafile$name, splt)[[1]][1]),
                    " -- ", input$subsetString)
    } else {
      name = paste0(trimws(strsplit(input$datafile$name, splt)[[1]][1]))
    }
    
    # Define set of raw quantity data
    return(list("choices" = choices,
                "conditions" = conditions,
                "errors" = errors,
                "genes" = genes,
                "name" = name,
                "qty" = qty))
    
  })
  
  observe({
    
    # Update checkbox group input for selection of housekeeping gene
    updateCheckboxGroupInput(session,
                             inputId = "housekeepingGenes",
                             choices = step1()$genes,
                             selected = step1()$genes[1])
    
    # Update selectize input for selection of control condition
    updateSelectizeInput(session,
                         inputId = "control",
                         label = "Select Control Condition(s)",
                         choices = step1()$conditions,
                         server = TRUE)
    
    # Update select input for selection of desired output
    updateSelectInput(session,
                      inputId = "outfile",
                      choices = step1()$choices,
                      selected = switch(input$method,
                                        "absolute" = "Raw Quantities",
                                        "relative" = "Raw Ct Values"))

  })

  output$fileUploaded = reactive({
    return(!is.null(input$datafile))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = F)
  
  step2 = reactive({
    
    genes = step1()$genes
    qty = step1()$qty
    ctrlCond = input$control
    hkGenes = input$housekeepingGenes
    if (length(hkGenes) == 0 | length(hkGenes) == length(genes)) {
      normalized = matrix(paste("Please select at least 1 and at most",
                                length(genes) - 1, "housekeeping gene(s)."))
      foldChange = normalized
    } else if (input$method == "absolute") {
      hkCols = which(colnames(qty) %in% hkGenes)
      nf = apply(as.matrix(qty[, hkCols]), 1, geometricMean)
      nf = nf / mean(nf, na.rm = T)
      nonHkCols = setdiff(seq(1, ncol(qty)), hkCols)
      normalized = qty
      normalized[, nonHkCols] = qty[, nonHkCols] / nf
      foldChange = data.frame()
    } else if (input$method == "relative") {
      hkCols = which(colnames(qty) %in% hkGenes)
      gMeanHkGenes = apply(as.matrix(qty[, hkCols]), 1, geometricMean)
      qty = cbind(qty, gMeanHkGenes)
      normalized = t(apply(qty, 1, function(x) x - x["gMeanHkGenes"]))
      controlCondition = apply(normalized, 2, function(x)
        mean(x[ctrlCond], na.rm = T))
      normalized = rbind(normalized, controlCondition)
      normalized = apply(normalized, 2, function(x) x - x["controlCondition"])
      normalized[abs(normalized) < 0.0000001] = 0
      if (length(hkGenes) == 1) {
        hkMean = geometricMean(qty[, hkCols])
        hkCol = qty[, hkCols] - hkMean
        normalized[, hkCols] = c(hkCol, 0)
      }
      foldChange = 2^(-normalized)
    }
    
    # Define list of useful processed data sets
    return(list("foldChange" = foldChange,
                "normalized" = normalized))
    
  })
  
  output$selectTreatments = renderUI({
    nTreatments = as.integer(input$nTreatments)
    lapply(1:(2 * nTreatments), function(i) {
      if (i %% 2 == 1) {
        j = ceiling(i / 2)
        textInput(inputId = paste0("trt", j, "name"),
                  label = paste("Enter Name for Treatment", j),
                  value = paste("Treatment", j))
      } else {
        j = i / 2
        selectizeInput(inputId = paste0("trt", j),
                       label = paste("Select Treatment", j, "Condition(s)"),
                       choices = step1()$conditions,
                       multiple = TRUE)
      }
    })
  })

  step3 = reactive({

    # Generate data to produce bar graph of results
    ctrlCond = input$control
    if (input$method == "absolute") {
      df = step2()$normalized
      var = "Quantity"
    } else {
      df = step2()$foldChange
      var = "Fold Change"
    }
    myData = data.frame("signal" = as.vector(df[ctrlCond, -ncol(df)]),
                        "condition" = "Control",
                        "detector" = rep(colnames(df)[-ncol(df)],
                                         each = length(ctrlCond)))
    for (i in 1:input$nTreatments) {
      treatments = input[[paste0("trt", i)]]
      myData = rbind(myData,
                     data.frame("signal" = as.vector(df[treatments, -ncol(df)]),
                                "condition" = input[[paste0("trt", i, "name")]],
                                "detector" = rep(colnames(df)[-ncol(df)],
                                                 each = length(treatments))))
    }
    plotData = aggregate(myData$signal,
                         by = list("condition" = myData$condition,
                                   "detector" = myData$detector),
                         FUN = function(x) c(mean = mean(x), sd = sd(x),
                                             n = length(x)))
    plotData = do.call(data.frame, plotData)
    plotData$se = plotData$x.sd / sqrt(plotData$x.n)
    colnames(plotData) = c("condition", "detector", "mean", "sd", "n", "se")
    dodge = position_dodge(width = 0.9)
    limits = aes(ymax = plotData$mean + plotData$se,
                 ymin = plotData$mean - plotData$se)
    p = ggplot(data = plotData, aes(x = factor(condition), y = mean,
                                    fill = factor(detector))) +
      geom_bar(stat = "identity",
               position = position_dodge(0.9)) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Condition", y = var) +
      scale_fill_discrete(name = "Detector")
    return(p)
    
  })
  
  output$tableDisplay = renderTable({
    switch(input$outfile,
           "ΔΔCt Values" = step2()$normalized,
           "Errors" = step1()$errors,
           "Fold Changes" = step2()$foldChange,
           "Normalized Quantities" = step2()$normalized,
           "Raw Ct Values" = step1()$qty,
           "Raw Quantities" = step1()$qty)
  }, include.rownames = T)
  
  output$plotDisplay = renderPlot({
    print(step3())
  })
  
  output$downloadData = downloadHandler(
    filename = function(con) paste0(step1()$name, " (", input$outfile, ").csv"),
    content = function(con) {
      write.csv(switch(input$outfile,
                       "ΔΔCt Values" = step2()$normalized,
                       "Errors" = step1()$errors,
                       "Fold Changes" = step2()$foldChange,
                       "Normalized Quantities" = step2()$normalized,
                       "Raw Ct Values" = step1()$qty,
                       "Raw Quantities" = step1()$qty),
                con)
    }
  )
  
  output$downloadPlot = downloadHandler(
    filename = function(con) paste(step1()$name, "(Plot).pdf"),
    content = function(con) ggsave(con, step3())
  )
  
  output$downloadAll = downloadHandler(
    filename = function(con) paste0(step1()$name, ".zip"),
    content = function(con) {
      if (input$method == "absolute") {
        files = paste(step1()$name, c("(Errors).csv",
                                      "(Normalized Quantities).csv",
                                      "(Raw Quantities).csv"))
        tmpdir = tempdir()
        setwd(tempdir())
        write.csv(step1()$errors, files[1])
        write.csv(step2()$normalized, files[2])
        ggsave(files[3], step3())
        write.csv(step1()$qty, files[4])
        zip(zipfile = con, files = files)
      } else {
        files = paste(step1()$name, c("(Errors).csv",
                                      "(Fold Changes).csv",
                                      "(ΔΔCt Values).csv",
                                      "(Plot).pdf",
                                      "(Raw Ct Values).csv"))
        tmpdir = tempdir()
        setwd(tempdir())
        write.csv(step1()$errors, files[1])
        write.csv(step2()$foldChange, files[2])
        write.csv(step2()$normalized, files[3])
        ggsave(files[4], step3())
        write.csv(step1()$qty, files[5])
        zip(zipfile = con, files = files)
      }
    }
  )
  
}
