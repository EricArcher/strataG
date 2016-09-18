save.fname <- function(x, vals) {
  if(is.null(vals$gtypes)) return(NULL)
  label <- make.names(description(vals$gtypes))
  output.dir <- paste0(label, "_", x)
  output.dir <- file.path(vals$wd, output.dir)
  if(!dir.exists(output.dir)) dir.create(output.dir)
  file.path(output.dir, paste0(label, "_"))
}

smry <- reactive({
  if(is.null(vals$gtypes)) NULL else summary(vals$gtypes)
})

seq.smry.df <- reactive({
  x <- smry()
  if(is.null(x)) NULL else {
    df <- round(data.frame(x$strata.smry[, -2]), 4)
    cbind(strata = rownames(df), df)
  }
})

hap.freq.df <- reactive({
  x <- smry()
  if(is.null(x)) NULL else {
    df <- data.frame(x$allele.freqs[[1]][, "freq", ])
    cbind(haplotype = rownames(df), df)
  }
})



output$seqSmryTable <- renderDataTable({
  df <- seq.smry.df()
  DT::datatable(
    df, rownames = FALSE, selection = "none",
    options = list(paging = FALSE, searching = FALSE, scrollX = TRUE)
  )
})

observeEvent(input$saveSeqSmry, {
  isolate({
    df <- seq.smry.df()
    if(!is.null(df)) {
      fname <- paste0(save.fname("QAQC", vals), "seq.smry.csv")
      write.csv(df, file = fname, row.names = FALSE)
      showNotification("Sequence summaries saved", duration = 2, type = "message")
    }
  })
})


output$hapFreqTable <- renderDataTable({
  df <- hap.freq.df()
  DT::datatable(
    df, rownames = FALSE, selection = "none",
    options = list(searching = FALSE, scrollX = TRUE)
  )
})

observeEvent(input$saveHapFreq, {
  isolate({
    df <- hap.freq.df()
    if(!is.null(df)) {
      fname <- paste0(save.fname("QAQC", vals), "hap.freq.csv")
      write.csv(df, file = fname, row.names = FALSE)
      showNotification("Haplotype frequencies saved", duration = 2, type = "message")
    }
  })
})



output$seqLikePlot <- renderPlot({
  if(is.null(vals$gtypes)) NULL else {
    sequenceLikelihoods(
      getSequences(sequences(vals$gtypes), 1), model = input$smrySubModel, 
      pairwise.deletion = input$smryPwsDelete,
      n = input$smrySeqLikeN
    )
  }
})

observeEvent(input$saveSeqLike, {
  if(!is.null(vals$gtypes)) {
    isolate({
      df <- sequenceLikelihoods(
        getSequences(sequences(vals$gtypes), 1), model = input$smrySubModel, 
        pairwise.deletion = input$smryPwsDelete,
        n = 0
      )
      if(!is.null(df)) {
        fname <- paste0(save.fname("QAQC", vals), "seq.likelihoods.csv")
        write.csv(df, file = fname, row.names = FALSE)
        showNotification("Sequence likelihoods saved", duration = 2, type = "message")
      }
    })
  }
})



output$lowFreqSubTable <- renderDataTable({
  if(is.null(vals$gtypes)) NULL else {
    lowFreqSubs(
      getSequences(sequences(vals$gtypes), 1), 
      min.freq = input$minFreq, motif.length = input$motifLength
    )
  }
})

observeEvent(input$saveLowFreqSubs, {
  if(!is.null(vals$gtypes)) {
    isolate({
      df <- lowFreqSubs(
        getSequences(sequences(vals$gtypes), 1), 
        min.freq = input$minFreq, motif.length = input$motifLength
      )
      if(!is.null(df)) {
        fname <- paste0(save.fname("QAQC", vals), "low.freq.subs.csv")
        write.csv(df, file = fname, row.names = FALSE)
        showNotification("Low frequency substitutions saved", duration = 2, type = "message")
      }
    })
  }
})
