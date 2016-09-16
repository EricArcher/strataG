smry <- reactive({
  if(is.null(vals$gtypes)) NULL else summary(vals$gtypes)
})

output$seqSmryTable <- renderDataTable({
  x <- smry()
  if(is.null(x)) NULL else {
    df <- round(data.frame(x$strata.smry[, -2]), 4)
    df <- cbind(strata = rownames(df), df)
    DT::datatable(
      df, rownames = FALSE, selection = "none",
      options = list(paging = FALSE, searching = FALSE, scrollX = TRUE)
    )
  }
})

output$hapFreqTable <- renderDataTable({
  x <- smry()
  if(is.null(x)) NULL else {
    df <- data.frame(x$allele.freqs[[1]][, "freq", ])
    df <- cbind(haplotype = rownames(df), df)
    DT::datatable(
      df, rownames = FALSE, selection = "none",
      options = list(searching = FALSE, scrollX = TRUE)
    )
  }
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

output$lowFreqSubTable <- renderDataTable({
  if(is.null(vals$gtypes)) NULL else {
    lowFreqSubs(
      getSequences(sequences(vals$gtypes), 1), 
      min.freq = input$minFreq, motif.length = input$motifLength
    )
  }
})
  