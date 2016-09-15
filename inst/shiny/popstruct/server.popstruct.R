source("server.popstruct.ui.R", local = TRUE)

getStats <- function() {
  stats <- if(input$ploidy == "1") {
    choice <- c(input$chi2, input$fst, input$phist)
    list(statChi2, statFst, statPhist)[choice]
  } else {
    choice <- c(
      input$chi2, input$fst, input$fst.prime, input$gst, 
      input$gst.prime, input$gst.dbl.prime, input$jost.d, input$fis
    )
    list(
      statChi2, statFst, statFstPrime, statGst, 
      statGstPrime, statGstDblPrime, statJostD, statFis
    )[choice]
  }
}

ovl.df <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    overallTest(
      user.data$current.g, stats = getStats(), 
      nrep = as.numeric(input$nrep), 
      pairwise.deletion = input$pairwise.deletion,
      quietly = TRUE
    )$result
  }
})

pws.df <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    pairwiseTest(
      user.data$current.g, stats = getStats(), 
      nrep = as.numeric(input$nrep), 
      pairwise.deletion = input$pairwise.deletion,
      quietly = TRUE
    )$result
  }
})

observeEvent(input$run.popstruct, {
  output$ovlResults <- renderDataTable({
    if(is.null(user.data$current.g) & input$ovl) NULL else {
      df <- round(ovl.df(), 4)
      DT::datatable(
        df, rownames = TRUE,
        options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
      )
    }
  })
  
  output$pwsResults <- renderDataTable({
    if(is.null(user.data$current.g) & input$ovl) NULL else {
      df <- round(pws.df(), 4)
      DT::datatable(
        df, rownames = FALSE,
        options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
      )
    }
  })
})