source("server.popstruct.ui.R", local = TRUE)

getStats <- function() {
  stats <- if(ploidy(vals$gtypes) == 1) {
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
  if(is.null(vals$gtypes)) NULL else {
    overallTest(
      vals$gtypes, stats = getStats(), 
      nrep = as.numeric(input$nrep), 
      pairwise.deletion = input$pairwise.deletion,
      model = input$subModel,
      quietly = TRUE
    )$result
  }
})

pws.df <- reactive({
  if(is.null(vals$gtypes)) NULL else {
    pairwiseTest(
      vals$gtypes, stats = getStats(), 
      nrep = as.numeric(input$nrep), 
      pairwise.deletion = input$pairwise.deletion,
      model = input$subModel,
      quietly = TRUE
    )$result
  }
})

observeEvent(input$run.popstruct, {
  output$ovlResults <- renderDataTable({
    if(is.null(vals$gtypes) | !input$ovl) NULL else {
      df <- round(ovl.df(), 4)
      DT::datatable(
        df, rownames = TRUE,
        options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
      )
    }
  })
  
  output$pwsResults <- renderDataTable({
    if(is.null(vals$gtypes) | !input$ovl) NULL else {
      df <- round(pws.df(), 4)
      DT::datatable(
        df, rownames = FALSE,
        options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
      )
    }
  })
})

observeEvent(input$savePopStructResults, {
  label <- make.names(description(vals$gtypes))
  output.dir <- paste0(label, "_PopStructResults")
  output.dir <- file.path(vals$wd, output.dir)
  if(!dir.exists(output.dir)) dir.create(output.dir)
  
  if(input$ovl) {
    df <- ovl.df()
    if(!is.null(df)) {
      fname <- file.path(output.dir, paste0(label, "_overall.results.csv"))
      write.csv(df, file = fname)
    }
  }
  
  if(input$pws) {
    df <- pws.df()
    if(!is.null(df)) {
      fname <- file.path(output.dir, paste0(label, "_pairwise.results.csv"))
      write.csv(df, file = fname, row.names = FALSE)
    }
  }
})