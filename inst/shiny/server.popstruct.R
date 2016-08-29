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

output$ovlResults <- renderDataTable({
  input$btn.run.pop.struct
  if(is.null(isolate(user.data$current.g))) NULL else {
    df <- overallTest(
      isolate(user.data$current.g), stats = getStats(), 
      nrep = as.numeric(input$permutation), 
      pairwise.deletion = input$pairwise.deletion
    )$result
    df <- round(df, 4)
  
    DT::datatable(
      df, rownames = FALSE,
      options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
    )
  }
})

output$pwsResults <- renderDataTable({
  input$btn.run.pop.struct
  if(is.null(isolate(user.data$current.g))) NULL else {
    df <- pairwiseTest(
      isolate(user.data$current.g), stats = getStats(), 
      nrep = as.numeric(input$permutation), 
      pairwise.deletion = input$pairwise.deletion
    )$result
    df <- round(df, 4)
    
    DT::datatable(
      df, rownames = FALSE,
      options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
    )
  }
})