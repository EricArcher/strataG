output$dt.by.sample <- renderDataTable({
  df <- by.sample()
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    selection = list(mode = "multiple", target = "row")
  )
})

output$dt.by.locus <- renderDataTable({
  df <- by.locus()
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    selection = list(mode = "multiple", target = "row")
  )
})

observeEvent(input$btn.remove.samples, {
  i <- input$dt.by.sample_rows_selected
  if(!is.null(i)) {
    df <- by.sample()
    id <- df$id[i]
    all.inds <- indNames(current.g)
    to.keep <- setdiff(all.inds, id)
    if(length(to.keep) > 0) {
      current.g <<- current.g[to.keep, , ]
      reloaded$count <- reloaded$count + 1
    }
  }
})

observeEvent(input$btn.remove.loci, {
  i <- input$dt.by.locus_rows_selected
  if(!is.null(i)) {
    df <- by.locus()
    loc <- df$locus[i]
    all.loci <- locNames(current.g)
    to.keep <- setdiff(all.loci, loc)
    if(length(to.keep) > 0) {
      current.g <<- current.g[ ,to.keep , ]
      reloaded$count <- reloaded$count + 1
    }
  }
})