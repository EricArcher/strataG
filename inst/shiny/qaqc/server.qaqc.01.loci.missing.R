ui.loci.missing <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput("sl.loci.missing", h4("Percent of missing loci"), 0, 1, 0.8),
      textOutput("txt.loci.missing"),
      plotOutput("plot.loci.missing")
    ),
    mainPanel(
      dataTableOutput("dt.loci.missing"),
      actionButton("btn.remove.loci.missing", label = "Remove samples")
    )
  )
}

output$txt.loci.missing <- renderPrint({
  if(is.null(isolate(vals$gtypes))) return()
  num.missing <- floor(input$sl.loci.missing * nLoc(isolate(vals$gtypes)))
  cat("Number of loci:", num.missing)
})

output$dt.loci.missing <- renderDataTable({
  df <- by.sample()
  if(!is.null(df)) {
    df <- df[, c("id", "strata", "pct.loci.missing.genotypes", "num.loci.missing.genotypes")]
    colnames(df) <- c("ID", "Strata", "% Missing", "# Missing")
    df <- df[order(df[, 3], decreasing = TRUE), ]
    df <- df[df[, 3] >= input$sl.loci.missing, ]
    df <- round(df, 4)
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
  )
})

output$plot.loci.missing <- renderPlot({
  df <- by.sample()
  if(is.null(df)) NULL else {
    df <- data.frame(x = 1:nrow(df), y = sort(df$pct.loci.missing.genotypes))
    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point() + geom_line() +
      xlab("Samples") + ylab("% Genotypes Missing") + 
      ylim(range(c(df[, "y"], input$sl.loci.missing))) +
      geom_hline(yintercept = input$sl.loci.missing, color = "red")
    print(p)
  }
})

observeEvent(input$btn.remove.loci.missing, {
  isolate({
    df <- by.sample()
    i <- which(df$pct.loci.missing.genotypes >= input$sl.loci.missing)
    if(length(i) > 0) {
      id <- as.character(df$id[i])
      all.inds <- indNames(vals$gtypes)
      to.keep <- setdiff(all.inds, id)
      if(length(to.keep) > 0) vals$gtypes <- vals$gtypes[to.keep, , ]
      vals$qaqc.reports$samples[id, "step.removed"] <- vals$qaqc.step
      vals$qaqc.reports$samples[id, "threshold"] <- input$sl.loci.missing
    }
  })
})