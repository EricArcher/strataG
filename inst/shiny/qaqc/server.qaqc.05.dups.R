ui.dups <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput("sl.dups", h4("Percent of shared loci"), 0, 1, 0.8),
      textOutput("txt.dups"),
      plotOutput("plot.dups")
    ),
    mainPanel(
      dataTableOutput("dt.dups"),
      uiOutput("dup.remove")
    )
  )
}

output$txt.dups <- renderPrint({
  if(is.null(vals$gtypes)) return()
  cat("Number of loci:", ceiling(input$sl.dups * nLoc(vals$gtypes)))
})

output$dt.dups <- renderDataTable({
  df <- dups()
  if(!is.null(df)) {
    df <- df[df$prop.loci.shared >= input$sl.dups, ]
    df <- round(df, 4)
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
  )
})

output$plot.dups <- renderPlot({
  df <- dups()
  if(is.null(df)) NULL else {
    df <- data.frame(x = 1:nrow(df), y = sort(df$prop.loci.shared))
    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point() + geom_line() +
      xlab("Pairs of samples") + ylab("% Loci shared") + 
      ylim(range(c(df[, "y"], input$sl.dups))) +
      geom_hline(yintercept = input$sl.dups, color = "red")
    print(p)
  }
})

output$dup.remove <- renderUI({
  df <- dups()
  if(is.null(df)) return(NULL)
  df <- df[df$prop.loci.shared >= input$sl.dups, ]
  ids <- unique(unlist(df[, 1:2]))
  ids <- ids[order(nchar(ids), ids)]
  verticalLayout(
    hr(),
    selectizeInput(
      "dup.samples.select", "Select samples to remove",
      choices = ids, multiple = TRUE
    ),
    actionButton("btn.remove.dup.samples", "Remove selected samples")
  )
})

observeEvent(input$btn.remove.dup.samples, {
  isolate({
    if(!is.null(input$dup.samples.select)) {
      id <- input$dup.samples.select
      all.inds <- indNames(vals$gtypes)
      to.keep <- setdiff(all.inds, id)
      if(length(to.keep) > 0) {
        vals$gtypes <- vals$gtypes[to.keep, , ]
        vals$qaqc.reports$samples[id, "step.removed"] <- vals$qaqc.step
        vals$qaqc.reports$samples[id, "threshold"] <- input$sl.dups
      }
    }
  })
})