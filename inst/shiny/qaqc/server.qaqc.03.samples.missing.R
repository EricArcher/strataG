ui.samples.missing <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput("sl.samples.missing", h4("Percent of missing samples"), 0, 1, 0.05),
      textOutput("txt.samples.missing"),
      plotOutput("plot.samples.missing")
    ),
    mainPanel(
      dataTableOutput("dt.samples.missing"),
      actionButton("btn.remove.samples.missing", label = "Remove loci")
    )
  )
}

output$txt.samples.missing <- renderPrint({
  if(is.null(vals$gtypes)) return()
  num.missing <- floor(input$sl.samples.missing * nInd(vals$gtypes))
  cat("Number of samples:", num.missing)
})

df.samples.missing <- reactive({  
  df <- by.locus()
  if(!is.null(df)) {
    df <- df[, c("locus", "pct.missing", "num.missing")]
    colnames(df) <- c("Locus", "% Missing", "# Missing")
    df <- df[order(df[, 2], decreasing = TRUE), ]
    df <- df[df[, 2] >= input$sl.samples.missing, ]
    df <- round(df, 4)
  }
  df
})

output$dt.samples.missing <- renderDataTable({
  df <- df.samples.missing()
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
  )
})

output$plot.samples.missing <- renderPlot({
  df <- by.locus()
  if(is.null(df)) NULL else {
    df <- data.frame(x = 1:nrow(df), y = sort(df$pct.missing))
    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point() + geom_line() +
      xlab("Loci") + ylab("% Samples Missing") + 
      ylim(range(c(df[, "y"], input$sl.samples.missing))) +
      geom_hline(yintercept = input$sl.samples.missing, color = "red")
    print(p)
  }
})

observeEvent(input$btn.remove.samples.missing, {
  isolate({
    df <- by.locus()
    i <- which(df$pct.missing >= input$sl.samples.missing)
    if(length(i) > 0) {
      loc <- as.character(df$locus[i])
      all.loci <- locNames(vals$gtypes)
      to.keep <- setdiff(all.loci, loc)
      if(length(to.keep) > 0) {
        vals$gtypes <- vals$gtypes[ ,to.keep , ]
        vals$qaqc.reports$loci[loc, "step.removed"] <<- vals$qaqc.step
        vals$qaqc.reports$loci[loc, "threshold"] <<- input$sl.samples.missing
      }
    }
  })
})