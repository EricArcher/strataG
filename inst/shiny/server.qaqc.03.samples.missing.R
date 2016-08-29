ui.samples.missing <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "sl.samples.missing", label = h4("3) Percent of missing samples"),
        min = 0, max = 1, value = 0.05
      ),
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
  if(is.null(user.data$current.g)) return()
  cat("Number of samples:", input$sl.samples.missing * nInd(user.data$current.g))
})

output$dt.samples.missing <- renderDataTable({
  df <- by.locus()
  if(!is.null(df)) {
    df <- df[, c("locus", "pct.missing", "num.missing")]
    colnames(df) <- c("Locus", "% Missing", "# Missing")
    df <- df[order(df[, 2], decreasing = TRUE), ]
    df <- df[df[, 2] >= input$sl.samples.missing, ]
    df <- round(df, 4)
  }
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
  df <- by.locus()
  i <- which(df$pct.missing >= input$sl.samples.missing)
  if(length(i) > 0) {
    loc <- df$locus[i]
    all.loci <- locNames(user.data$current.g)
    to.keep <- setdiff(all.loci, loc)
    if(length(to.keep) > 0) {
      user.data$current.g <- user.data$current.g[ ,to.keep , ]
    }
  }
})