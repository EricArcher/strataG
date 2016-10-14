ui.samples.hmzg <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput("sl.samples.hmzg", h4("Percent of homozygous samples"), 0, 1, 0.8),
      textOutput("txt.samples.hmzg"),
      plotOutput("plot.samples.hmzg")
    ),
    mainPanel(
      dataTableOutput("dt.samples.hmzg"),
      actionButton("btn.remove.samples.hmzg", label = "Remove loci")
    )
  )
}

output$txt.samples.hmzg <- renderPrint({
  if(is.null(vals$gtypes)) return()
  num.hmzg <- floor(input$sl.samples.hmzg * nInd(vals$gtypes))
  cat("Number of samples:", num.hmzg)
})

df.samples.hmzg <- reactive({
  df <- by.locus()
  if(!is.null(df)) {
    df <- df[, c("locus", "pct.hmzg", "num.hmzg")]
    colnames(df) <- c("Locus", "% Homozygous", "# Homozygous")
    df <- df[order(df[, 2], decreasing = TRUE), ]
    df <- df[df[, 2] >= input$sl.samples.hmzg, ]
    df <- round(df, 4)
  }
  df
})

output$dt.samples.hmzg <- renderDataTable({
  df <- df.samples.hmzg()
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
  )
})

output$plot.samples.hmzg <- renderPlot({
  df <- by.locus()
  if(is.null(df)) NULL else {
    df <- data.frame(x = 1:nrow(df), y = sort(df$pct.hmzg))
    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point() + geom_line() +
      xlab("Loci") + ylab("% Samples Homozygous") + 
      ylim(range(c(df[, "y"], input$sl.samples.hmzg))) +
      geom_hline(yintercept = input$sl.samples.hmzg, color = "red")
    print(p)
  }
})

observeEvent(input$btn.remove.samples.hmzg, {
  isolate({
    df <- by.locus()
    i <- which(df$pct.hmzg >= input$sl.samples.hmzg)
    if(length(i) > 0) {
      loc <- as.character(df$locus[i])
      all.loci <- locNames(vals$gtypes)
      to.keep <- setdiff(all.loci, loc)
      if(length(to.keep) > 0) {
        vals$gtypes <- vals$gtypes[ ,to.keep , ]
        vals$qaqc.reports$loci[loc, "step.removed"] <<- vals$qaqc.step
        vals$qaqc.reports$loci[loc, "threshold"] <<- input$sl.samples.hmzg
      }
    }
  })
})