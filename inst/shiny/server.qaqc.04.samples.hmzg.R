ui.samples.hmzg <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "sl.samples.hmzg", label = h4("4) Percent of homozygous loci"),
        min = 0, max = 1, value = 0.8
      ),
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
  if(is.null(user.data$current.g)) return()
  cat("Number of samples:", input$sl.samples.hmzg * nInd(user.data$current.g))
})

output$dt.samples.hmzg <- renderDataTable({
  df <- by.locus()
  if(!is.null(df)) {
    df <- df[, c("locus", "pct.hmzg", "num.hmzg")]
    colnames(df) <- c("Locus", "% Homozygous", "# Homozygous")
    df <- df[order(df[, 2], decreasing = TRUE), ]
    df <- df[df[, 2] >= input$sl.samples.hmzg, ]
    df <- round(df, 4)
  }
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
  df <- by.locus()
  i <- which(df$pct.hmzg >= input$sl.samples.hmzg)
  if(length(i) > 0) {
    loc <- df$locus[i]
    all.loci <- locNames(user.data$current.g)
    to.keep <- setdiff(all.loci, loc)
    if(length(to.keep) > 0) {
      user.data$current.g <- user.data$current.g[ ,to.keep , ]
    }
  }
})