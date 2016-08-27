ui.samples.hmzg <- function() {
  fluidPage(
    splitLayout(
      cellWidths = c("50%", "50%"),
      wellPanel(
        titlePanel("4) Percent of homozygous loci"),
        actionButton("btn.run.samples.hmzg", label = "Refresh"),
        sliderInput(
          "sl.samples.hmzg", label = NULL,
          min = 0, max = 1, value = 0.8
        ),
        textOutput("txt.samples.hmzg"),
        plotOutput("plot.samples.hmzg")
      ),
      verticalLayout(
        actionButton("btn.remove.samples.hmzg", label = "Remove loci"),
        dataTableOutput("dt.samples.hmzg")
      )
    )
  )
}

updateSliderInput(session, "sl.samples.hmzg", step = 1 / nInd(current.g))

output$txt.samples.hmzg <- renderPrint({
  cat("Number of samples:", input$sl.samples.hmzg * nInd(current.g))
})

output$dt.samples.hmzg <- renderDataTable({
  df <- by.locus()
  if(!is.null(df)) {
    df <- df[, c("locus", "pct.hmzg", "num.hmzg")]
    colnames(df) <- c("Locus", "% Homozygous", "# Homozygous")
    df <- df[order(df[, 2], decreasing = TRUE), ]
    df <- df[df[, 2] >= input$sl.samples.hmzg, ]
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE)
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
    all.loci <- locNames(current.g)
    to.keep <- setdiff(all.loci, loc)
    if(length(to.keep) > 0) {
      current.g <<- current.g[ ,to.keep , ]
      reloaded$count <- reloaded$count + 1
    }
  }
})