ui.dups <- function() {
  fluidPage(
    splitLayout(
      cellWidths = c("50%", "50%"),
      wellPanel(
        titlePanel("5) Percent of loci shared"),
        actionButton("btn.run.dups", label = "Refresh"),
        sliderInput(
          "sl.dups", label = NULL,
          min = 0, max = 1, value = 0.8
        ),
        textOutput("txt.dups"),
        plotOutput("plot.dups")
      ),
      verticalLayout(
        dataTableOutput("dt.dups")
      )
    )
  )
}

updateSliderInput(session, "sl.dups", step = 1 / nLoc(current.g))

output$txt.dups <- renderPrint({
  cat("Number of loci:", input$sl.dups * nLoc(current.g))
})

output$dt.dups <- renderDataTable({
  df <- dups()
  if(!is.null(df)) df <- df[df$prop.loci.shared >= input$sl.dups, ]
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE)
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
