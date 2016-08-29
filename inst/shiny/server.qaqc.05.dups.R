ui.dups <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "sl.dups", label = h4("5) Percent of loci shared"),
        min = 0, max = 1, value = 0.8
      ),
      textOutput("txt.dups"),
      plotOutput("plot.dups")
    ),
    mainPanel(dataTableOutput("dt.dups"))
  )
}

output$txt.dups <- renderPrint({
  if(is.null(user.data$current.g)) return()
  cat("Number of loci:", input$sl.dups * nLoc(user.data$current.g))
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