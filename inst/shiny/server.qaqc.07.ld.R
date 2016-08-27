ui.ld <- function() {
  fluidPage(
    splitLayout(
      cellWidths = c("50%", "50%"),
      wellPanel(
        titlePanel("7) Linkage disequilibrium p-values"),
        actionButton("btn.run.ld", label = "Refresh"),
        wellPanel(
          sliderInput(
            "sl.ld", label = NULL,
            min = 0, max = 0.2, value = 0.05
          ),
          plotOutput("plot.ld")
        )
      ),
      dataTableOutput("dt.ld")
    )
  )
}
            
output$dt.ld <- renderDataTable({
  df <- ld()
  if(!is.null(df)) {
    df <- round(df, 3)
    df <- df[df$p.value <= input$sl.ld, ]
    df <- df[order(df$p.value, df$Locus.1, df$Locus.2), 1:4]
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE)
  )
})

output$plot.ld <- renderPlot({
  df <- ld()
  if(is.null(df)) NULL else {
    df <- data.frame(x = 1:nrow(df), y = sort(df$p.value, decreasing = TRUE))
    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point() + geom_line() +
      xlab("Pairs") + ylab("p-value") + 
      ylim(range(c(df[, "y"], input$sl.ld))) +
      geom_hline(yintercept = input$sl.ld, color = "red")
    print(p)
  }
})
