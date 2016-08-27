ui.hwe.jack <- function() {
  fluidPage(
    wellPanel(
      titlePanel("6) Hardy-Weinberg p-values"),
      actionButton("btn.run.hwe", label = "Refresh"),
      dataTableOutput("dt.hwe")
    ),
    splitLayout(
      cellWidths = c("50%", "50%"),
      wellPanel(
        titlePanel("Hardy-Weinberg jackknife critical alpha"),
        actionButton("btn.run.hwe.jack", label = "Refresh"),
        wellPanel(
          sliderInput(
            "sl.hwe.jack", label = NULL,
            min = 0, max = 0.2, value = 0.05
          ),
          plotOutput("plot.hwe.jack")
        )
      ),
      dataTableOutput("dt.hwe.jack")
    )
  )
}

output$dt.hwe <- renderDataTable({
  df <- hwe()
  if(!is.null(df)) {
    df <- t(sapply(p.adjust.methods, function(m) p.adjust(df, method = m)))
    df <- data.frame(method = rownames(df), df)
    df <- round(df, 4)
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE)
  )
})
            
output$dt.hwe.jack <- renderDataTable({
  infl <- jackInfluential(hw.jack(), alpha = input$sl.hwe.jack)
  df <- infl$influential
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 25, searching = FALSE)
  )
})

output$plot.hwe.jack <- renderPlot({
  infl <- jackInfluential(hw.jack(), alpha = input$sl.hwe.jack)
  plot(infl)
})
