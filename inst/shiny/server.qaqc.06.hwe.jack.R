ui.hwe.jack <- function() {
  sidebarLayout(
    sidebarPanel(
      titlePanel(h4("6) Hardy-Weinberg p-values")),
      dataTableOutput("dt.hwe"),
      hr(),
      sliderInput(
        "sl.hwe.jack", label = h4("Hardy-Weinberg jackknife critical alpha"),
        min = 0, max = 0.2, value = 0.05
      ),
      plotOutput("plot.hwe.jack")
    ),
    mainPanel(dataTableOutput("dt.hwe.jack"))
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
    options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
  )
})
            
output$dt.hwe.jack <- renderDataTable({
  jack.result <- hw.jack()
  if(is.null(jack.result)) return()
  infl <- jackInfluential(jack.result, alpha = input$sl.hwe.jack)
  df <- infl$influential
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 25, searching = FALSE, scrollX = TRUE)
  )
})

output$plot.hwe.jack <- renderPlot({
  jack.result <- hw.jack()
  if(is.null(jack.result)) return()
  infl <- jackInfluential(jack.result, alpha = input$sl.hwe.jack)
  plot(infl)
})
