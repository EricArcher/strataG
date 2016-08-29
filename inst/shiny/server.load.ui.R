source("ui.load.sidepanel.R", local = TRUE)

output$loadSource <- renderUI({
  list(ui.load.csv, ui.load.rdata)[[as.numeric(input$data.type)]]()
})

output$loadedData <- renderUI({
  if(length(input$ploidy) == 0 | is.null(input$ploidy)) return()
  if(input$ploidy == "1") {
    tabsetPanel(
      tabPanel("Genetic data", dataTableOutput("file.gen.data")),
      tabPanel("Stratification schemes", dataTableOutput("file.schemes")),
      tabPanel("Sequences", verbatimTextOutput("file.fasta"))
    )
  } else {
    tabsetPanel(
      tabPanel("Genetic data", dataTableOutput("file.gen.data")),
      tabPanel("Stratification schemes", dataTableOutput("file.schemes"))
    )
  }
})

