source("server.load.ui.sidepanel.csv.R", local = TRUE)
source("server.load.ui.sidepanel.rdata.R", local = TRUE)
source("server.load.display.gtypes.R", local = TRUE)

output$loadSidepanel <- renderUI({
  list(uiLoadCsv, uiLoadRdata)[[as.numeric(input$loadDataSource)]]()
})

output$loadedData <- renderUI({
  if(!is.null(vals$gtypes)) {
    verticalLayout(
      titlePanel(h4("Loaded gtypes")),
      selectInput("selectedScheme", NULL, NULL),
      verbatimTextOutput("gtypeSmry"),
      hr(),
      fluidRow(
        column(4, textInput("gtypesName", "Name for gtypes object in .rdata file")),
        column(3, actionButton("save.gtypes", label = "Save gtypes to .rdata file"))
      )
    )
  } else if(length(input$ploidy) == 0 | is.null(input$ploidy)) {
    NULL
  } else if(input$ploidy == "1") {
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