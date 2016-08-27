ui.gtypes <- function() {
  tabPanel(
    title = "Load data",
    sidebarLayout(
      
      sidebarPanel(
        radioButtons(
          "data.type",
          label = h4("Data source"),
          choices = list("Comma-separated file (.csv)" = 1, "Existing gtypes object (.rdata)" = 2),
          selected = 1
        ),
        uiOutput("gtypesLoadData")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("Genetic data", dataTableOutput("file.gen.data")),
          tabPanel("Stratification schemes", dataTableOutput("file.schemes")),
          tabPanel("FASTA", verbatimTextOutput("file.fasta"))
        ),
        uiOutput("loadedGtype")
      )
      
    )
  )
}