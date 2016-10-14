ui.popstruct <- function() {
  tabPanel(
    title = "Population structure",
    sidebarLayout(
      sidebarPanel(
        titlePanel(h4("Type")),
        fluidRow(
          column(3, checkboxInput("ovl", "Overall", value = TRUE)),
          column(3, checkboxInput("pws", "Pairwise", value = TRUE))
        ),
        hr(),
        titlePanel(h4("Metrics")),
        uiOutput("metrics"),
        numericInput(
          "nrep", label = h4("Number of permutations"), 
          min = 0, value = 1000
        ),
        actionButton("run.popstruct", label = "Run")
      ),
      mainPanel(
        actionButton("savePopStructResults", label = "Save results"),
        hr(),
        titlePanel(h4("Overall")), 
        dataTableOutput("ovlResults"),
        hr(),
        titlePanel(h4("Pairwise")),
        dataTableOutput("pwsResults")
      )
    )
  )
}