ui.setup <- function() {
  tabPanel(
    title = "Setup",
    fluidRow(
      column(
        width = 2,
        shinyDirButton(
          id = "workingDirBtn",
          label = "Choose working directory",
          title = "Choose working directory"
        )
      ),
      column(
        width = 6,
        verbatimTextOutput("workingDir")
      )
    ),
    hr(),
    radioButtons(
      inputId = "loadDataSource",
      label = "Choose data source",
      choices = c(
        "Comma-delimited text files (.csv)" = 1,
        "Existing gtypes object (.rdata)" = 2
      )
    )
  )
}