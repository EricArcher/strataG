ui.setup <- function() {
  tabPanel(
    title = "Setup",
    shinyDirButton(
      id = "workingDirBtn",
      label = "Choose working directory",
      title = "Choose working directory"
    ),
    verbatimTextOutput("workingDir"),
    helpText("This is the directory that contains your data and where output will be saved."),
    hr(),
    radioButtons(
      inputId = "loadDataSource",
      label = "Choose the type of data you'll be loading",
      choices = c(
        "Comma-delimited text files (.csv)" = 1,
        "Existing gtypes object (.rdata)" = 2
      )
    )
  )
}