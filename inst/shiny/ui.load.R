ui.load.sidebar <- function() {
  sidebarPanel(
    shinyDirButton("wd", "Select working directory", "Select a directory"),
    verbatimTextOutput("wdPath"),
    hr(),
    radioButtons(
      "data.type", label = h4("Data source"),
      choices = list(
        "Comma-separated file (.csv)" = 1, 
        "Existing gtypes object (.rdata)" = 2
      ),
      selected = 1
    ),
    uiOutput("loadSource")
  )
}
  
ui.load <- function() {
  tabPanel(
    title = "Load data",
    sidebarLayout(ui.load.sidebar(), mainPanel(uiOutput("loadedData")))
  )
}