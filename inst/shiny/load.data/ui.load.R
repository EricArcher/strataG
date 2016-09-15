ui.load <- function() {
  tabPanel(
    title = "Load data",
    sidebarLayout(
      sidebarPanel(uiOutput("loadSidepanel")), 
      mainPanel(uiOutput("loadedData"))
    )
  )
}