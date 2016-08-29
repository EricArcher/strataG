ui.gtypes <- function() {
  tabPanel(
    title = "Loaded gtypes",
    actionButton("save.gtypes", label = "Save gtypes to .rdata file"),
    hr(),
    helpText("Select stratification scheme"),
    uiOutput("schemeMenu"),
    hr(),
    verbatimTextOutput("gtypeSmry")
  )
}