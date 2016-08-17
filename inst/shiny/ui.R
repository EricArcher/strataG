shinyUI(fluidPage(tabsetPanel(
  source("ui.gtypes.R", local = TRUE),
  source("ui.qaqc.R", local = TRUE),
  source("ui.popstruct.R", local = TRUE)
)))