library(shiny)
library(ggplot2)
library(reshape2)
source("ui.gtypes.R", local = TRUE)
source("ui.qaqc.R", local = TRUE)
source("ui.popstruct.R", local = TRUE)

shinyUI(
  fluidPage(
    titlePanel("strataG GUI"),
    tabsetPanel(ui.gtypes(), ui.qaqc(), ui.popstruct())
  )
)