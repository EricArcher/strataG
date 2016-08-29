library(shiny)
library(shinyFiles)
library(DT)
library(ggplot2)
library(reshape2)
library(strataG)
library(parallel)

source("ui.load.R", local = TRUE)
source("ui.gtypes.R", local = TRUE)
source("ui.qaqc.R", local = TRUE)
source("ui.popstruct.R", local = TRUE)

shinyApp(
  ui = shinyUI(fluidPage(
    titlePanel("strataG GUI"),
    tabsetPanel(
      id = "main.tab",
      ui.load(), ui.gtypes(), ui.qaqc(), ui.popstruct()
    ) 
  )),
  server = shinyServer(function(input, output, session) {
    user.data <- reactiveValues(current.g = NULL)
    qaqc.flow <- reactiveValues(step = 1)
    volumes <- getVolumes()
    source("server.load.R", local = TRUE)
    source("server.gtypes.R", local = TRUE)
    source("server.qaqc.R", local = TRUE)
    source("server.popstruct.R", local = TRUE)
  })
)