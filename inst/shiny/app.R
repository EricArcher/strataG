library(shiny)
library(shinyFiles)
library(DT)
library(ggplot2)
library(reshape2)
library(strataG)
library(parallel)

source(file.path("setup", "ui.setup.R"), local = TRUE, chdir = TRUE)
source(file.path("load.data", "ui.load.R"), local = TRUE, chdir = TRUE)
source(file.path("qaqc", "ui.qaqc.R"), local = TRUE, chdir = TRUE)
source(file.path("popstruct", "ui.popstruct.R"), local = TRUE, chdir = TRUE)

shinyApp(
  ui = fluidPage(
    verticalLayout(
      titlePanel("strataG GUI"),
      hr(),
      tabsetPanel(
        id = "main.tab", ui.setup(), ui.load(), ui.qaqc(), ui.popstruct()
      )
    )
  ),
  
  server = function(input, output, session) {
    vals <- reactiveValues(
      wd = NULL,
      scratch.env = new.env(),
      gtypes = NULL,
      qaqc.step = 1,
      qaqc.reports = NULL
    )
    id <- NULL
    
    source(file.path("setup", "server.setup.R"), local = TRUE, chdir = TRUE)
    source(file.path("load.data", "server.load.R"), local = TRUE, chdir = TRUE)
    source(file.path("qaqc", "server.qaqc.R"), local = TRUE, chdir = TRUE)
    source(file.path("qaqc", "server.qaqc.haploid.R"), local = TRUE, chdir = TRUE)
    source(file.path("popstruct", "server.popstruct.R"), local = TRUE, chdir = TRUE)
  }
)