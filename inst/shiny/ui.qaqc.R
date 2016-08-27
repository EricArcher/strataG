source("ui.qaqc.haploid.R", local = TRUE)
source("ui.qaqc.diploid.R", local = TRUE)

ui.qaqc <- function() {
  tabPanel(
    "QA/QC",
    qaqcDiploid(),
    titlePanel("Current gtypes object"),
    actionButton("btn.save.gtypes", "Save gtypes object"),
    actionButton("btn.write.gtypes", "Write .csv file"),
    mainPanel(verbatimTextOutput("gtypes.smry"))
  )
}