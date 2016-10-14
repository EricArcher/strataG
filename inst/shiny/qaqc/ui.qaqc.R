source("ui.qaqc.haploid.R", local = TRUE)
source("ui.qaqc.diploid.R", local = TRUE)

ui.qaqc <- function() tabPanel("QA/QC", uiOutput("uiQAQC"))
