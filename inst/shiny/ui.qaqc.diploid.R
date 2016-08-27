qaqcDiploid <- function() {
  tabsetPanel(
    tabPanel(
      "Workflow",
      fluidRow(
        actionButton("btn.qaqc.next", label = "Next"),
        actionButton("btn.qaqc.reset", label = "Reset")
      ),
      uiOutput("qaqcDiploidStep")
    ),
    tabPanel(
      "Samples",  
      actionButton("btn.remove.samples", label = "Remove selected"),
      dataTableOutput("dt.by.sample")
    ),
    tabPanel(
      "Loci",
      actionButton("btn.remove.loci", label = "Remove selected"),
      dataTableOutput("dt.by.locus")
    )
  )
}