qaqcDiploid <- function() {
  tabsetPanel(
    tabPanel(
      "Workflow",
      uiOutput("qaqcDiploidStep")
    ),
    tabPanel(
      "Samples", 
      sidebarLayout(
        sidebarPanel(
          radioButtons(
            "sampleSelect", "Select samples to view",
            choices = c("All" = "all", "Not removed" = "not.removed", "Removed" = "removed"),
            selected = "all"
          ),
          actionButton("btn.remove.samples", label = "Remove selected samples"),
          width = 3
        ), 
        mainPanel(dataTableOutput("dt.by.sample"))
      )
    ),
    tabPanel(
      "Loci",     
      sidebarLayout(
        sidebarPanel(
          radioButtons(
            "locusSelect", "Select loci to view",
            choices = c("All" = "all", "Not removed" = "not.removed", "Removed" = "removed"),
            selected = "all"
          ),
          actionButton("btn.remove.loci", label = "Remove selected loci"),
          width = 3
        ), 
        mainPanel(dataTableOutput("dt.by.locus"))
      )
    )
  )
}