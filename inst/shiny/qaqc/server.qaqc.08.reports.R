ui.reports <- function() {
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput(
        "reportGroup", label = h4("Choose reports to save"),
        choices = list(
          "By-sample summary" = "by.sample",
          "By-locus summary" = "by.locus",
          "Duplicates" = "dups",
          "Hardy-Weinberg Equilibrium" = "hwe",
          "Linkage Disequilibrium" = "ld"
        ), 
        selected = c("by.sample", "by.locus", "dups", "hwe", "ld")
      ),
      actionButton("saveReports", label = "Save selected reports"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Samples", dataTableOutput("sampleReport")),
        tabPanel("Loci", dataTableOutput("locusReport")),
        tabPanel("Duplicates", dataTableOutput("dupReport")),
        tabPanel("Hardy-Weinberg Equilibrium", uiOutput("hweReport")),
        tabPanel("Linkage Disequilibrium", uiOutput("ldReport"))
      )
    )
  )
}


output$sampleReport <- renderDataTable({
  df <- qaqc.reports$samples
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, scrollX = TRUE)
  )
})

output$locusReport <- renderDataTable({
  df <- qaqc.reports$loci
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, scrollX = TRUE)
  )
})

output$dupReport <- renderDataTable({
  df <- qaqc.reports$dups
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, scrollX = TRUE)
  )
})

output$hweReport <- renderUI({
  tabs <- lapply(strataNames(vals$gtypes), function(x) {
    tabPanel(
      title = x,
      verticalLayout(
        titlePanel(h4("HWE p-values")),
        dataTableOutput(strata.id["hwe.report", x]),
        hr(),
        titlePanel(h4("HWE Jackknife results")),
        dataTableOutput(strata.id["hwe.jack.report", x])
      )
    )
  })
  
  do.call(tabsetPanel, tabs)
})

output$ldReport <-renderUI({
  tabs <- lapply(strataNames(vals$gtypes), function(x) {
    tabPanel(
      title = x,
      dataTableOutput(strata.id["ld.report", x])
    )
  })
  
  do.call(tabsetPanel, tabs)
})