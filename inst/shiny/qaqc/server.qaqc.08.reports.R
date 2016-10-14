ui.reports <- function() {
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput(
        "reportsToSave", label = h4("Choose reports to save"),
        choices = c(
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
  df <- vals$qaqc.reports$samples
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    selection = "none",
    options = list(paging = nrow(df) > 10, scrollX = TRUE)
  )
})

output$locusReport <- renderDataTable({
  df <- vals$qaqc.reports$loci
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    selection = "none",
    options = list(paging = nrow(df) > 10, scrollX = TRUE)
  )
})

output$dupReport <- renderDataTable({
  df <- vals$qaqc.reports$dups
  if(!is.null(df)) df <- round(df, 4)
  DT::datatable(
    df, rownames = FALSE,
    selection = "none",
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

output$ldReport <- renderUI({
  tabs <- lapply(strataNames(vals$gtypes), function(x) {
    tabPanel(
      title = x,
      dataTableOutput(strata.id["ld.report", x])
    )
  })
  
  do.call(tabsetPanel, tabs)
})

observeEvent(input$saveReports, {
  isolate({
    label <- make.names(description(vals$gtypes))
    output.dir <- paste0(label, "_QAQCreports")
    output.dir <- file.path(vals$wd, output.dir)
    if(!dir.exists(output.dir)) dir.create(output.dir)
    
    if("by.sample" %in% input$reportsToSave) {
      fname <- file.path(output.dir, paste0(label, "_by.sample.report.csv"))
      write.csv(vals$qaqc.reports$samples, file = fname, row.names = FALSE)
    }
    
    if("by.locus" %in% input$reportsToSave) {
      fname <- file.path(output.dir, paste0(label, "_by.locus.report.csv"))
      write.csv(vals$qaqc.reports$loci, file = fname, row.names = FALSE)
    }
    
    if("dups" %in% input$reportsToSave) {
      fname <- file.path(output.dir, paste0(label, "_duplicates.report.csv"))
      write.csv(vals$qaqc.reports$dups, file = fname, row.names = FALSE)
    }
    
    if("hwe" %in% input$reportsToSave & !is.null(vals$qaqc.reports$hwe.jack)) {
      for(x in names(vals$qaqc.reports$hwe.jack)) {
        fname <- file.path(output.dir, paste0(label, "_HWE.pvalues_", x, ".csv"))
        write.csv(vals$qaqc.reports$hwe.jack[[x]]$p.mat, file = fname, row.names = FALSE)
        if(!is.null(vals$qaqc.reports$hwe.jack[[x]]$inf)) {
          fname <- file.path(output.dir, paste0(label, "_HWE.jack.infl_", x, ".csv"))
          write.csv(vals$qaqc.reports$hwe.jack[[x]]$inf, file = fname, row.names = FALSE)
        }
      }
    }
    
    if("ld" %in% input$reportsToSave & !is.null(vals$qaqc.reports$ld)) {
      for(x in names(vals$qaqc.reports$ld)) {
        fname <- file.path(output.dir, paste0(label, "_lnkg.diseq_", x, ".csv"))
        write.csv(vals$qaqc.reports$ld[[x]], file = fname, row.names = FALSE)
      }
    }
    
    showNotification("Reports written", duration = 2, type = "message")
  })
})