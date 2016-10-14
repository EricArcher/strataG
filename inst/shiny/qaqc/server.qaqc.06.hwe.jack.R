ui.hwe.jack <- function() {
  tabs <- lapply(strataNames(vals$gtypes), function(x) {
    tabPanel(
      title = x,
      verticalLayout(
        titlePanel(h4("HWE p-values")),
        dataTableOutput(strata.id["dt.hwe", x]),
        hr(),
        titlePanel(h4("HWE Jackknife results")),
        dataTableOutput(strata.id["dt.hwe.jack", x]),
        hr(),
        plotOutput(strata.id["plot.hwe.jack", x])
      )
    )
  })
  
  sidebarLayout(
    sidebarPanel(
      NULL,
      radioButtons(
        "hwe.source", h4("Choose source for HWE calculations"),
        list("Exact test (pegas)" = 1, "Heterozygote deficiency (Genepop)" = 2), 
        selected = 1
      ),
      sliderInput("sl.hwe.jack", h4("HWE critical alpha"), 0, 0.2, 0.05),
      actionButton("calc.hwe", "Calculate HWE and Jackknife"),
      width = 3
    ),
    mainPanel(
      do.call(tabsetPanel, tabs),
      uiOutput("hwe.remove")
    )
  )
}

observeEvent(input$calc.hwe, {
  jack.results <- hw.jack()
  
  loadTabs <- function(x) {
    # format HWE p value matrix
    df <- jack.results[[x]]$obs
    if(!is.null(df)) df <- round(hwe.p.mat(df), 4)
    output[[strata.id["dt.hwe", x]]] <- renderDataTable({ 
      DT::datatable(
        df, rownames = FALSE,
        options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
      )
    })
    
    # format jackknife table
    infl <- jackInfluential(jack.results[[x]], alpha = input$sl.hwe.jack)
    infl.df <- infl$influential
    if(!is.null(infl.df)) infl.df <- round(infl.df, 4)
    output[[strata.id["dt.hwe.jack", x]]] <- renderDataTable({
      DT::datatable(
        infl.df, rownames = FALSE,
        options = list(paging = nrow(infl.df) > 10, searching = FALSE, scrollX = TRUE)
      )
    })
    
    # plot jackknife odds-ratio distribution
    if(!is.null(jack.results[[x]])) {
      output[[strata.id["plot.hwe.jack", x]]] <- renderPlot(plot(infl))
    }
  }
  
  for(x in colnames(strata.id)) loadTabs(x)
})

output$hwe.remove <- renderUI({
  if(first.run) return(NULL)
  df <- do.call(rbind, lapply(hw.jack(), function(x) {
    jackInfluential(x, alpha = input$sl.hwe.jack)$influential
  }))
  if(is.null(df)) return(NULL)
  ids <- unique(df$excluded)
  loci <- unique(df$locus)
  ids <- ids[order(nchar(ids), ids)]
  loci <- loci[order(nchar(loci), loci)]
  verticalLayout(
    hr(),
    splitLayout(
      verticalLayout(
        selectizeInput(
          "hw.samples.select", "Select samples to remove",
          choices = ids, multiple = TRUE
        ),
        actionButton("btn.remove.hw.samples", "Remove selected samples")
      ),
      verticalLayout(
        selectizeInput(
          "hw.loci.select", "Select loci to remove",
          choices = loci, multiple = TRUE
        ),
        actionButton("btn.remove.hw.loci", "Remove selected loci")
      )
    )
  )
})

observeEvent(input$btn.remove.hw.samples, {
  isolate({
    if(!is.null(input$hw.samples.select)) {
      id <- input$hw.samples.select
      all.inds <- indNames(vals$gtypes)
      to.keep <- setdiff(all.inds, id)
      if(length(to.keep) > 0) {
        vals$gtypes <- vals$gtypes[to.keep, , ]
        vals$qaqc.reports$samples[id, "step.removed"] <- vals$qaqc.step
        vals$qaqc.reports$samples[id, "threshold"] <- input$sl.hwe.jack
      }
    }
  })
})

observeEvent(input$btn.remove.hw.loci, {
  isolate({
    if(!is.null(input$hw.loci.select)) {
      loc <- input$hw.loci.select
      all.loci <- locNames(vals$gtypes)
      to.keep <- setdiff(all.loci, loc)
      if(length(to.keep) > 0) {
        vals$gtypes <- vals$gtypes[ ,to.keep , ]
        vals$qaqc.reports$loci[loc, "step.removed"] <- vals$qaqc.step
        vals$qaqc.reports$loci[loc, "threshold"] <- input$sl.hwe.jack
      }
    }
  })
})