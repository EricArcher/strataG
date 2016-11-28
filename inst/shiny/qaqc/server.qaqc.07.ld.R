ui.ld <- function() {
  tabs <- lapply(strataNames(vals$gtypes), function(x) {
    tabPanel(
      title = x,
      splitLayout(
        plotOutput(strata.id["plot.ld", x]),
        dataTableOutput(strata.id["dt.ld", x])
      )
    )
  })
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("sl.ld", h4("Linkage disequilibrium p-value"), 0, 0.2, 0.05),
      actionButton("calc.ld", "Calculate Linkage Disequilibrium"),
      width = 3
    ),
    mainPanel(
      do.call(tabsetPanel, tabs),
      uiOutput("ld.remove")
    )
  )
}

observeEvent(input$calc.ld, {
  ld.results <- ld()
  
  loadTabs <- function(x) {
    df <- ld.results[[x]]
    # format and load table
    if(!is.null(df)) {
      df <- round(df, 4)
      df <- df[df$p.value <= input$sl.ld, ]
      df <- df[order(df$p.value, df$Locus.1, df$Locus.2), 1:4]
    }
    output[[strata.id["dt.ld", x]]] <- renderDataTable({
      DT::datatable(
        df, rownames = FALSE,
        options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
      )
    })
    
    # plot p-value distribution
    output[[strata.id["plot.ld", x]]] <- renderPlot({
      if(is.null(df)) NULL else if(nrow(df) == 0) NULL else {
        df <- data.frame(x = 1:nrow(df), y = sort(df$p.value, na.last = TRUE))
        p <- ggplot(df, aes_string(x = "x", y = "y")) + 
          geom_point() + geom_line() +
          xlab("Pairs") + ylab("p-value") + 
          ylim(range(c(df[, "y"], input$sl.ld))) +
          geom_hline(yintercept = input$sl.ld, color = "red")
        print(p)
      }
    })
  }
  
  for(x in colnames(strata.id)) loadTabs(x)
})

output$ld.remove <- renderUI({
  if(first.run) return(NULL)
  df <- do.call(rbind, lapply(ld(), function(x) {
    if(is.null(x)) return(NULL)
    x[x$p.value <= input$sl.ld, ]
  }))
  if(is.null(df)) return(NULL)
  loci <- unique(unlist(df[, 1:2]))
  verticalLayout(
    hr(),
    selectizeInput(
      "ld.loci.select", "Select loci to remove",
      choices = loci, multiple = TRUE
    ),
    actionButton("btn.remove.ld.loci", "Remove selected loci")
  )
})

observeEvent(input$btn.remove.ld.loci, {
  isolate({
    if(!is.null(input$ld.loci.select)) {
      loc <- input$ld.loci.select
      all.loci <- locNames(vals$gtypes)
      to.keep <- setdiff(all.loci, loc)
      if(length(to.keep) > 0) {
        vals$gtypes <- vals$gtypes[ ,to.keep , ]
        vals$qaqc.reports$loci[loc, "step.removed"] <- vals$qaqc.step
        vals$qaqc.reports$loci[loc, "threshold"] <- input$sl.ld
      }
    }
  })
})