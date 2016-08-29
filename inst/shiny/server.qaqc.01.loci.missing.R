ui.loci.missing <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "sl.loci.missing", label = h4("1) Percent of missing genotypes"),
        min = 0, max = 1, value = 0.8
      ),
      textOutput("txt.loci.missing"),
      plotOutput("plot.loci.missing")
    ),
    mainPanel(
      dataTableOutput("dt.loci.missing"),
      actionButton("btn.remove.loci.missing", label = "Remove samples")
    )
  )
}

output$txt.loci.missing <- renderPrint({
  if(is.null(user.data$current.g)) return()
  cat("Number of loci:", input$sl.loci.missing * nLoc(user.data$current.g))
})

output$dt.loci.missing <- renderDataTable({
  df <- by.sample()
  if(!is.null(df)) {
    df <- df[, c("id", "strata", "pct.loci.missing.genotypes", "num.loci.missing.genotypes")]
    colnames(df) <- c("ID", "Strata", "% Missing", "# Missing")
    df <- df[order(df[, 3], decreasing = TRUE), ]
    df <- df[df[, 3] >= input$sl.loci.missing, ]
    df <- round(df, 4)
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
  )
})

output$plot.loci.missing <- renderPlot({
  df <- by.sample()
  if(is.null(df)) NULL else {
    df <- data.frame(x = 1:nrow(df), y = sort(df$pct.loci.missing.genotypes))
    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point() + geom_line() +
      xlab("Samples") + ylab("% Genotypes Missing") + 
      ylim(range(c(df[, "y"], input$sl.loci.missing))) +
      geom_hline(yintercept = input$sl.loci.missing, color = "red")
    print(p)
  }
})

observeEvent(input$btn.remove.loci.missing, {
  df <- by.sample()
  i <- which(df$pct.loci.missing.genotypes >= input$sl.loci.missing)
  if(length(i) > 0) {
    id <- df$id[i]
    all.inds <- indNames(current.g)
    to.keep <- setdiff(all.inds, id)
    if(length(to.keep) > 0) {
      user.data$current.g <- user.data$current.g[to.keep, , ]
    }
  }
})