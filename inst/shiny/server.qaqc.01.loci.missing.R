ui.loci.missing <- function() {
  fluidPage(
    splitLayout(
      cellWidths = c("50%", "50%"),
      wellPanel(
        titlePanel("1) Percent of missing genotypes"),
        actionButton("btn.run.loci.missing", label = "Refresh"),
        sliderInput(
          "sl.loci.missing", label = NULL,
          min = 0, max = 1, value = 0.8
        ),
        textOutput("txt.loci.missing"),
        plotOutput("plot.loci.missing")
      ),
      verticalLayout(
        actionButton("btn.remove.loci.missing", label = "Remove samples"),
        dataTableOutput("dt.loci.missing")
      )
    )
  )
}

updateSliderInput(session, "sl.loci.missing", step = 1 / nLoc(current.g))

output$txt.loci.missing <- renderPrint({
  cat("Number of loci:", input$sl.loci.missing * nLoc(current.g))
})

output$dt.loci.missing <- renderDataTable({
  df <- by.sample()
  if(!is.null(df)) {
    df <- df[, c("id", "strata", "pct.loci.missing.genotypes", "num.loci.missing.genotypes")]
    colnames(df) <- c("ID", "Strata", "% Missing", "# Missing")
    df <- df[order(df[, 3], decreasing = TRUE), ]
    df <- df[df[, 3] >= input$sl.loci.missing, ]
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE)
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
      current.g <<- current.g[to.keep, , ]
      reloaded$count <- reloaded$count + 1
    }
  }
})