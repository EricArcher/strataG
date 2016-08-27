ui.loci.hmzg <- function() {
  fluidPage(
    splitLayout(
      cellWidths = c("50%", "50%"),
      wellPanel(
        titlePanel("2) Percent of homozygous genotypes"),
        actionButton("btn.run.loci.hmzg", label = "Refresh"),
        sliderInput(
          "sl.loci.hmzg", label = NULL,
          min = 0, max = 1, value = 0.8
        ),
        textOutput("txt.loci.hmzg"),
        plotOutput("plot.loci.hmzg")
      ),
      verticalLayout(
        actionButton("btn.remove.loci.hmzg", label = "Remove samples"),
        dataTableOutput("dt.loci.hmzg")
      )
    )
  )
}

updateSliderInput(session, "sl.loci.hmzg", step = 1 / nLoc(current.g))

output$txt.loci.hmzg <- renderPrint({
  cat("Number of loci:", input$sl.loci.hmzg * nLoc(current.g))
})

output$dt.loci.hmzg <- renderDataTable({
  df <- by.sample()
  if(!is.null(df)) {
    num.hmzg <- with(df, pct.loci.homozygous * (nLoc(current.g) - num.loci.missing.genotypes))
    df <- df[, c("id", "strata", "pct.loci.homozygous")]
    df <- cbind(df, '# Homozygous' = num.hmzg)
    colnames(df) <- c("ID", "Strata", "% Homozygous", "# Homozygous")
    df <- df[order(df[, 3], decreasing = TRUE), ]
    df <- df[df[, 3] >= input$sl.loci.hmzg, ]
  }
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE)
  )
})

output$plot.loci.hmzg <- renderPlot({
  df <- by.sample()
  if(is.null(df)) NULL else {
    df <- data.frame(x = 1:nrow(df), y = sort(df$pct.loci.homozygous))
    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point() + geom_line() +
      xlab("Samples") + ylab("% Genotypes Homozygous") + 
      ylim(range(c(df[, "y"], input$sl.loci.hmzg))) +
      geom_hline(yintercept = input$sl.loci.hmzg, color = "red")
    print(p)
  }
})

observeEvent(input$btn.remove.loci.hmzg, {
  df <- by.sample()
  i <- which(df$pct.loci.homozygous >= input$sl.loci.hmzg)
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