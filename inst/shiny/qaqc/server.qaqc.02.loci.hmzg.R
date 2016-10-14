ui.loci.hmzg <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput("sl.loci.hmzg", h4("Percent of homozygous loci"), 0, 1, 0.8),
      textOutput("txt.loci.hmzg"),
      plotOutput("plot.loci.hmzg")
    ),
    mainPanel(
      dataTableOutput("dt.loci.hmzg"),
      actionButton("btn.remove.loci.hmzg", label = "Remove samples")
    )
  )
}

output$txt.loci.hmzg <- renderPrint({
  if(is.null(vals$gtypes)) return()
  num.hmzg <- floor(input$sl.loci.hmzg * nLoc(vals$gtypes))
  cat("Number of loci:", num.hmzg)
})

df.loci.hmzg <- reactive({
  df <- by.sample()
  if(!is.null(df)) {
    num.genotyped <- nLoc(vals$gtypes) - df$num.loci.missing.genotypes
    num.hmzg <- df$pct.loci.homozygous * num.genotyped
    df <- df[, c("id", "strata", "pct.loci.homozygous")]
    df <- cbind(df, '# Homozygous' = num.hmzg)
    colnames(df) <- c("ID", "Strata", "% Homozygous", "# Homozygous")
    df <- df[order(df[, 3], decreasing = TRUE), ]
    df <- df[df[, 3] >= input$sl.loci.hmzg, ]
    df <- round(df, 4)
  }
  df
})

output$dt.loci.hmzg <- renderDataTable({
  df <- df.loci.hmzg()
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
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
  isolate({
    df <- by.sample()
    i <- which(df$pct.loci.homozygous >= input$sl.loci.hmzg)
    if(length(i) > 0) {
      id <- as.character(df$id[i])
      all.inds <- indNames(vals$gtypes)
      to.keep <- setdiff(all.inds, id)
      if(length(to.keep) > 0) {
        vals$gtypes <- vals$gtypes[to.keep, , ]
        vals$qaqc.reports$samples[id, "step.removed"] <<- vals$qaqc.step
        vals$qaqc.reports$samples[id, "threshold"] <<- input$sl.loci.hmzg
      }
    }
  })
})