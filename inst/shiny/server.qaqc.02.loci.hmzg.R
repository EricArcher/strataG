ui.loci.hmzg <- function() {
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "sl.loci.hmzg", label = h4("2) Percent of homozygous genotypes"),
        min = 0, max = 1, value = 0.8
      ),
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
  if(is.null(user.data$current.g)) return()
  cat("Number of loci:", input$sl.loci.hmzg * nLoc(user.data$current.g))
})

output$dt.loci.hmzg <- renderDataTable({
  df <- by.sample()
  if(!is.null(df)) {
    num.genotyped <- nLoc(user.data$current.g) - df$num.loci.missing.genotypes
    num.hmzg <- df$pct.loci.homozygous * num.genotyped
    df <- df[, c("id", "strata", "pct.loci.homozygous")]
    df <- cbind(df, '# Homozygous' = num.hmzg)
    colnames(df) <- c("ID", "Strata", "% Homozygous", "# Homozygous")
    df <- df[order(df[, 3], decreasing = TRUE), ]
    df <- df[df[, 3] >= input$sl.loci.hmzg, ]
    df <- round(df, 4)
  }
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
  df <- by.sample()
  i <- which(df$pct.loci.homozygous >= input$sl.loci.hmzg)
  if(length(i) > 0) {
    id <- df$id[i]
    all.inds <- indNames(user.data$current.g)
    to.keep <- setdiff(all.inds, id)
    if(length(to.keep) > 0) {
      user.data$current.g <- user.data$current.g[to.keep, , ]
    }
  }
})