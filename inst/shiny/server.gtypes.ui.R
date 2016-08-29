output$schemeMenu <- renderUI({
  items <- if(is.null(schemes.df())) NULL else {
    x <- c("Do not stratify", colnames(schemes.df()))
    names(x) <- x
    x
  }
  fluidRow(
    column(3, selectInput("selected.scheme", NULL, items)),
    column(1, actionButton("stratify", label = "Stratify"))
  )
})
