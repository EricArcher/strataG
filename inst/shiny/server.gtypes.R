source("server.gtypes.ui.R", local = TRUE)

# view current gtypes - updates when user.data$current.g is updated
output$gtypeSmry <- renderPrint({
  options(digits = 3)
  current.g <<- user.data$current.g
  user.data$current.g
})

# stratify current gtypes - responds to 'stratify' actionButton
observe({
  input$stratify
  isolate({
    if(!is.null(user.data$current.g)) {
      x <- input$selected.scheme
      if(x == "Do not stratify") {
        strata(user.data$current.g) <- "Default"
      } else if(x != "") {
        user.data$current.g <- stratify(user.data$current.g, x)
      }
      current.g <- user.data$current.g
    }
  })
})

# save gtypes object to working directory
observe({
  input$save.gtypes
  isolate({
    save.dir <- parseDirPath(volumes, input$wd)
    if(!is.null(user.data$current.g) & !is.null(save.dir)) {
      if(save.dir != "") {
        g <- user.data$current.g
        fname <- paste0(make.names(description(g)), " gtypes.rdata")
        save(g, file = file.path(save.dir, fname))
      }
    }
  })
})