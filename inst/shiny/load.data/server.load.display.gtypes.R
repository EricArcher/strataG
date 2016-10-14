# view current gtypes - updates when vals$gtypes is updated
output$gtypeSmry <- renderPrint({
  items <- if(is.null(vals$gtypes)) {
    NULL
  } else if(is.null(schemes(vals$gtypes))) NULL else {
    x <- c("Do not stratify", colnames(schemes(vals$gtypes)))
    names(x) <- x
    x <- c("Select a stratification scheme" = "", x)
    x
  }
  
  updateSelectInput(session, "selectedScheme", NULL, items)
  updateTextInput(
    session, "gtypesName", 
    value = if(is.null(vals$gtypes)) NULL else {
      make.names(description(vals$gtypes))
    }
  )
  options(digits = 3)
  vals$gtypes
})

# stratify current gtypes - responds to 'stratify' actionButton
observeEvent(input$selectedScheme, {
  isolate({
    if(!is.null(vals$gtypes)) {
      x <- input$selectedScheme
      if(is.null(x)) return()
      if(x == "Do not stratify") {
        strata(vals$gtypes) <- "Default"
      } else if(x != "") {
        vals$gtypes <- stratify(vals$gtypes, x, drop = FALSE)
      }
      gtypesLoaded()
    }
  })
})

# responds to changes in vals$gtypes
# create ids for controls on QA/QC HWE and LD strata tabs
observe({
  if(!is.null(vals$gtypes)) {
    sn <- strataNames(vals$gtypes)
    strata.id <<- sapply(sn, function(x) {
      c(
        dt.hwe = paste0("dt.hwe.", x),
        plot.hwe.jack = paste0("plot.hwe.jack.", x),
        dt.hwe.jack = paste0("dt.hwe.jack.", x),
        plot.ld = paste0("plot.ld.", x),
        dt.ld = paste0("dt.ld.", x), 
        hwe.report = paste0("hwe.", x),
        hwe.jack.report = paste0("hwe.jack.", x),
        ld.report = paste0("ld.", x)
      )
    })
    colnames(strata.id) <- sn
  } else strata.id <<- NULL
})

# save gtypes object to working directory
observeEvent(input$save.gtypes, {
  isolate({
    if(!is.null(vals$gtypes) & !is.null(vals$wd) & !is.null(input$gtypesName)) {
      g.name <- input$gtypesName
      if(g.name != "") {
        g.name <- make.names(g.name)
        assign(g.name, vals$gtypes)
        fname <- paste0(g.name, "_gtypes.rdata")
        save(list = g.name, file = file.path(vals$wd, fname))
        showNotification("File saved", duration = 2, type = "message")
      }
    }
  })
})