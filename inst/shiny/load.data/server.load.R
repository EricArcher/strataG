source("server.load.ui.R", local = TRUE)

gen.data.df <- reactive({
  x <- input$genDataFile
  if(x == "") NULL else readGenData(file.path(vals$wd, x))
})

schemes.df <- reactive({
  x <- input$schemeFile
  if(x == "") NULL else {
    df <- readGenData(file.path(vals$wd, x))
    rownames(df) <- df[, 1]
    df
  }
})

seqs.DNAbin <- reactive({
  x <- input$fastaFile
  if(x == "") NULL else read.fasta(file.path(vals$wd, x))
})


# tab 1: view genetic data in table
output$file.gen.data <- renderDataTable({
  df <- gen.data.df()
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, scrollX = TRUE)
  )
})

# tab 2: view stratification schemes in table
output$file.schemes <- renderDataTable({
  df <- schemes.df()
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, scrollX = TRUE)
  )
})

# tab 3: view sequences in table
output$file.fasta <- renderPrint(seqs.DNAbin())

gtypesLoaded <- function() {      
  vals$qaqc.step <- 1
  first.run <<- TRUE
  if(ploidy(vals$gtypes) > 1) {
    vals$qaqc.reports <<- list()
    vals$qaqc.reports$samples <<- by.sample()
    vals$qaqc.reports$samples$threshold <<- vals$qaqc.reports$samples$step.removed <<- NA
    vals$qaqc.reports$loci <<- by.locus()
    vals$qaqc.reports$loci$threshold <<- vals$qaqc.reports$loci$step.removed <<- NA 
  }
}

# load gtypes - responds to loadCsvGtypes actionButton
observeEvent(input$loadCsvGtypes, {
  isolate({
    gen.data <- gen.data.df()
    vals$gtypes <- if(is.null(gen.data)) NULL else {
      g <- df2gtypes(
        gen.data,
        ploidy = as.numeric(input$ploidy),
        id.col = if(as.numeric(input$idCol) == 0) NULL else as.numeric(input$idCol),
        strata.col = if(as.numeric(input$strataCol) == 0) NULL else as.numeric(input$strataCol),
        loc.col = as.numeric(input$lociCol),
        sequences = if(as.numeric(input$ploidy) > 1) NULL else seqs.DNAbin(),
        schemes = schemes.df(),
        description = if(input$description == "") NULL else input$description
      )
      if(ploidy(g) == 1 & !is.null(sequences(g)) & input$labelHaps) {
        labelHaplotypes(g)$gtypes
      } else g
    }
    if(!is.null(vals$gtypes)) gtypesLoaded()
  })
})



observe({
  rm(list = ls(envir = vals$scratch.env), envir = vals$scratch.env)
  fname <- as.character(input$slctRdataFile)
  if(length(fname) == 1) {
    if(fname != "") {
      load(file.path(vals$wd, fname), envir = vals$scratch.env)
      objs <- ls(envir = vals$scratch.env)
      valid.gtypes <- sapply(objs, function(x) {
        is.gtypes(eval(parse(text = x), envir = vals$scratch.env))
      })
      rm(list = objs[!valid.gtypes], envir = vals$scratch.env)
      objs <- objs[valid.gtypes]
      if(length(objs) > 0) {
        obj.list <- as.list(objs)
        names(obj.list) <- objs
        objs <- c("Select a gtypes object" = "", objs)
        updateSelectInput(session, "slctGtypes", choices = objs)
      }
    }
  }
})

observeEvent(input$slctGtypes, {
  isolate({
    if(!is.null(input$slctGtypes)) {
      if(input$slctGtypes != "") {
        vals$gtypes <- get(input$slctGtypes, envir = vals$scratch.env)
        gtypesLoaded()
      }
    }
  })
})