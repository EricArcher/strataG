source("server.load.ui.R", local = TRUE)

shinyDirChoose(
  input, "wd", roots = volumes, 
  session = session, restrictions = system.file(package = "base")
)
output$wdPath <- renderPrint(cat(parseDirPath(volumes, input$wd)))

gen.data.df <- reactive({
  x <- input$gen.data
  if(is.null(x)) NULL else readGenData(x$datapath)
})

schemes.df <- reactive({
  x <- input$schemes
  if(is.null(x)) NULL else {
    df <- readGenData(x$datapath)
    rownames(df) <- df[, 1]
    df
  }
})

seqs.DNAbin <- reactive({
  x <- input$fasta
  if(is.null(x)) NULL else read.fasta(x$datapath)
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

# load gtypes - responds to loadGtypes actionButton
observe({
  input$loadGtypes
  isolate({
    gen.data <- gen.data.df()
    user.data$current.g <- if(is.null(gen.data)) NULL else {
      df2gtypes(
        gen.data,
        ploidy = as.numeric(input$ploidy),
        id.col = if(input$idCol == 0) NULL else input$idCol,
        strata.col = if(input$strataCol == 0) NULL else input$strataCol,
        loc.col = input$lociCol,
        sequences = if(as.numeric(input$ploidy) > 1) NULL else seqs.DNAbin(),
        schemes = schemes.df(),
        description = if(input$description == "") NULL else input$description
      )
    }
    if(!is.null(user.data$current.g)) {
      updateTabsetPanel(session, "main.tab", "Loaded gtypes")
    }
  }) 
})