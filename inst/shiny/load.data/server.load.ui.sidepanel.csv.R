uiLoadCsv <- function() {
  verticalLayout(
    fileInput(
      "gen.data", label = h4("Choose a .csv file of genetic data"),
      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
    ),
    uiOutput("col.select"),
    
    selectInput(
      "ploidy", label = h4("Select the ploidy"),
      choices = list("haploid" = 1, "diploid" = 2), selected = 2
    ),
    
    fileInput(
      "schemes", 
      label = h4("Choose a .csv file of stratification schemes (optional)"),
      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
    ),
    
    conditionalPanel(
      "input.ploidy == '1'",
      fileInput(
        "fasta", label = h4("Choose a FASTA formatted file of sequences"),
        accept = c(".fasta", ".fas", ".txt")
      )
    ),
    
    textInput(
      "description", label = h4("Enter a description of the data"), value = ""
    ),
    
    actionButton("loadCsvGtypes", label = "Load gtypes object")
  )
}

output$col.select <- renderUI({
  df <- gen.data.df()
  cols <- NULL
  if(!is.null(df)) {
    cols <- 1:length(colnames(df))
    names(cols) <- colnames(df)
  }
  fluidRow(
    column(
      4,
      selectInput("idCol", label = h5("Sample ID column"), c("N/A" = 0, cols), selected = cols[1])
    ),
    column(
      4,
      selectInput("strataCol", label = h5("Stratification column"), c("N/A" = 0, cols), selected = cols[2])
    ),
    column(
      4,
      selectInput("lociCol", label = h5("First locus column"), cols, selected = cols[3])
    )
  )
})