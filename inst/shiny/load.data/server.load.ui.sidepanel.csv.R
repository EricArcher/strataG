uiLoadCsv <- function() {
  csv.fnames <- if(is.null(vals$wd)) NULL else {
    fnames <- dir(vals$wd, pattern = "(*.csv$)", ignore.case = TRUE)
    c("Choose a file" = "", fnames)
  }
  fasta.fnames <- if(is.null(vals$wd)) NULL else {
    fnames <- dir(vals$wd, pattern = "(*.fasta$)|(*.fas$)|(*.txt)", ignore.case = TRUE)
    c("Choose a file" = "", fnames)
  }
  
  verticalLayout(
    selectInput(
      "genDataFile", label = h4("Choose a .csv file of genetic data"),
      choices = csv.fnames
    ),
    uiOutput("col.select"),
    
    selectInput(
      "ploidy", label = h4("Select the ploidy"),
      choices = list("haploid" = 1, "diploid" = 2), selected = 2
    ),
    
    selectInput(
      "schemeFile", label = h4("Choose a .csv file of stratification schemes (optional)"),
      choices = csv.fnames
    ),
    
    conditionalPanel(
      "input.ploidy == '1'",
      selectInput(
        "fastaFile", label = h4("Choose a FASTA formatted file of sequences"),
        choices = fasta.fnames
      ),
      conditionalPanel(
        "input.fastaFile != ''",
        checkboxInput("labelHaps", "Label haplotypes", value = TRUE)
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