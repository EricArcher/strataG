ui.load.csv <- function() {
  wellPanel(
    fileInput(
      "gen.data", label = h4("Choose a .csv file of genetic data"),
      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
    ),
    helpText(h5("Select columns from genetic data. If column is not present, select 0.")),
    fluidRow(
      column(
        4,
        titlePanel(h5("Sample IDs")),
        numericInput("idCol", label = NULL, min = 0, value = 1)
      ),
      column(
        4,
        titlePanel(h5("Strata")),
        numericInput("strataCol", label = NULL, min = 0, value = 2)
      ),
      column(
        4,
        titlePanel(h5("Loci")),
        numericInput("lociCol", label = NULL, min = 1, value = 3)
      )
    ),
    hr(),
    
    selectInput(
      "ploidy", label = h4("Select the ploidy"),
      choices = list("haploid" = 1, "diploid" = 2), selected = 2
    ),
    hr(),
    
    fileInput(
      "schemes", 
      label = h4("Choose a .csv file of stratification schemes (optional)"),
      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
    ),
    hr(),
    
    conditionalPanel(
      "input.ploidy == '1'",
      fileInput(
        "fasta", label = h4("Choose a FASTA formatted file of sequences"),
        accept = c(".fasta", ".fas", ".txt")
      ),
      hr()
    ),
    
    textInput(
      "description", label = h4("Enter a description of the data"), value = ""
    ),
    hr(),
    
    actionButton("loadGtypes", label = "Load gtypes object")
  )
}

ui.load.rdata <- function() {
  wellPanel(
    fileInput(
      "gtypesR", 
      label = h4("Load an existing gtypes object"), 
      accept = c('.rds', '.rda', '.rdata')
    )
  )
}