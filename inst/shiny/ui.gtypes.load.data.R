ui.gtypes.csv <- function() {
  wellPanel(
    fileInput(
      "gen.data", label = h4("Choose a .csv file of genetic data"),
      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
    ),
    
    selectInput(
      "ploidy", label = h4("Select the ploidy"),
      choices = list("haploid" = 1, "diploid" = 2), selected = 1
    ),
    
    fluidRow(
      helpText(h5("Select columns from genetic data. If column is not present, select 0.")),
      column(
        3,
        titlePanel(h4("Sample IDs")),
        numericInput("idCol", label = NULL, min = 0, value = 1)
      ),
      column(
        3,
        titlePanel(h4("Strata")),
        numericInput("strataCol", label = NULL, min = 0, value = 2)
      ),
      column(
        3,
        titlePanel(h4("Loci")),
        numericInput("lociCol", label = NULL, min = 1, value = 3)
      )
    ),
    
    fileInput(
      "schemes", label = h4("Choose a .csv file of stratification schemes (optional)"),
      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
    ),
    
    uiOutput("gtypesFasta"),
    
    textInput(
      "description",
      label = h4("Enter a description of the data"), value = ""
    ),
    
    textInput(
      "save.directory", 
      label = h4("Select a folder to save output to"), 
      value = "~/Desktop"
    ),
    
    actionButton("loadGtypes", label = "Load gtypes object")
  )
}

ui.gtypes.rdata <- function() {
  wellPanel(
    fileInput(
      "gtypesR", 
      label = h4("Load an existing gtypes object"), 
      accept = c('.rds', '.rda', '.rdata')
    )
  )
}