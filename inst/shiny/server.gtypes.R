# directory for saved file
save.directory <- reactive({
  save.directory <- input$save.directory
  if (is.null(save.directory)) 
    return (NULL)
  save.directory 
})

# reads/stores genetic data file input 
gen.data <- reactive({
  gen.data <- input$genetic.data
  if(is.null(gen.data))
    return (NULL)
  gen.data <- read.csv(gen.data$datapath)
  gen.data
})

output$gen.data.head <- renderPrint({
  head(gen.data())
})
###End###

###makes the visibility of the file input for strata data conditional on dropdown menu###
output$ui_strata.upload <- renderUI({
  if (is.null(input$select.strata))
    return()
  
  switch(input$select.strata,
         "1" = fileInput("strata", label = "Upload separate strata file"),
         "2" = NULL)
})
###End###

###reads/stores strata data file input###  
strata.data <- reactive({
  strata.data <- input$strata
  if(is.null(strata.data))
    return (NULL)
  strata.data <- read.csv(strata.data$datapath)
  
  strata.schemes <- strata.data[, -1]
  rownames(strata.schemes) <- strata.data[, 1]
  strata.schemes
  
})

output$stats.strata <- renderPrint({
  head(strata.data())
})
###End###

###make the fasta file input conditional on the ploidy value###   
output$ui_1 <- renderUI({
  if (is.null(input$select.ploidy))
    return()
  
  switch(input$select.ploidy,
         "1" = fileInput("fasta", label = "If you have a FASTA file, you may upload it here"),
         "2" = NULL)
})
###End###

###reads/stores/prints fasta data input###
fasta.data <- reactive({
  fasta.data <- input$fasta
  if(is.null(fasta.data))
    return (NULL)
  fasta.data <- read.csv(fasta.data$datapath)
  fasta.data
})

output$stats.fasta <- renderPrint({
  head(fasta.data())
})
###End###

###reactive element that creats/displays gtypes object when gtypes button is pushed###
gtypes.object <- eventReactive(input$loadGtypes, {
  if(is.null(input$gtypesR)){
    genetic.data <- gen.data()
    FASTA.File <- if(is.null(input$fasta)==FALSE & input$select.ploidy==1)
    {fasta.data()} 
    
    STRATA.File <- if(is.null(input$strata)==FALSE & input$select.strata==1)
    {strata.data()}
    
    IDcol <- if(input$idCol>0){input$idCol}
    STRATAcol <- if(input$strataCol>0){input$strataCol}
    ploidy.num <- ifelse(input$select.ploidy==1, 1, 2)
    
    
    GTYPES <- df2gtypes(genetic.data, ploidy = ploidy.num, strata.col = STRATAcol,
                        id.col = IDcol, loc.col = input$lociCol, schemes = STRATA.File,
                        sequences = FASTA.File, description = input$description)
    
    if(is.null(schemes(GTYPES))==FALSE) {
      strata.options<- if(input$strataOptions.sep == "Do Not Stratify" | input$select.strata == 2)
      {NULL} else {input$strataOptions.sep}
      
      GTYPES.stratified <- stratify(GTYPES, strata.options)
      GTYPES.stratified
    } else {GTYPES}
  }
  else{
    gtypesRDS <- input$gtypesR
    gtypesRDS <- readRDS(gtypesRDS$datapath)
    
    if(is.null(schemes(gtypesRDS))==FALSE) {
      strata.options <- if(input$strataOptions.rds == "Do Not Stratify" | is.null(schemes(gtypesRDS)))
      {NULL} else {input$strataOptions.rds}
      
      gtypesRDS.stratified <- stratify(gtypesRDS, strata.options)
      gtypesRDS.stratified
    } else {gtypesRDS}
  }
})
#End###

###save gtypes as .rds file### 
save.gtypes <- observeEvent(input$save.gtypes, {
  if(!(is.null(gtypes.object()))){
    filename.gtypes <- paste(description(gtypes.object()),"rds", sep = ".")
    filename.gtypes <- paste(save.directory(), filename.gtypes, sep = "/")
    
    saveRDS(gtypes.object(), file = filename.gtypes)
  } else {NULL}
})
###End###

###Prints gtypes object###
output$stats_gtype <- renderPrint({  
  gtypes.object()
})
###End###

###option to stratify uploaded gtypes object that has not already been stratified### 
output$stratif.rds <- renderUI({
  if (is.null(input$gtypesR))
    return()
  
  uiOutput("strataColMenu.rds")
}) 

schemesdata.rds <- reactive({
  gtypesFile <- input$gtypesR
  gtypesFile <- readRDS(gtypesFile$datapath)
  schemesdata.rds <- schemes(gtypesFile)
  if (is.null(schemesdata.rds)) {
    return(NULL)
  }
  schemesdata.rds 
})

output$strataColMenu.rds <- renderUI({
  schemes.data <-schemesdata.rds()
  if (is.null(schemes.data)) return(NULL)
  
  items=names(schemes.data)
  first.item <- c("Do Not Stratify")
  names(items)=items
  strat.options <- c(first.item, items)
  selectInput("strataOptions.rds", "Stratification schemes",strat.options)
})
###End###

###option to stratify gtypes object when seperate strata file is uploaded###
output$stratif.sep <- renderUI({
  if (is.null(input$strata) | input$select.strata==2)
    return()
  
  uiOutput("strataColMenu.sep")
}) 

output$strataColMenu.sep <- renderUI({
  stratadata <- strata.data()
  if (is.null(stratadata) | input$select.strata==2) return(NULL)
  
  items=names(stratadata)
  first.item <- c("Do Not Stratify")
  names(items)=items
  strat.options <- c(first.item, items)
  selectInput("strataOptions.sep", "Stratification schemes",strat.options)
})
###End###
