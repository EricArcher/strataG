shinyServer(function(input, output, session) {
  
  ###reads input for directory for saved files### 
  save.directory <- reactive({
    save.directory <- input$save.directory
    if (is.null(save.directory)) 
      return (NULL)
    save.directory 
  })
  ###End###
  
  ###reads/stores genetic data file input###   
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
  
  ###reactive interface for qaqc tab###
  #selection input for minimum percentage of shared loci for dup gen
  output$ui_dup.gen <- renderUI({
    if (is.null(input$checkbox.dup.gen))
      return()
    
    switch(input$checkbox.dup.gen,
           "TRUE" = numericInput("percent", label = h5("Choose Minimum Percentage of Shared Loci"), min = 0,
                                 value = 0.66, max = 1, step = 0.01),
           "FALSE" = NULL)
  })
  
  #selection input for minimum frequency for low freq sub
  output$ui_low.freq <- renderUI({
    if (is.null(input$checkbox6))
      return()
    
    switch(input$checkbox.low.freq,
           "TRUE" = numericInput("minimum.freq", label = h5("Minimum Frequency"), value = 1, min = 1),
           "FALSE" = NULL)
  })
  ###End###
  
  ###qaqc functions###
  #by sample summary 
  by.sample <- eventReactive(input$run_qaqc, {
    if(input$checkbox.by.sample==TRUE){
      sortStrata <- if(nStrata(gtypes.object()) > 1) {TRUE} else {FALSE}
      
      smry.sample <- summarizeSamples(gtypes.object(), sort.by.strata = sortStrata)
      smry.sample}
  })
  
  output$stats.by.sample <- renderPrint({  
    head(by.sample())
  })
  
  #by locus summary 
  by.locus <- eventReactive(input$run_qaqc, {
    if(input$checkbox.by.locus==TRUE){
      sortStrata <- if(nStrata(gtypes.object()) > 1) {TRUE} else {FALSE}
      
      if(sortStrata) {
        smry <- summarizeLoci(gtypes.object(), by.strata = sortStrata)
        hwe.p <- lapply(strataSplit(gtypes.object()), hweTest)
        sapply(names(smry), function(x) {
          cbind(smry[[x]], hwe.p = hwe.p[[x]])
        }, USE.NAMES = TRUE, simplify = FALSE)
      } else {
        smry <- summarizeLoci(gtypes.object(), by.strata = sortStrata)
        hwe.p <- hweTest(gtypes.object())
        cbind(smry, hwe.p = hwe.p)
      }
    }
  })
  
  output$stats.by.locus <- renderPrint({  
    head(by.locus())
  })

  #Duplicate genotypes
  duplicate.genotypes <- eventReactive(input$run_qaqc, {
    if(input$checkbox.dup.gen==TRUE){
      #Find samples that share alleles at 2/3rds of the loci
      dupGenotypes(gtypes.object(), num.shared = input$percent)}
  })
  
  output$stats.dup.gen <- renderPrint({  
    duplicate.genotypes()
  })
  
  #sequence Summary 
  output$ui_checkbox.seq.sum<-renderUI({
    if(is.null(sequences(gtypes.object()))==FALSE)
      return(checkboxInput("checkbox.seq.sum", "Sequence summary"))
    else (NULL)
  })
  
  sequence.summary <- eventReactive(input$run_qaqc, {
    if(input$checkboxseq.sum==TRUE){
      
      seq.smry <- summarizeSeqs(sequences(gtypes.object()))
      seq.smry} else {NULL}
  })
  
  output$stats.seq.sum <- renderPrint({  
    head(sequence.summary())
  })
  
  #Haplotype likelihood
  output$ui_checkbox.hap.like<-renderUI({
    if(is.null(sequences(gtypes.object()))==FALSE)
      return(checkboxInput("checkbox.hap.like", "Haplotype likelihood"))
    else (NULL)
  })
  
  haplotype.likelihood <- eventReactive(input$run_qaqc, {
    if(input$checkbox.hap.like==TRUE){
      
      haplotypeLikelihoods(sequences(gtypes.object()))
    } else {NULL}
  })
  
  output$plot.hap.like <-renderPlot({
    haplotype.likelihood()
  })
  
  output$stats.hap.like <-renderPrint({
    haplotype.likelihood()
  })
  
  #Low frequency substitutions
  output$ui_checkbox.low.freq<-renderUI({
    if(is.null(sequences(gtypes.object()))==FALSE)
      return(checkboxInput("checkbox.low.freq", "Low frequency substitutions"))
    else (NULL)
  })
  
  low.frequency.substitutions <- eventReactive(input$run_qaqc, {
    if(input$checkbox.low.freq==TRUE){
      
      lowFreqSubs(sequences(gtypes.object()), min.freq = input$minimum.freq)
    } else {NULL}
  })
  
  output$stats.low.freq <- renderPrint({  
    low.frequency.substitutions()
  })
  ###END###
  
  ###save button for qaqc###
  save.qaqc <- observeEvent(input$save.qaqc, {
    filename.by.sample <- paste(description(gtypes.object()),"bysample.csv", sep = ".")
    filename.by.sample <- paste(save.directory(), filename.by.sample, sep = "/")
    filename.by.locus <- paste(description(gtypes.object()),"bylocus.csv", sep = ".")
    filename.by.locus <- paste(save.directory(), filename.by.locus, sep = "/")
    
    filename.dup.gen <- paste(description(gtypes.object()),"dupgen.csv", sep = ".")
    filename.dup.gen <- paste(save.directory(), filename.dup.gen, sep = "/")
    filename.seqsum <- paste(description(gtypes.object()),"seqsum.csv", sep = ".")
    filename.seqsum <- paste(save.directory(), filename.seqsum, sep = "/")
    filename.haplike <- paste(description(gtypes.object()),"haplike.csv", sep = ".")
    filename.haplike <- paste(save.directory(), filename.haplike, sep = "/")
    filename.lfs <- paste(description(gtypes.object()),"lowfreq.csv", sep = ".")
    filename.lfs <- paste(save.directory(), filename.lfs, sep = "/")
    
    
    checkbox.seq.sum <- if(is.null(sequences(gtypes.object()))) {FALSE} else {input$checkbox4}
    checkbox.hap.like <- if(is.null(sequences(gtypes.object()))) {FALSE} else {input$checkbox5}
    checkbox.low.freq <- if(is.null(sequences(gtypes.object()))) {FALSE} else {input$checkbox6}
    
    if(input$checkbox.by.sample==TRUE){
      write.csv(by.sample(), file = filename.by.sample)}
    if(input$checkbox.by.locus==TRUE){
      write.csv(by.locus(), file = filename.by.locus)}
    if(input$checkbox.dup.gen==TRUE){
      write.csv(duplicate.genotypes(), file= filename.dup.gen)}
    if(checkbox.seq.sum==TRUE){
      write.csv(sequence.summary(), file = filename.seqsum)}
    if(checkbox.hap.like==TRUE){
      write.csv(haplotype.likelihood(), file = filename.haplike)}
    if(checkbox.low.freq==TRUE){
      write.csv(low.frequency.substitutions(), file = filename.lfs)}
  })
  ###End###
  
  ###displays for the output metrics on population structure tab###
  output$ui_phist<-renderUI({
    if(ploidy(gtypes.object())==1)
      return(checkboxInput("phist", "statPhist"))
  })
  
  output$ui_d<-renderUI({
    if(ploidy(gtypes.object())==2)
      return(checkboxInput("d", "statJostD"))
  })
  
  output$ui_fst.prime<-renderUI({
    if(ploidy(gtypes.object())==2)
      return(checkboxInput("fst.prime", "statFstPrime"))
  })
  
  output$ui_fis<-renderUI({
    if(ploidy(gtypes.object())==2)
      return(checkboxInput("fis", "statFis"))
  })
  
  output$ui_gst<-renderUI({
    if(ploidy(gtypes.object())==2)
      return(checkboxInput("gst", "statGst"))
  })
  
  output$ui_gst.prime<-renderUI({
    if(ploidy(gtypes.object())==2)
      return(checkboxInput("gst.prime", "statGstPrime"))
  })
  
  output$ui_gst.dbl.prime<-renderUI({
    if(ploidy(gtypes.object())==2)
      return(checkboxInput("gst.dbl.prime", "statGstDblPrime"))
  })
  ###End###
  
  ###Display the option to set pairwise deletion true/false###
  pairwise.deletion <- reactive({
    pairwise.deletion <- input$pairwise.deletions
    if (is.null(pairwise.deletion)) 
      return (NULL)
    pairwise.deletion <- if(pairwise.deletion==1){TRUE} else {FALSE}
    pairwise.deletion
  })
  ###End###
  
  ###Calcuates overall population structure###
  overall.pop <- eventReactive(input$run_pop.struc, {
    if(input$overall.checkbox==TRUE & ploidy(gtypes.object())==2){
      chi2.list <- if(input$chi2==TRUE){"Chi2"} else{NULL}
      fst.list <- if(input$fst==TRUE){"Fst"} else{NULL}
      d.list <- if(input$d==TRUE){"D"} else{NULL}
      fst.prime.list <- if(input$fst.prime==TRUE){"Fst.Prime"} else{NULL}
      fis.list <- if(input$fis==TRUE){"Fis"} else{NULL}
      gst.list <- if(input$gst==TRUE){"Gst"} else{NULL}
      gst.prime.list <- if(input$gst.prime==TRUE){"Gst.Prime"} else{NULL}
      gst.dbl.prime.list <- if(input$gst.dbl.prime==TRUE){"Gst.Dbl.Prime"} else{NULL}
      
      function.list.2 <- c(chi2.list, fst.list,
                           d.list, fst.prime.list, fis.list,
                           gst.list, gst.prime.list, gst.dbl.prime.list)
      
      overall.test<<- overallTest(gtypes.object(), stats = function.list.2, nrep = input$permutation, 
                                  pairwise.deletion = pairwise.deletion)}
    
    if(input$overall.checkbox==TRUE & ploidy(gtypes.object())==1){
      chi2.list <- if(input$chi2==TRUE){"Chi2"} else{NULL}
      fst.list <- if(input$fst==TRUE){"Fst"} else{NULL}
      phist.list <- if(input$phist==TRUE){"PHIst"} else{NULL}
      
      function.list.1 <- c(chi2.list, fst.list, phist.list)
      
      overall.test <<- overallTest(gtypes.object(), stats = function.list.1, nrep = input$permutation,
                                   pairwise.deletion = pairwise.deletion)}
    
    if(input$overall.checkbox==FALSE) {
      overall.test <<- NULL}
  })
  
  output$stats_overall.pop<-renderPrint({
    overall.pop()
  })
  ###End###
  
  ###Calculates pairwise population structure### 
  pairwise.pop <-eventReactive(input$run_pop.struc, {
    if(input$pairwise.checkbox==TRUE & ploidy(gtypes.object())==2){
      chi2.list <- if(input$chi2==TRUE){"Chi2"} else{NULL}
      fst.list <- if(input$fst==TRUE){"Fst"} else{NULL}
      d.list <- if(input$d==TRUE){"D"} else{NULL}
      fst.prime.list <- if(input$fst.prime==TRUE){"Fst.Prime"} else{NULL}
      fis.list <- if(input$fis==TRUE){"Fis"} else{NULL}
      gst.list <- if(input$gst==TRUE){"Gst"} else{NULL}
      gst.prime.list <- if(input$gst.prime==TRUE){"Gst.Prime"} else{NULL}
      gst.dbl.prime.list <- if(input$gst.dbl.prime==TRUE){"Gst.Dbl.Prime"} else{NULL}
      
      function.list.2 <- c(chi2.list, fst.list,
                           d.list, fst.prime.list, fis.list,
                           gst.list, gst.prime.list, gst.dbl.prime.list)
      
      pairwise.test <<- pairwiseTest(gtypes.object(), function.list.2, nrep = input$permutation,
                                     pairwise.deletion = pairwise.deletion)}
    
    if(input$pairwise.checkbox==TRUE & ploidy(gtypes.object())==1){
      chi2.list <- if(input$chi2==TRUE){"Chi2"} else{NULL}
      fst.list <- if(input$fst==TRUE){"Fst"} else{NULL}
      phist.list <- if(input$phist==TRUE){"PHIst"} else{NULL}
      
      function.list.1 <- c(chi2.list, fst.list, phist.list)
      
      pairwise.test <<- pairwiseTest(gtypes.object(), function.list.1, nrep = input$permutation,
                                     pairwise.deletion = pairwise.deletion)}
    
    if(input$pairwise.checkbox==FALSE) {
      pairwise.test <<- NULL}
  })
  
  output$stats_pairwise.pop<-renderPrint({
    pairwise.pop()
  })
  ###End###
  
  ###save button for pop structure###
  save.pop.structure <- observeEvent(input$save.pop.structure, {
    filename.overall <- paste(description(gtypes.object()),"overall.csv", sep = ".")
    filename.overall <- paste(save.directory(), filename.overall, sep = "/")
    filename.pairwise <- paste(description(gtypes.object()),"pairwise.csv", sep = ".")
    filename.pairwise <- paste(save.directory(), filename.pairwise, sep = "/")
    
    if(input$overall.checkbox==TRUE){
      write.csv(overall.test$result, file = filename.overall)}
    if(input$pairwise.checkbox==TRUE){
      write.csv(pairwise.test$result, file = filename.pairwise)}
  })
  ###End###
  
})