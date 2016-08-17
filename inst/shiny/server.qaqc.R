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
