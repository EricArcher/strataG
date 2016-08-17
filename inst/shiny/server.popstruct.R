
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
