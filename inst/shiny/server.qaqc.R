output$qaqcDiploidStep <- renderUI({
  list(
    ui.loci.missing, ui.loci.hmzg, ui.samples.missing, ui.samples.hmzg,
    ui.dups, ui.hwe.jack, ui.ld
  )[[qaqc.flow$step]]()
})

qaqc.flow <- reactiveValues(step = 1)
observeEvent(input$btn.qaqc.next, {
  if(qaqc.flow$step < 7) qaqc.flow$step <- qaqc.flow$step + 1
})
observeEvent(input$btn.qaqc.reset, qaqc.flow$step <- 1)

reloaded <- reactiveValues(count = 0)

output$gtypes.smry <- renderPrint({
  reloaded$count
  options(digits = 2)
  current.g
}, width = 120)

by.sample <- reactive({
  reloaded$count
  input$btn.run.loci.missing
  input$btn.run.loci.hmzg
  if(is.null(current.g)) NULL else summarizeSamples(current.g)
})

by.locus <- reactive({
  reloaded$count
  input$btn.run.samples.missing
  input$btn.run.samples.hmzg
  if(is.null(current.g)) NULL else {
    df <- data.frame(summarizeLoci(current.g))
    df$num.missing <- nInd(current.g) - df$num.genotyped
    df$pct.missing <- df$num.missing / nInd(current.g)
    df$pct.hmzg <- 1 - df$obsvd.heterozygosity
    df$num.hmzg <- df$pct.hmzg * df$num.genotyped
    cbind(locus = rownames(df), df)
  }
})

dups <- reactive({
  input$btn.run.dups
  if(is.null(current.g)) NULL else dupGenotypes(current.g, 0, detectCores(logical = F) - 1)
})

hw.jack <- reactive({
  input$btn.run.hwe.jack
  if(is.null(current.g)) NULL else jackHWE(current.g, show.progress = FALSE)
})

hwe <- reactive({
  input$btn.run.hwe
  if(is.null(current.g)) NULL else hweTest(current.g)
})

ld <- reactive({
  input$btn.run.ld
  if(is.null(current.g)) NULL else LDgenepop(current.g)
})
  
source("server.qaqc.01.loci.missing.R", local = TRUE)
source("server.qaqc.02.loci.hmzg.R", local = TRUE)
source("server.qaqc.03.samples.missing.R", local = TRUE)
source("server.qaqc.04.samples.hmzg.R", local = TRUE)
source("server.qaqc.05.dups.R", local = TRUE)
source("server.qaqc.06.hwe.jack.R", local = TRUE)
source("server.qaqc.07.ld.R", local = TRUE)
source("server.qaqc.remove.samples.loci.R", local = TRUE)
