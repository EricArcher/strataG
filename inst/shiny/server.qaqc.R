output$qaqcDiploidStep <- renderUI({
  list(
    ui.loci.missing, ui.loci.hmzg, ui.samples.missing, ui.samples.hmzg,
    ui.dups, ui.hwe.jack, ui.ld
  )[[qaqc.flow$step]]()
})

observeEvent(input$btn.qaqc.next, {
  if(qaqc.flow$step < 7) qaqc.flow$step <- qaqc.flow$step + 1
})
observeEvent(input$btn.qaqc.reset, qaqc.flow$step <- 1)

by.sample <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    summarizeSamples(user.data$current.g)
  }
})

by.locus <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    df <- data.frame(summarizeLoci(user.data$current.g))
    df$num.missing <- nInd(user.data$current.g) - df$num.genotyped
    df$pct.missing <- df$num.missing / nInd(user.data$current.g)
    df$pct.hmzg <- 1 - df$obsvd.heterozygosity
    df$num.hmzg <- df$pct.hmzg * df$num.genotyped
    cbind(locus = rownames(df), df)
  }
})

dups <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    dupGenotypes(user.data$current.g, 0, detectCores(logical = F) - 1)
  }
})

hw.jack <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    jackHWE(user.data$current.g, show.progress = FALSE)
  }
})

hwe <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    hweTest(user.data$current.g)
  }
})

ld <- reactive({
  if(is.null(user.data$current.g)) NULL else {
    LDgenepop(user.data$current.g)
  }
})
  
source("server.qaqc.01.loci.missing.R", local = TRUE)
source("server.qaqc.02.loci.hmzg.R", local = TRUE)
source("server.qaqc.03.samples.missing.R", local = TRUE)
source("server.qaqc.04.samples.hmzg.R", local = TRUE)
source("server.qaqc.05.dups.R", local = TRUE)
source("server.qaqc.06.hwe.jack.R", local = TRUE)
source("server.qaqc.07.ld.R", local = TRUE)
source("server.qaqc.remove.samples.loci.R", local = TRUE)