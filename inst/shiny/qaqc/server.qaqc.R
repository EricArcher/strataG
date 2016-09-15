output$flowStep <- reactive(vals$qaqc.step)

output$stepLabel <- renderUI({
  step.title <- switch(
    vals$qaqc.step,
    "1) Percent of Missing Loci per Sample",
    "2) Percent of Homozygous Loci per Sample",
    "3) Percent of Missing Samples per Locus",
    "4) Percent of Homozygous Samples per Locus",
    "5) Duplicate Check",
    "6) Hardy-Weinberg Equilibrium",
    "7) Linkage Disequilibrium",
    "8) QA/QC Summary Reports"
  )
  
  step.help <- switch(
    vals$qaqc.step,
    "Identify samples with too many missing loci by moving the slider to set a 
    threshold. Pressing the [Remove samples] button will remove samples that are 
    at or exceed this threshold and appear in the table.",
    "Identify samples with too many homozygous loci by moving the slider to set a 
    threshold. Pressing the [Remove samples] button will remove samples that are 
    at or exceed this threshold and appear in the table.",
    "Identify loci with too many missing samples by moving the slider to set a 
    threshold. Pressing the [Remove loci] button will remove loci that are 
    at or exceed this threshold and appear in the table.",
    "Identify loci with too many homozygous samples by moving the slider to set a 
    threshold. Pressing the [Remove loci] button will remove loci that are 
    at or exceed this threshold and appear in the table.",
    "Identify samples that share too many loci by moving the slider to set a 
    threshold.",
    "Identify samples that are influential to Hardy-Weinberg Equilibrium (HWE). 
    Use the slider to set a critical alpha level below which loci are considered to 
    be out of HWE.",
    "Identify pairs of loci that are out of linkage disequilibrium (LD). Use the slider 
    to set a critical alpha level below which two loci are considered to be out of LD.",
    "View summary reports from the QA/QC steps and save them as .csv files"
  )
  
  verticalLayout(
    wellPanel(
      verticalLayout(
        titlePanel(h4(step.title)),
        helpText(step.help),
        uiOutput("nextButton")
      )
    ),
    hr()
  )
})

output$qaqcHaploidStep <- renderUI({
  if(is.null(vals$qaqc.step)) return()
})
  
output$qaqcDiploidStep <- renderUI({
  if(is.null(vals$qaqc.step)) return()
  switch(
    vals$qaqc.step,
    ui.loci.missing, ui.loci.hmzg, ui.samples.missing, ui.samples.hmzg,
    ui.dups, ui.hwe.jack, ui.ld, ui.reports
  )()
})

output$nextButton <- renderUI({
  if(vals$qaqc.step == 8) return()
  actionButton("btn.qaqc.next", label = "Next")
})

observeEvent(input$btn.qaqc.next, {
  isolate({
    if(vals$qaqc.step < 8) {
      vals$qaqc.step <- vals$qaqc.step + 1
      # if(vals$qaqc.step == 5) qaqc.reports$dup <<- NULL
      # if(vals$qaqc.step == 6) qaqc.reports$hwe.jack <<- NULL
      # if(vals$qaqc.step == 7) qaqc.reports$ld <<- NULL
      first.run <<- TRUE
    }
  })
})

output$dt.by.sample <- renderDataTable({
  df <- qaqc.reports$samples
  if(!is.null(df)) df <- round(df, 4)
  df <- switch(
    input$sampleSelect,
    all = df,
    not.removed = df[is.na(df$step.removed), ],
    removed = df[!is.na(df$step.removed), ]
  )
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, mode = "multiple", target = "row", scrollX = TRUE)
  )
})

output$dt.by.locus <- renderDataTable({
  df <- qaqc.reports$loci
  if(!is.null(df)) df <- round(df, 4)  
  df <- switch(
    input$locusSelect,
    all = df,
    not.removed = df[is.na(df$step.removed), ],
    removed = df[!is.na(df$step.removed), ]
  )
  DT::datatable(
    df, rownames = FALSE,
    options = list(paging = nrow(df) > 10, mode = "multiple", target = "row", scrollX = TRUE)
  )
})

observeEvent(input$btn.remove.samples, {
  i <- input$dt.by.sample_rows_selected
  if(!is.null(i)) {
    df <- by.sample()
    id <- df$id[i]
    all.inds <- indNames(vals$gtypes)
    to.keep <- setdiff(all.inds, id)
    if(length(to.keep) > 0) {
      vals$gtypes <- vals$gtypes[to.keep, , ]
    }
  }
})

observeEvent(input$btn.remove.loci, {
  i <- input$dt.by.locus_rows_selected
  if(!is.null(i)) {
    df <- by.locus()
    loc <- df$locus[i]
    all.loci <- locNames(vals$gtypes)
    to.keep <- setdiff(all.loci, loc)
    if(length(to.keep) > 0) {
      vals$gtypes <- vals$gtypes[ ,to.keep , ]
    }
  }
})

source("server.qaqc.summaries.R", local = TRUE)
source("server.qaqc.01.loci.missing.R", local = TRUE)
source("server.qaqc.02.loci.hmzg.R", local = TRUE)
source("server.qaqc.03.samples.missing.R", local = TRUE)
source("server.qaqc.04.samples.hmzg.R", local = TRUE)
source("server.qaqc.05.dups.R", local = TRUE)
source("server.qaqc.06.hwe.jack.R", local = TRUE)
source("server.qaqc.07.ld.R", local = TRUE)
source("server.qaqc.08.reports.R", local = TRUE)