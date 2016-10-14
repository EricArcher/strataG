by.sample <- reactive({
  if(is.null(vals$gtypes)) NULL else summarizeSamples(vals$gtypes)
})

by.locus <- reactive({
  if(is.null(vals$gtypes)) NULL else {
    df <- data.frame(summarizeLoci(vals$gtypes))
    df$num.missing <- nInd(vals$gtypes) - df$num.genotyped
    df$pct.missing <- df$num.missing / nInd(vals$gtypes)
    df$pct.hmzg <- 1 - df$obsvd.heterozygosity
    df$num.hmzg <- df$pct.hmzg * df$num.genotyped
    cbind(locus = rownames(df), df)
  }
})

dups <- reactive({
  if(is.null(vals$gtypes)) NULL else {
    if(!is.null(id)) removeNotification(id)
    id <<- showNotification(
      "Looking for duplicates...", duration = NULL, closeButton = FALSE,
      type = "message"
    )
    df <- dupGenotypes(vals$gtypes, 0, detectCores(logical = F) - 1)
    removeNotification(id)
    if(first.run) {
      vals$qaqc.reports$dups <<- df
      first.run <<- FALSE
    }
    df
  }
})

hwe.p.mat <- function(x) {
  df <- t(sapply(rev(p.adjust.methods), function(m) p.adjust(x, method = m)))
  data.frame(method = rownames(df), df)
}

loadHweReportTabs <- function(x) {
  df <- vals$qaqc.reports$hwe.jack[[x]]$p.mat
  output[[strata.id["hwe.report", x]]] <- renderDataTable({ 
    DT::datatable(
      df, rownames = FALSE,
      selection = "none",
      options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
    )
  })
  
  # format jackknife table
  infl.df <- vals$qaqc.reports$hwe.jack[[x]]$inf
  output[[strata.id["hwe.jack.report", x]]] <- renderDataTable({
    DT::datatable(
      infl.df, rownames = FALSE,
      selection = "none",
      options = list(paging = nrow(infl.df) > 10, searching = FALSE, scrollX = TRUE)
    )
  })
}

hw.jack <- reactive({
  isolate({
    if(is.null(vals$gtypes)) NULL else {
      if(!is.null(id)) removeNotification(id)
      id <<- showNotification(
        "Running HWE jackknife...", duration = NULL, closeButton = FALSE,
        type = "message"
      )
      jack.list <- lapply(
        strataSplit(vals$gtypes), jackHWE, 
        use.genepop = input$hwe.source == 2, show.progress = FALSE
      )
      removeNotification(id)
      if(first.run) {
        vals$qaqc.reports$hwe.jack <- sapply(jack.list, function(x) {
          list(p.mat = hwe.p.mat(x$obs), infl = jackInfluential(x)$influential)
        }, simplify = FALSE, USE.NAMES = TRUE)
        for(x in colnames(strata.id)) loadHweReportTabs(x)
        first.run <<- FALSE
      }
      jack.list
    }
  })
})

loadLdReportTabs <- function(x) {
  output[[strata.id["ld.report", x]]] <- renderDataTable({ 
    DT::datatable(
      vals$qaqc.reports$ld[[x]], rownames = FALSE,
      selection = "none",
      options = list(paging = nrow(df) > 10, searching = FALSE, scrollX = TRUE)
    )
  })
}

ld <- reactive({
  isolate({
    if(is.null(vals$gtypes)) NULL else {
      if(!is.null(id)) removeNotification(id)
      id <<- showNotification(
        "Calculating linkage disequilibrium...", duration = NULL, 
        closeButton = FALSE, type = "message"
      )
      ld.list <- lapply(strataSplit(vals$gtypes), LDgenepop)
      removeNotification(id)
      if(first.run) {
        vals$qaqc.reports$ld <- ld.list
        for(x in colnames(strata.id)) loadLdReportTabs(x)
        first.run <<- FALSE
      }
      ld.list
    }
  })
})