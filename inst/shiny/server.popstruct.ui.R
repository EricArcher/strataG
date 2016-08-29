ui.haploid.metrics <- function() {
  wellPanel(
    checkboxInput("chi2", "Chi2", value = TRUE),
    checkboxInput("fst", "Fst", value = TRUE),
    checkboxInput("phist", "Phist", value = TRUE),
    conditionalPanel(
      "input.phist == TRUE",
      hr(),
      checkboxInput("pairwise.deletion", "Pairwise deletion", value = TRUE)
    )
  )
}

ui.diploid.metrics <- function() {  
  wellPanel(
    checkboxInput("chi2", "Chi2", value = TRUE),
    checkboxInput("fst", "Fst", value = TRUE),
    checkboxInput("fst.prime", "F'st", value = TRUE),
    checkboxInput("gst", "Gst", value = TRUE),
    checkboxInput("gst.prime", "G'st", value = TRUE),
    checkboxInput("gst.dbl.prime", "G''st", value = TRUE),
    checkboxInput("jost.d", "Jost's D", value = TRUE),
    checkboxInput("fis", "Fis", value = TRUE)
  )
}

output$metrics <- renderUI({
  if(length(input$ploidy) == 0) return()
  if(input$ploidy == "1") ui.haploid.metrics() else ui.diploid.metrics()
})