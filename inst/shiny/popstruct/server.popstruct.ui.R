ui.haploid.metrics <- function() {
  wellPanel(
    checkboxInput("chi2", "Chi2", value = TRUE),
    checkboxInput("fst", "Fst", value = TRUE),
    checkboxInput("phist", "Phist", value = TRUE),
    conditionalPanel(
      "input.phist",
      hr(),
      checkboxInput("pairwise.deletion", "Pairwise deletion", value = TRUE),
      selectInput(
        "subModel",
        "Choose a substitution model",
        choices = c(
          "raw", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87",
          "T92", "TN93", "GG95", "logdet", "paralin"
        ),
        selected = "TN93"
      )
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
  if(is.null(vals$gtypes)) return()
  if(ploidy(vals$gtypes) == 1) ui.haploid.metrics() else ui.diploid.metrics()
})