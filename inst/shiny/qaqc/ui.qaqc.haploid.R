qaqcHaploid <- function() {
  tabsetPanel(
    tabPanel(
      "Sequence summary",
      splitLayout(
        verticalLayout(
          titlePanel(h4("By strata summary")), 
          dataTableOutput("seqSmryTable"),
          actionButton("saveSeqSmry", "Save table")
        ),
        verticalLayout(
          titlePanel(h4("Frequency table")), 
          dataTableOutput("hapFreqTable"),
          actionButton("saveHapFreq", "Save table")
        )
      )
    ),
    tabPanel(
      "Sequence likelihoods",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("smryPwsDelete", "Pairwise deletion", value = TRUE),
          selectInput(
            "smrySubModel",
            "Choose a substitution model",
            choices = c(
              "raw", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87",
              "T92", "TN93", "GG95", "logdet", "paralin"
            ),
            selected = "raw"
          ),
          numericInput(
            "smrySeqLikeN", "Number of sequences to plot", 
            value = 20, min = 5, step = 1
          ),
          actionButton("saveSeqLike", "Save likelihood table")
        ),
        mainPanel(plotOutput("seqLikePlot"))
      )
    ),
    tabPanel(
      "Low frequency substitutions",
      sidebarLayout(
        sidebarPanel(
          numericInput(
            "minFreq", "Minimum frequency of base to be flagged",
            value = 3, min = 1, step = 1
          ),
          numericInput(
            "motifLength", "Length of motif around low frequency base",
            value = 10, min = 4, step = 2
          ),
          actionButton("saveLowFreqSubs", "Save table")
        ),
        mainPanel(dataTableOutput("lowFreqSubTable"))
      )
    )
  )
}