###Sidebar panel for population structure tab###
tabPanel(
  title = "Population Structure",
  titlePanel("Population Structure"),
  helpText("This tab allows you perform population structure analyses. First choose to run the population structure
            tests in an overall and/or pairwise manner by checking the box or boxes in the input panel. 
            Then select the specific metric(s) to be computed and specify the number of permutations in the box below.
            Click the Run button to analyze the data and display the outputs in the viewer tabs. You may save the 
            outputs to your specified directory in .csv format by clicking the Save button."),
  sidebarLayout(
    sidebarPanel(
      h4("Types"),
      checkboxInput("overall.checkbox", "Overall"),
      checkboxInput("pairwise.checkbox", "Pairwise"),
      titlePanel(h4("Metrics")),
      checkboxInput("chi2", "statChi2"),
      checkboxInput("fst", "statFst"),
      uiOutput("ui_phist"),
      uiOutput("ui_d"),
      uiOutput("ui_fst.prime"),
      uiOutput("ui_fis"),
      uiOutput("ui_gst"),
      uiOutput("ui_gst.prime"),
      uiOutput("ui_gst.dbl.prime"),
      selectInput(
        "pairwise.deletions",
        label = h5("Pairwise Deletion"),
        choices = list("Yes" = 1, "No" = 2),
        selected = 1),
      numericInput(
        "permutation",
        label = h5("Number of permutations"), 
        min = 0, value = 1000),
      actionButton("run_pop.struc", label = "Run"),
      actionButton("save.pop.structure", label = "Save")),
           
      # Main panel for the population structure tab that display pop structure data###
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel("Overall", verbatimTextOutput("stats_overall.pop")),
          tabPanel("Pairwise", verbatimTextOutput("stats_pairwise.pop"))
        )
      )
    )
)
