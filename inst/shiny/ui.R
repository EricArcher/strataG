library(shiny)

shinyUI(fluidPage(
  
  tabsetPanel(
    
    ###Upper sections of gtypes panel that gives overall directions and asks user to upload###
    ###already existing gtypes object###
    tabPanel("gtypes",
             
             titlePanel("gtypes Object"),
             
             helpText("This application will allow you to analyze your data with functions of the R package strataG.",
                      "To use these functions, you must first upload or create a gtypes object."),
             
             helpText("If you already have a gtypes object saved as a RData file, you may upload it here. If a 
                      seperate strata file was used to created your uploaded gtypes object, you will
                      be prompted to select a stratification scheme. Click the Load gtypes button at the buttom of the page 
                      to complete the upload and display your object in the gtypes tab below. 
                      Please note that you will have to refresh the app if you wish to make a new
                      gtypes object after loading an existing one here."),
             
             fileInput("gtypesR", label = "Choose an existing gtypes object saved as an RData file", accept = c('.rds', '.rda', '.RData')),
             
             #display different stratification schemes if a seperate strata file was used to created uploaded gtypes
             uiOutput("stratif.rds"),
             
             
             helpText("If you do not yet have a gtypes object, input your raw data below to create a new gtypes object.
                      It is important here to note a few formatting requirements for the raw data:"),
             helpText("1. Raw genetic data files and strata files must be in .csv format and sequence data must be in FASTA format"),
             helpText("2. If there are strata columns in the genetic data file, the strata columns must come before the loci columns"),
             helpText("Summaries of the uploaded raw data will display in the tabs below to guide you as you select the appropriate columns
                      for your data. Click the Load gtypes object button at the end of the input panel to make your gtypes 
                      object appear in the gtypes tab."),
             br(),
             ###End###
             
             
             ###Side bar panel that allows user to create a gtypes object from uploaded data###
             sidebarLayout(
               sidebarPanel(
                 
                 
                 
                 textInput("save.directory", label = h5("Input the datapath for the directory where you would like to save any ouputted files. This will default
                                                        to your Desktop if you do not specify a different datapath."), 
                           value = "~/Desktop"),
                 
                 textInput("description", label = h5("Provide a description of the inputted data. This description will serve as",
                                                     "the title of the gtypes object and of any subsequent files created"), 
                           value = ""),
                 
                 
                 
                 fileInput('genetic.data', 'Choose a .csv file of your genetic data',
                           accept=c('text/csv', 
                                    'text/comma-separated-values,text/plain', 
                                    '.csv')),
                 helpText("Look at your genetic data file in the viewer to the right and select the appropriate column numbers for your data"),
                 
                 numericInput("idCol", label = h5("Select the column containing your data's IDs. If there is no ID column, select 0."), 
                              min = 0, value = 1),
                 
                 
                 numericInput("lociCol", label = h5("Select the first column containing loci"), 
                              min = 1, value = 1),
                 
                 
                 numericInput("strataCol", 
                              label = h5("Select the column in your genetic data file containing the stratification scheme you wish to examine",
                                         "in your gtypes object. If there are no strata columns, select 0"),
                              min = 0, value = 1),
                 
                 
                 selectInput("select.strata", label = h5("Do you have a separate strata file?"), 
                             choices = list("Yes" = 1, "No" = 2), selected = 2),
                 
                 uiOutput("ui_strata.upload"),
                 
                 uiOutput("stratif.sep"),
                 
                 selectInput("select.ploidy", label = h5("Indicate the ploidy of your data here"),
                             choices = list("1" = 1, "2" = 2), selected = 1),
                 
                 uiOutput("ui_1"),
                 
                 
                 actionButton("loadGtypes", label = "Load gtypes object")),
               ###End###
               
               ###Panel that displays the data from the gtypes tab###
               mainPanel(
                 tabsetPanel(type = "tabs", 
                             
                             tabPanel("Genetic Data", verbatimTextOutput('gen.data.head')), 
                             
                             tabPanel("FASTA", verbatimTextOutput('stats.fasta')),
                             
                             tabPanel("Strata", verbatimTextOutput('stats.strata')),
                             
                             tabPanel("gtypes object", 
                                      verbatimTextOutput('stats_gtype'),
                                      helpText("Click the button below to save your gtypes object as an RData file"),
                                      actionButton("save.gtypes", label = "Save"))
                             
                 )))),
    ###End###
    
    ###Sidebar of the qaqc tab###
    tabPanel("QA/QC",
             titlePanel("QA/QC"),
             helpText("Select the QA/QC function(s) you would like to run and click the Run 
                      button to perform the analyses and display the output(s) in the tabs below. 
                      You may save the outputs to your specified directory in .csv format by 
                      clicking the Save button"),
             
             sidebarLayout(
               sidebarPanel(h4("QA/QC Summary functions"), 
                            
                            checkboxInput("checkbox.by.sample", "by sample summary"),
                            checkboxInput("checkbox.by.locus", "by locus summary"),
                            checkboxInput("checkbox.dup.gen", "Duplicate genotypes"),
                            
                            uiOutput("ui_dup.gen"),
                            
                            uiOutput("ui_checkbox.seq.sum"),
                            uiOutput("ui_checkbox.hap.like"),
                            uiOutput("ui_checkbox.low.freq"),
                            uiOutput("ui_low.freq"),
                            
                            actionButton("run_qaqc", label="Run"),
                            actionButton("save.qaqc", label = "Save")
               ),
               ###End###
               
               ###Main panel for qaqc that displays the data###
               mainPanel(
                 tabsetPanel(type ="tabs",
                             tabPanel("by sample summary", verbatimTextOutput('stats.by.sample')),
                             
                             tabPanel("by locus summary", verbatimTextOutput('stats.by.locus')),
                             
                             tabPanel("Duplicate genotypes", verbatimTextOutput('stats.dup.gen')),
                             
                             tabPanel("Sequence summary",verbatimTextOutput('stats.seq.sum')),
                             
                             tabPanel("Haplotype likelihood",
                                      tabsetPanel(type="tabs",
                                                  tabPanel("Haplotype likelihood plot",
                                                           plotOutput("plot.hap.like")),
                                                  tabPanel("Haplotype likelihood summary",
                                                           verbatimTextOutput('stats.hap.like')))),
                             
                             tabPanel("Low frequency substitutions", verbatimTextOutput('stats.low.freq'))
                 )))),
    ###End###
    
    ###Sidebar panel for population structure tab###
    tabPanel("Population Structure",
             titlePanel("Population Structure"),
             helpText("This tab allows you perform population structure analyses. First choose to run the population structure
                      tests in an overall and/or pairwise manner by checking the box or boxes in the input panel. 
                      Then select the specific metric(s) to be computed and specify the number of permutations in the box below.
                      Click the Run button to analyze the data and display the outputs in the viewer tabs. You may save the 
                      outputs to your specified directory in .csv format by clicking the Save button."),
             
             sidebarLayout(
               sidebarPanel(h4("Types"),
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
                            
                            selectInput("pairwise.deletions", label = h5("Pairwise Deletion"),
                                        choices = list("Yes" = 1, "No" = 2), selected = 1),
                            numericInput("permutation", label = h5("Number of permutations"), 
                                         min = 0, value = 1000),
                            
                            actionButton("run_pop.struc", label = "Run"),
                            actionButton("save.pop.structure", label = "Save")),
               ###End###
               
               ###Main panel for the population structure tab that display pop structure data###
               mainPanel(
                 tabsetPanel(type="tabs",
                             tabPanel("Overall", 
                                      verbatimTextOutput("stats_overall.pop")),
                             tabPanel("Pairwise", 
                                      verbatimTextOutput("stats_pairwise.pop")))
               )
               ###End###
             )))
  
)

)