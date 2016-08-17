###Upper sections of gtypes panel that gives overall directions and asks user to upload###
###already existing gtypes object###
    tabPanel(
      title = "gtypes",
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
           
      ###Side bar panel that allows user to create a gtypes object from uploaded data###
      sidebarLayout(
        sidebarPanel(
          textInput("save.directory", 
                     label = h5("Input the datapath for the directory where you would like to save any ouputted files. This will default
                                 to your Desktop if you do not specify a different datapath."), 
                     value = "~/Desktop"),
           textInput("description",
                     label = h5("Provide a description of the inputted data. This description will serve as",
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
                           
               ))))
  ###End###