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
                        uiOutput("ui_checkbox.seq.like"),
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
                         
                         tabPanel("Sequence likelihood",
                                  tabsetPanel(type="tabs",
                                              tabPanel("Sequence likelihood plot",
                                                       plotOutput("plot.seq.like")),
                                              tabPanel("Sequence likelihood summary",
                                                       verbatimTextOutput('stats.seq.like')))),
                         
                         tabPanel("Low frequency substitutions", verbatimTextOutput('stats.low.freq'))
             ))))
###End###