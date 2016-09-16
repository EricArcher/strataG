VolumeRoots <- c(HOME = Sys.getenv("HOME"), getVolumes()())
shinyDirChoose(input, "workingDirBtn", session = session, roots = VolumeRoots)
wd <- reactive({
  if(!is.null(input$workingDirBtn)) {
    vals$wd <- parseDirPath(VolumeRoots, input$workingDirBtn)
    vals$wd
  } else NULL
})
output$workingDir <- renderText(wd())