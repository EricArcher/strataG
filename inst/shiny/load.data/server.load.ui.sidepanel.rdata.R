uiLoadRdata <- function() {    
  fnames <- if(is.null(vals$wd)) NULL else {
    dir(vals$wd, pattern = "(*.rdata$)|(*.rda$)", ignore.case = TRUE)
  }
  verticalLayout(
    selectInput(
      "slctRdataFile",
      label = h5("Select an .rdata file"),
      choices = fnames
    ),
    selectInput(
      "slctGtypes",
      label = h5("Select a gtypes object"),
      choices = NULL
    )
  )
}