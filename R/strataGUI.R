#' @title strataG GUI
#' @description A graphical user interface for creating gtypes objects, conducting 
#'   quality control analyses, and analyses of population structure.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom shiny runApp
#' @importFrom shinyFiles shinyDirChoose
#' @importFrom DT datatable
#' @export
#' 
strataGUI <- function() {
  appDir <- system.file("shiny", package = "strataG")
  if(appDir == "") {
    stop("shiny app folder could not be found. Try re-installing 'strataG'", call. = FALSE)
  }
  runApp(appDir, display.mode = "normal")
}