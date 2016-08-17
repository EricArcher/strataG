shinyServer(function(input, output, session) {
  source("server.gtypes.R", local = TRUE)
  source("server.qaqc.R", local = TRUE)
  source("server.popstruct.R", local = TRUE)
})