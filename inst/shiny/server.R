shinyServer(
  function(input, output, session) {
    library(strataG)
    current.g <- msats.g
    source("server.gtypes.R", local = TRUE)
    source("server.qaqc.R", local = TRUE)
  }
)