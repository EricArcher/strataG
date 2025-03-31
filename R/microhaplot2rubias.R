#' @title Convert MicroHaplot file to rubias format
#' @description Converts a reference or mixture file from the MicroHaplot 
#'   package to the proper format for the rubias package
#'   
#' @param x data frame from MicroHaplot or the name a of .csv file.
#' @param sample.type determines type of MicroHaplot sample. Can be any 
#'   unambiguous substring of \code{"reference"} or \code{"mixture"}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} Converted from code by
#'   John Horne (\url{https://johnbhorne.wordpress.com/2019/03/22/snp-haplotypes-for-mixed-stock-fishery-analysis-microhaplot-to-rubias-conversion-in-r/})
#'
#' @export
#' 
microhaplot2rubias <- function(x, sample.type = c("reference", "mixture")) {
  x <- if(is.character(x)) {
    if(!file.exists(x)) stop("The file '", x, "' cannot be found.")
    utils::read.csv(x)
  } else if(is.data.frame(x)) x else {
    stop("'x' is not a character or data.frame")
  }
  
  rubias <- split(x, list(x$group, x$indiv.ID)) |> 
    lapply(function(id.df) {
      id.df <- id.df[order(id.df$locus), ]
      
      # count number of heterozygotes per locus (duplicated locus entries)
      hets <- duplicated(id.df$locus)
      hets[which(hets) - 1] <- TRUE
      rep.times <- as.numeric(hets)   
      # loci w/ 0 duplicates are homozygous, still need 2 alleles
      rep.times[rep.times == 0] <- 2
      
      # create new data frame
      id.df <- id.df[rep(1:nrow(id.df), rep.times), ]
      i <- rep(1:2, length.out = nrow(id.df))
      id.df$locus <- paste(id.df$locus, i, sep = ".")
      id.df |> 
        tidyr::pivot_wider(
          id_cols = c("group", "indiv.ID"), 
          names_from = "locus",
          values_from = "haplo"
        )
    }) |> 
    dplyr::bind_rows() |> 
    dplyr::rename(
      collection = "group",
      indiv = "indiv.ID"
    ) |> 
    as.data.frame()
  
  sample.type <- match.arg(sample.type)
  cbind(
    data.frame(
      sample_type = rep(sample.type, nrow(rubias)), 
      repunit = if(sample.type == "reference") {
        rubias$collection
      } else {
        rep(NA, nrow(rubias))
      },        
      stringsAsFactors = FALSE
    ),
    rubias
  )
}
