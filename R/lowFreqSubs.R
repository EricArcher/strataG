#' @title Low Frequency Substitutions
#' @description Check nucleotide sites for low frequency substitutions.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' @param min.freq minimum frequency of base to be flagged.
#' @param motif.length length of motif around low frequency base to output.
#' @param simplify if there is a single locus, return result in a simplified
#'   form? If \code{FALSE} a list will be returned wth one element per locus.
#' 
#' @return data.frame listing id, site number, and motif around low frequency 
#'   base call.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.haps)
#' 
#' lowFreqSubs(dolph.haps)
#' 
#' @export
#' 
lowFreqSubs <- function(x, min.freq = 3, motif.length = 10, simplify = TRUE) {  
  result <- sapply(
    apex::getSequences(as.multidna(x), simplify = FALSE), 
    function(dna) {
      motif.half <- max(1, round(motif.length / 2, 0))
      var.sites <- variableSites(dna)$site.freq
      has.min.freq <- apply(var.sites, 2, function(freq) {
        sum(freq > 0 & freq < min.freq) > 0
      })
      if(!any(has.min.freq)) return(NULL)
      sites.w.min.freq <- var.sites[, has.min.freq, drop = FALSE]
      dna <- as.character(as.matrix(dna))
      sites.to.check <- purrr::map(colnames(sites.w.min.freq), function(col) {
        position <- as.numeric(col)
        site.freq <- sites.w.min.freq[, col]
        bases <- names(site.freq)[which(site.freq > 0 & site.freq < min.freq)]
        tibble::tibble(
          id = names(dna[, position])[dna[, position] %in% bases], 
          site = position
        ) %>% 
          dplyr::mutate(
            base = dna[.data$id, position],
            freq = sapply(.data$base, function(i) var.sites[i, col]),
            motif = sapply(.data$id, function(i) {
              start.bp <- max(1, position - motif.half)
              end.bp <- min(ncol(dna), position + motif.half)
              paste(dna[i, start.bp:end.bp], collapse = "")
            })
          )
      }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::arrange(.data$id, .data$site) %>% 
        as.data.frame()
      rownames(sites.to.check) <- NULL
      sites.to.check
    },
    simplify = FALSE
  )
  
  if(length(result) == 1 & simplify) result[[1]] else result
}