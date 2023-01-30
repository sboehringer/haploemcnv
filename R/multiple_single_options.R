#' @title Handling the genotype results for multiple single gene copy genotypes
#' 
#' @description Support function of \code{\link{all_options}}. Function can only correctly handle single gene copies with allelic (`/`) or genotypic ambiguities (`|`) and cannot correctly handle arguments with a `+`.
#' 
#' @param xx A character. The possible genotype option of the individual. The genotype of the i-th donor, can not handle "+".
#' @param sign A character. On which sigh the genotype must be split, can be either "/" (allelic ambiguities) or "|" (genotypic ambiguities).
#' @param NEGs An optional character. How the `NEG` allele should be called, by default it will be called `NEG`.
#'
#' @return All possible genotype for the single gene copy genotype.
#'
# #' @examples
# #' multiple_single_options("001/002/003", sign = "/")
# #' multiple_single_options("004|005", sign = "|")
# #' \dontrun{multiple_single_options("004|005", sign = "/")}
#' 
multiple_single_options <- function(xx, sign, NEGs = "NEG"){
  splt_sign <- paste(c("\\", sign), collapse = "")
  splitted <- unlist(strsplit(xx, splt_sign))
  
  mult_opt <- NULL
  for(j in 1:length(splitted)){
    mult_opt[j] <- paste(c(splitted[j], NEGs), collapse = "+")
  }
  
  return(mult_opt)
}
