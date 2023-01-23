#' @title Handling all possible gene copy genotype of a genotype with at least one "+", including the multi-allele haplotypes
#' 
#' @description Support function of \code{\link{all_options}}. Function can only correctly handle an argument that has a `+` and cannot correctly handle arguments that contain allelic (`/`) or genotypic ambiguities (`|`).
#'
#' @param xx A character. The possible genotype option of the individual.
#'
#' @return All possible genotypes with the multi-allele haplotype
#' 
# #' @examples
# #' \dontrun{
# #' multiple_gene_copies("001+002")
# #' multiple_gene_copies("001+002+003")
# #' }
#' 
multiple_gene_copies <- function(xx, NEGs = NULL){
  splitted <- sort(unlist(strsplit(xx, "\\+")))
  len_splitted <- length(splitted)
  
  out <- NULL
  if(len_splitted == 1){
    out <- paste(c(xx, NEGs), collapse = "+")
    
  } else {
    
    if(len_splitted == 2){
      out <- c(xx, paste(sort(c(paste(splitted[1], splitted[2], sep = "^"), NEGs)), collapse = "+"))  # A+B or AB+NEG
      
    } else {
      
      if(len_splitted == 3){
        out <- c(paste(sort(c(paste(splitted[1], splitted[2], sep = "^"), splitted[3])), collapse = "+"),  # AB+C
                 paste(sort(c(paste(splitted[1], splitted[3], sep = "^"), splitted[2])), collapse = "+"),  # AC+B
                 paste(sort(c(paste(splitted[2], splitted[3], sep = "^"), splitted[1])), collapse = "+"),  # BC+A
                 paste(sort(c(paste(splitted[1], splitted[2], splitted[3], sep = "^"), NEGs)), collapse = "+"))  # ABC+NEG
        
      } else {
        
        if(len_splitted == 4){
          out <- c(paste(sort(c(paste(splitted[1], splitted[2], sep = "^"), paste(splitted[3], splitted[4], sep = "^"))), collapse = "+"),  # AB+CD
                   paste(sort(c(paste(splitted[1], splitted[3], sep = "^"), paste(splitted[2], splitted[4], sep = "^"))), collapse = "+"),  # AC+BD
                   paste(sort(c(paste(splitted[1], splitted[4], sep = "^"), paste(splitted[2], splitted[3], sep = "^"))), collapse = "+"),  # AD+BC
                   paste(sort(c(paste(splitted[1], splitted[2], splitted[3], sep = "^"), splitted[4])), collapse = "+"),  # ABC+D
                   paste(sort(c(paste(splitted[1], splitted[2], splitted[4], sep = "^"), splitted[3])), collapse = "+"),  # ABD+C
                   paste(sort(c(paste(splitted[1], splitted[3], splitted[4], sep = "^"), splitted[2])), collapse = "+"),  # ACD+B
                   paste(sort(c(paste(splitted[2], splitted[3], splitted[4], sep = "^"), splitted[1])), collapse = "+"),  # BCD+A
                   paste(sort(c(paste(splitted[1], splitted[2], splitted[3], splitted[4], sep = "^"), NEGs)), collapse = "+"))  # ABCD+NEG
        }
      }
    }
  }
  
  return(out)
}
