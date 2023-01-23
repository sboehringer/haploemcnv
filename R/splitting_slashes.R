#' @title Splitting on the slashes ("/")
#' 
#' @description Support function of \code{\link{all_options}}. Function can only correctly handle genotype calls with allelic ambiguities (`/`) and cannot correctly handle arguments with a `+` or genotypic ambiguities (`|`).
#' 
#' @param xx A character. The genotype call of an individual.
#' 
#' @return All possible genotypes compatible with the genotype call of the individual.
#' 
# #' @examples 
# #' \dontrun{splitting_slashes("001/002+003/004/005")}
#' 
splitting_slashes = function(xx){
  splttd = unlist(strsplit(xx, "\\+")); len_splttd = length(splttd)
  
  splitted = as.list(NULL)
  for(j in 1:len_splttd){
    splitted[[j]] = sort(unlist(strsplit(splttd[j], "\\/")))
  }
  
## Select all possible combinations and keep only the options which have no NA
  splitted.grid = expand.grid(as.data.frame(t(make_data_frame(splitted))))
  good.grid = as.data.frame(splitted.grid[which(complete.cases(splitted.grid)), ])
  len_good.grid = nrow(good.grid)
  
  mult_geno = NULL
  for(j in 1:len_good.grid){
    mult_geno = c(mult_geno, paste(sort(as.matrix(good.grid[j, ])), collapse = "+"))
  }
  
  return(mult_geno)
}
