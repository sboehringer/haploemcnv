#' @title Reducing haplotype or diplotype list elements
#' 
#' @description To reduce the hapltoype or diplotype lists to a specified set of genes.
#' 
#' @param lst A list. Data containing the possible diplotypes of each donor as a separate list element, can be obtained via \code{\link{all_options}}.
#' @param nr_gene An integer or a vector. Which genes \code{lst} needs to be reduced to.
#' 
#' @return A list with the haplotypes or diplotypes of the specified genes only.
#'
reducing_list = function(lst, nr_gene){

  lst_out = lapply(lst, function(p){unlist(lapply(lapply(strsplit(p, "\\+"), function(x){lapply(strsplit(x, "\\-"), function(y){paste(y[nr_gene], collapse = "-")})}), function(z){paste(sort(unlist(z)), collapse = "+")}))})
  
  return(lst_out)
}
