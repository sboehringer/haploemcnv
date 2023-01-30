#' @title Make a table from multiple vectors
#' 
#' @description To visualize the effect of different methods on the estimates, the estimates need to be compared against each other. With \code{table_maker} sorts these vectors and pastes them next to each other. 
#' 
#' @param ... The vectors to be put next to each other based on their element names.
#' @param freqs_lst An optional list. Each list element contains a separate vector. If \code{NULL}, than vectors need to be supplied via \code{...}.
#' @param freq_order An optional vector. The names and order of the haplotypes to be plotted. If \code{NULL} all unique haplotypes will be plotted from high to low for the first vector.
#' @param decreasing An optional logical scalar. Whether or not the frequencies should be increasing or decreasing (default = \code{TRUE}). Only relevant when freq_order = \code{NULL}.
#' @param cnames An optional vector. Names of each of the vectors, must be supplied in the order in which the vectors are supplied.
#' @param total A logical scalar. Whether or not the last row contains a overview of how many elements each vector had.
#' @param rounding An integer. To how many digits all values need to be rounded. 
#' 
#' @return A table with each vector as a separate vector and each row a separate haplotype.
#' 
table_maker <- function(..., freqs_lst = NULL, freq_order = NULL, decreasing = TRUE, cnames = NULL, total = FALSE, rounding = 5){

  if(is.null(freqs_lst)){
    freqs_lst <- list(...)
  }
  
  if(is.null(cnames)){
    cnames = names(freqs_lst)
  }
  len_cnames <- length(freqs_lst)
  
  if(is.null(freq_order)){
    freqs_lst_unlist = NULL
    for(i in 1:len_cnames){
      freqs_lst_unlist = c(freqs_lst_unlist, freqs_lst[[i]]) 
    }
    freqs_lst_unlist = sort(freqs_lst_unlist, decreasing = TRUE)
    freqs_lst_unlist = freqs_lst_unlist[!duplicated(names(freqs_lst_unlist))]
    
    freq_order = names(sort(freqs_lst_unlist, decreasing = decreasing))
    len_freq_order = length(freq_order)
  } else {
    len_freq_order = length(freq_order)
  }

  tab <- matrix(NA, nrow = len_freq_order, ncol = len_cnames, dimnames = list(freq_order, cnames))
  
  for(i in 1:len_cnames){
    for(j in 1:len_freq_order){
      haplo_name <- freq_order[j]
      
      tab[j, i] <- freqs_lst[[i]][haplo_name]
    }
  }
  
  if(!is.null(rounding)){
    tab = round(tab, rounding)  
  }
  
  if(total){
    tab = rbind(tab, "#H" = lengths(freqs_lst))
  }
  
  return(tab)
}
