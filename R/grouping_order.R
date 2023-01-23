#' @title Determines which gene is the best choice to add into the profile EM-algorithm reconstruction
#' 
#' @description Estimates the mean squared difference between the allele frequencies of the starting genes with and without the addition of the the new genes. The gene that is added into the reconstruction maximizes the mean squared difference. Support function of \code{grouping_order_all}. 
#' 
#' @param freqs1 A vector. Allele frequencies of the analysis of the starting gene without other genes.
#' @param freqs2 A vector. The allele frequencies of the starting gene in a grouping configuration.
#' @param comb_group An optional vector. Containing the names of the alleles grouped in the compound haplotype. If not \code{NULL} the frequency of the compound haplotype is equally distributed over all alleles. 
#' @param freqs1_no_comb A logical scalar. Whether or not the compound haplotype frequency of \code{freqs2} is equally distributed over all alleles that are not grouped in \code{freqs1}, but are grouped in \code{freqs2} (\code{FALSE} is default).
#'
#' @return A list with the mean difference between the alleles and the names of the alleles which are used for calculating the difference.
#' 
# #' @examples 
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' xy = combining_genes(x, y)
# #' 
# #' out_x = EM_algorithm(x)
# #' out_xy = EM_algorithm(xy)
# #' \dontrun{
# #' grouping_order(freqs1 = out_x$Freq[nrow(out_x$Freq), ], freqs2 = Haplo2AF(out_xy)[[1]])
# #' }
#' 
grouping_order = function(freqs1, freqs2, comb_group = NULL, freqs1_no_comb = FALSE){
  names_freqs1 = names(freqs1); names_freqs2 = names(freqs2)
  
  if(!is.null(comb_group)){
    len_comb_group = length(comb_group)
    
    names_subset1 = c(names_freqs1, comb_group); names_subset1 = names_subset1[!(str_detect(names_subset1, "Comb"))]
    names_subset2 = names_subset1
    len_names_subset = length(names_subset1)
    
    subset1 = vector("numeric", length = len_names_subset)
    names(subset1) = names_subset1
    
    comb_index1 = which(str_detect(names_freqs1, "Comb"))
    comb_value1 = sum(freqs1[comb_index1]) / len_comb_group
    
    
    subset2 = vector("numeric", length = len_names_subset)
    names(subset2) = names_subset2
    
    comb_index2 = which(str_detect(names_freqs2, "Comb"))
    
    comb_value2 = sum(freqs2[comb_index2]) / len_comb_group
    
    for(i in 1:len_names_subset){
      name = names_subset1[i]
      
      if(name %in% names_freqs1){
        index = which(names_freqs1 == name)
        subset1[i] = freqs1[index]
      } else {
        subset1[i] = comb_value1
      }
      
      if(name %in% names_freqs2){
        index = which(names_freqs2 == name)
        subset2[i] = freqs2[index]
      } else {
        subset2[i] = comb_value2
      }
    }
    
  } else {
    if(freqs1_no_comb == TRUE){
      len_names_subset = length(names_freqs1)
      subset1 = freqs1
      
      names_subset1 = names_freqs1
      names_subset2 = names_freqs1
      
      subset2 = vector("numeric", length = len_names_subset)
      names(subset2) = names_subset2
      
      for(i in 1:len_names_subset){
        name = names_freqs1[i]
        
        if(name %in% names_freqs2){
          subset2[name] = freqs2[name]
        }
      }
      
      subset2_combs = which(subset2 == 0)
      len_subset2_combs = length(subset2_combs)
      
      comb_index = which(str_detect(names_freqs2, "Comb"))
      
      comb_value = sum(freqs2[comb_index]) / length(subset2_combs)
      subset2[subset2_combs] = comb_value
      
    } else {
      index = which(names_freqs1 %in% names_freqs2)
      
      subset1 = freqs1[index]
      names_subset1 = names(subset1)
      len_names_subset = length(names_subset1)
      
      subset2 = freqs2
      names_subset2 = names(subset2)
    }
  }
  
  
  diff = NULL
  for(i in 1:len_names_subset){
    name = names(subset1[i])
    value1 = subset1[i]
    value2 = subset2[which(names_subset2 == name)]
    
    diff[i] = (value1 - value2) ** 2
  }
  
  if(!is.null(comb_group)){
    out = list("Diff" = mean(diff), "Remain" = names_freqs1, "From_comb" = comb_group)
  } else {
    if(freqs1_no_comb == TRUE){
      out = list("Diff" = mean(diff), "Remain" = names_subset1)
    } else {
      out = list("Diff" = mean(diff), "Remain" = names_subset1) 
    }
  }
  
  return(out)
}
