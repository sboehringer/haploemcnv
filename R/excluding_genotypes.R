#' @title Excluding specific diplotypes
#' 
#' @description To exclude specific diplotypes from the diplotype lists of individuals. Either by fully specifying the diplotypes to discard, or by specifying a allele that, if present in a diplotype, will discard all diplotypes. 
#' 
#' @param lst A list. All possible diplotypes of a donor, as obtained via \code{\link{all_options}}
#' @param remove_alleles An optional vector. If a diplotype contains any of these allele(s), they will be discarded.
#' @param remove_genotypes An optional vector. Which diplotypes need to be discarded.
#' @param no_hat3plus A logical scalar. Whether or not all diplotypes that have a haplotype that consists of three or more genes need to be discarded (\code{FALSE} is default). 
#' 
#' @return A list with the remaining possible diplotypes of a donor. There is also a vector with all removed diplotypes and a vector with the donors that have no diplotypes left.
#' 
# #' @examples 
# #' lst = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' \dontrun{
# #' excluding_genotypes(lst, remove_alleles = "002")
# #' }
#' 
excluding_genotypes = function(lst, remove_alleles = NULL, remove_genotypes = NULL, no_hat3plus = FALSE){
  
  rm_geno = NULL
  if(!is.null(remove_alleles)){
    
    un_geno = unique(unlist(lst))
    len_un_geno = length(un_geno)
    
    splt = strsplit(un_geno, "\\+")
    
    len_remove_alleles = length(remove_alleles)
    for(i in 1:len_remove_alleles){
      present = unlist(lapply(splt, function(x){remove_alleles[i] %in% x}))
      rm_geno = c(rm_geno, un_geno[present])
    }
  }
  
  
  if(no_hat3plus == TRUE){
    un_genos = unique(unlist(lst))
    
    splt = strsplit(un_genos, "\\+")
    len_splt = length(splt)
    
    index = NULL
    for(i in 1:len_splt){
      lens_splt = lengths(strsplit(splt[[i]], "\\^"))
      
      if(TRUE %in% (lens_splt > 2)){
        index = c(index, i)
      }
    }
    
    rm_geno = c(rm_geno, un_genos[index])
  }
  
  
  all_genos_remove = unique(c(rm_geno, remove_genotypes))
  
  empty_donors = NULL
  
  len_lst = length(lst)
  for(i in 1:len_lst){
    xx = lst[[i]]
    
    if(sum(all_genos_remove %in% xx) != 0){
      index = which(xx %in% all_genos_remove)
      
      lst[[i]] = lst[[i]][-index]
      
      if(is_empty(lst[[i]])){
        empty_donors = c(empty_donors, i)
      }
    }
  }
  
  return(list("lst" = lst, "removed_genotypes" = all_genos_remove, "empty_donors" = empty_donors))
}
