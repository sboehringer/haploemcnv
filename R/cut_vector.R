#' @title Split a numeric vector in subvectors
#' 
#' @description To get all K-1 subvectors for a K vector. Support function for various functions.
#' 
#' @param vect An integer vector. Must range from 1 till K.
#' @param vect_list An optional list. Each list element contains a integer (sub)vector to which a new integer needs to be added. Can be obtained via \code{cut_vector}.
#'
#' @return A list with each list element contains the K-1 subvectors.
#' 
# #' @examples 
# #' vect_lst = cut_vector(vect = 1:3)
# #' 
# #' \dontrun{
# #' cut_vector(vect = 1:4, vect_list = vect_lst)
# #' }
#' 
cut_vector = function(vect, vect_list = NULL){
  
  nr_genes = length(vect)
  nr_genes_p1 = nr_genes + 1; nr_genes_m1 = nr_genes - 1
  
  vect_reduced = list(NULL); reduced_index = 0
  help_val = nr_genes_m1
  
  if(is.null(vect_list)){
    for(i in 1:help_val){
      while(help_val < nr_genes_p1){
        
        if(i != help_val){
          reduced_index = reduced_index + 1
          vect_reduced[[reduced_index]] = c(i, help_val)
        }
        
        help_val = help_val + 1
      }
      help_val = nr_genes_m1
    }
    
    
  } else {
    len_vect_list = length(vect_list)
    for(i in 1:len_vect_list){
      while(help_val < nr_genes_p1){
        
        if(vect_list[[i]][length(vect_list[[i]])] < help_val){  # vect_list[[i]][nr_genes - 2] < help_val){
          reduced_index = reduced_index + 1
          vect_reduced[[reduced_index]] = c(vect_list[[i]], help_val) 
        }
        
        help_val = help_val + 1
      }
      help_val = nr_genes_m1
    }
    
  }
  
  return(vect_reduced)
}
