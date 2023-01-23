#' @title Calculates frequencies for the naive analysis
#' 
#' @description To calculate haplotype and diplotype frequencies for the naive analysis. 
#'
#' @param lst A list. Data containing the possible diplotypes of each donor as a separate list element, can be obtained via \code{\link{all_options}}.
#' @param haplo A logical scalar. Whether or not diplotype frequencies are estimated by first estimating haplotype frequenices, or immediately start by estimating diplotype frequencies (\code{TRUE} is default). 
#' @param adj_comb_vals An optional vector. \lars{Currently unknown}.
#' @param obs_diplos A logical scalar. Whether or not diplotype frequencies are estimated for only the observed diplotypes, or also for all theoretically possible diplotypes (\code{TRUE} is default).
#'
#' @return 
#' 
# #' @examples
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' \dontrun{
# #' frequencies_naive(lst = x)
# #' }
#' 
frequencies_naive = function(lst, haplo = TRUE, adj_comb_vals = NULL, obs_diplos = TRUE){
  
  len_lst = length(lst)
  
  if(!haplo & !obs_diplos){
    obs_diplos = TRUE
    warning("Diplotype frequencies are only estimated for the actual observed diplotypes")
  }
  
  if(obs_diplos){
    all_diplos = unique(unlist(lst)); len_all_diplos = length(all_diplos)
    diplo_freq_vect = rep(0, len_all_diplos); names(diplo_freq_vect) = all_diplos
  }
  
  if(haplo){
    all_haplos = unique(unlist(strsplit(unique(unlist(lst)), "\\+"))); len_all_haplos = length(all_haplos)
    haplo_freq_vect = rep(0, len_all_haplos); names(haplo_freq_vect) = all_haplos
    
    if(!obs_diplos){
      all_diplos = unique(unlist(apply(expand.grid(all_haplos, all_haplos), 1, function(x){paste(sort(x), collapse = "+")}))); len_all_diplos = length(all_diplos)
      diplo_freq_vect = rep(0, len_all_diplos); names(diplo_freq_vect) = all_diplos
    }
    
    for(i in 1:len_lst){
      xx = unlist(strsplit(lst[[i]], "\\+")); len_xx = length(xx); add_value = 1 / len_xx
      
      for(j in 1:len_xx){
        index = xx[j]
        
        haplo_freq_vect[index] = haplo_freq_vect[index] + add_value
      }
    }
    haplo_freq = haplo_freq_vect / sum(haplo_freq_vect)
    
    if(!is.null(adj_comb_vals)){
      comb_names = names(adj_comb_vals)
      len_comb_names = length(comb_names)
      
      present_comb = comb_names[sapply(comb_names, function(x){TRUE %in% str_detect(all_haplos, x)})]
      present_comb = present_comb[length(present_comb)]
      
      haplo_index = which(str_detect(all_haplos, present_comb))
      
      nr_genes = length(unlist(strsplit(str_sub(present_comb, start = 6), "")))
      haplo_freq[haplo_index] = (adj_comb_vals[present_comb]^(1 / nr_genes))^nr_genes
      

      haplo_freq[-haplo_index] = haplo_freq[-haplo_index] / sum(haplo_freq[-haplo_index]) * (1 - haplo_freq[haplo_index])
    }

        
    diplo_freq_vect = HF_2_DF(Diplotypes = all_diplos, haplo_freq = haplo_freq)
#    for(i in 1:len_all_diplos){
#      diplo_splt = unlist(strsplit(all_diplos[i], "\\+"))
#      
#      diplo_splt_one = diplo_splt[1]
#      diplo_splt_two = diplo_splt[2]
#      
#      if(diplo_splt_one == diplo_splt_two){
#        diplo_freq_vect[i] = haplo_freq[diplo_splt_one] * haplo_freq[diplo_splt_two]
#      } else {
#        diplo_freq_vect[i] = 2 * haplo_freq[diplo_splt_one] * haplo_freq[diplo_splt_two]
#      }
#    }
    
  } else {
    
    for(i in 1:len_lst){
      xx = lst[[i]]
      len_xx = length(xx)
      
      add_value = 1 / len_xx
      
      for(j in 1:len_xx){
        index = xx[j]
        diplo_freq_vect[index] = diplo_freq_vect[index] + add_value
      }
    }
  }
  
  diplo_freq = diplo_freq_vect / sum(diplo_freq_vect)


  if(haplo == TRUE){
    return(list("Diplo_freq" = sort(diplo_freq, decreasing = TRUE), "Haplo_freq" = sort(haplo_freq, decreasing = TRUE)))
  } else {
    return(list("Diplo_freq" = sort(diplo_freq, decreasing = TRUE)))
  }
}
