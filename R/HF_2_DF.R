#' @title Calculate frequencies for all diplotypes
#' 
#' @description Calculation of diplotype frequencies based on haplotype frequencies. 
#' 
#' @param Diplotypes A vector. Diplotype names for which the frequency needs to be estimated.
#' @param haplo_freq A vector. Haplotypes and their frequency.
#'
#' @return A vector with diplotype frequencies.
#' 
HF_2_DF = function(Diplotypes, haplo_freq){
  
  len_Diplotypes = length(Diplotypes)
  
  diplo_freq_vect <- rep(0, len_Diplotypes)
  names(diplo_freq_vect) <- Diplotypes
  
  for(i in 1:len_Diplotypes){
    diplo_splt <- unlist(strsplit(Diplotypes[i], "\\+"))
    
    diplo_splt_one <- diplo_splt[1]
    diplo_splt_two <- diplo_splt[2]
    
    if(diplo_splt_one == diplo_splt_two){
      diplo_freq_vect[i] <- haplo_freq[diplo_splt_one] * haplo_freq[diplo_splt_two]
    } else {
      diplo_freq_vect[i] <- 2 * haplo_freq[diplo_splt_one] * haplo_freq[diplo_splt_two]
      
    }
  }
  
  return(diplo_freq_vect)
}
