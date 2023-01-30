#' @title Naive HT reconstruction analysis with the EM-algorithm updated diplotype list
#' 
#' @description Runs the naive analysis with the EM-algorithm updated diplotype list, to be compared to the analysis of the complete EM-algorithm based \code{\link{HT_reconstruction}}.
#' 
#' @param lst A list. Data containing the possible diplotypes of each donor as a separate list element, can be obtained via \code{\link{all_options}}.
#' @param haplo A logical scalar. Whether or not diplotype frequencies are estimated by first estimating haplotype frequencies, or immediately start by estimating diplotype frequencies (\code{TRUE} is default). 
# #' @param lst_freqs_pre An optional vector. \lars{Currently unknown}.
# #' @param adj_comb_vals An optional vector. \lars{Currently unknown}.
#' 
Naive_EMbased_analysis = function(lst, haplo = TRUE){  # , lst_freqs_pre = NULL, adj_comb_vals = NULL){
  
  len_lst = length(lst)
  
  # if(is.null(lst_freqs_pre)){
    lst_freqs_pre = frequencies_naive(lst, haplo = haplo)
  # } 
  
  naive_choice = NULL
  for(i in 1:len_lst){
    xx <- lst[[i]]
    len_xx <- length(xx)
    
    if(len_xx == 1){
      naive_choice[i] <- xx
    } else {
      
      freqs_pres <- lst_freqs_pre$Diplo_freq[xx]
      
      max_freqs <- which(freqs_pres == max(freqs_pres))
      len_max_freqs <- length(max_freqs)
      
      if(len_max_freqs == 1){
        naive_choice[i] <- names(max_freqs)
      } else {
        naive_choice[i] <- names(sample(max_freqs, 1))
      }
    }
  }
  
  lst_freqs_post = frequencies_naive(as.list(naive_choice), haplo = haplo)
  
  
  if(haplo){
    out = list("Choice_vector" = naive_choice, "HFs_pre" = lst_freqs_pre$Haplo_freq, "HFs_post" = lst_freqs_post$Haplo_freq)
  } else {
    out = list("Choice_vector" = naive_choice, "DFs_pre" = lst_freqs_pre, "DFs_post", lst_freqs_post) 
  }
  
  return(out)
}
