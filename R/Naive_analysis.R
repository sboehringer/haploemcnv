#' @title Naive HT reconstruction analysis
#' 
#' @description Runs two naive analyses to be compared to the analysis with the EM-algorithm based \code{\link{HT_reconstruction}}.
#' 
#' @param lst A list. Data containing the possible diplotypes of each donor as a separate list element, can be obtained via \code{\link{all_options}}.
#' @param haplo A logical scalar. Whether or not diplotype frequencies are estimated by first estimating haplotype frequencies, or immediately start by estimating diplotype frequencies (\code{TRUE} is default).
#' @param CO_thresh A numeric value \in \{0, 1\}. The cut-off threshold, only haplotypes lower than this threshold are be collapsed. 
#' 
#' @return A choice vector for the dictatorial naive analysis and a probability matrix for the democratic naive analysis. For both methods also haplotype frequencies are returned if \code{haplo} is \code{TRUE}
#' 
# #' @examples 
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' \dontrun{
# #' Naive_analysis(lst = list(x, y))
# #' }
#' 
Naive_analysis = function(lst, haplo = TRUE, CO_thresh = 1e-7){
  
  nr_genes = length(lst)
  
  naive_dictatorial = function(x, HTFs = NULL, DTFs){
    y = DTFs[names(DTFs) %in% unlist(x)]; len_y = length(y)
    
    if(len_y == 0){
      if(!is.null(HTFs)){
        z = sample(x, size = 1)
        
      } else {
        xx = unlist(lapply(strsplit(x, "[\\+^]"), function(x){prod(HTFs[x])}))
        z = x[sample(which(xx == max(xx)), size = 1)]
      }
      
      
    } else if(len_y == 1){
      z = names(y)
      
    } else {
      ymax = which(y == max(y)); len_ymax = length(ymax)
      if(len_ymax == 1){
        z = names(y[ymax])
      } else {
        z = names(sample(y, size = 1))
      }
    }
    
    return(z)
  }
  
  # naive_democratic = function(x, DTFs, CO_thresh){
  #   comb_x = combining_AF(lst = x, haplotypes = HTFs, CO_thresh = CO_thresh)
  #   
  #   x = comb_x$lst
  #   y = DTFs[names(DTFs) %in% unlist(x)]; len_y = length(y)
  # 
  #   return(list(x = x, y = y))
  # }
  
  
  
  
  ## Dictatorial way
  dict = rep(list(NULL), nr_genes)
  for(i in 1:nr_genes){
    lst_i = lst[[i]]; unamb_index = (lengths(lst_i) == 1)
    
    lst_i_freqs = frequencies_naive(lst = lst_i[unamb_index], haplo = haplo, obs_diplos = FALSE)
    HTFs = lst_i_freqs$Haplo_freq; DTFs = lst_i_freqs$Diplo_freq
    
    dict[[i]] = lapply(lst_i, function(x){naive_dictatorial(x, HTFs = HTFs, DTFs = DTFs)})
  }

  dict_lst = dict[[1]]
  for(i in 2:nr_genes){
    dict_lst = combining_genes(dict_lst, dict[[i]])
  }

  dict_lst_DTFs = frequencies_naive(lst = dict_lst, haplo = haplo, obs_diplos = TRUE)$Diplo_freq
  dict_choices_lst = lapply(dict_lst, naive_dictatorial, NULL, dict_lst_DTFs)
  
  dict_HTFs = frequencies_naive(lst = dict_choices_lst, haplo = TRUE)$Haplo_freq
  
  

## Democratic way
  # demo = rep(list(NULL), nr_genes)
  # for(i in 1:nr_genes){
  #   lst_i = lst[[i]]
  #   
  #   lst_i_freqs = freqs[[i]] = frequencies_naive(lst = lst_i, haplo = haplo, obs_diplos = TRUE)
  # 
  #   demo[[i]] = lapply(lst_i, function(x){naive_democratic(x, DTFs = lst_i_freqs$Diplo_freq, CO_thresh)})
  # }
  # 
  # demo_lst = demo[[1]]
  # for(i in 2:nr_genes){
  #   curr = lapply(demo_lst, function(y){y$y})
  #   new = lapply(demo[[i]], function(y){y$y})
  #   
  #   demo_lst = combining_genes(demo_lst, lapply(demo[[i]], function(y){y$y}))
  # }
  
  # demo_choices_lst = lapply(dict_lst, naive_democratic, dict_lst_DTFs)
  
  demo_init_mat = EM_init_mat(dict_lst)
  dict_lst_DTFs = dict_lst_DTFs[colnames(demo_init_mat)]
  
  demo_mat = t(apply(demo_init_mat, 1, function(x){x * dict_lst_DTFs})); demo_mat = demo_mat / rowSums(demo_mat)
  
  diplos = colnames(demo_mat); D = length(diplos); diplos_splt = strsplit(diplos, "\\+")
  haplos = unique(unlist(diplos_splt)); M = length(haplos)
  
  Repeat_mat = matrix(0, nrow = D, ncol = M, dimnames = list(diplos, haplos))
  for(j in 1:M){
    Repeat_mat[, j] = unlist(lapply(diplos_splt, function(x){sum(x == haplos[j])}))
  }
  
  demo_HTFs = colSums(demo_mat %*% Repeat_mat); demo_HTFs = demo_HTFs / sum(demo_HTFs); demo_HTFs = sort(demo_HTFs, decreasing = TRUE)
  
  
  return(list("dictatorial_data" = dict_choices_lst, "Democratic_mat" = demo_mat, "dictatorial_freqs" = dict_HTFs, "Democratic_freqs" = demo_HTFs))
}
