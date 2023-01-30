#' @title Produces an initial value matrix 
#' 
#' @description Initial value matrix consists a 1 for all diplotypes that are compatible with the data, and a 0 for all diplotypes that are not compatible. Support function of \code{\link{EM_algorithm}}.
#'
#' @param lst A list. Each possible diplotype of an individual has a separate argument, can be obtained via \code{\link{all_options}}.
#' @param init_equal A logical scalar. Whether or not the probabilities in the initial value matrix are equal (a 0 when diplotype is not compatible with individuals diplotype call, and a 1 when the diplotype is compatible) or that the number of times a compatible diplotype is stated in \code{lst} equals the value the initial value matrix should have (\code{TRUE} is default).
#' @param start_HTFs An optional vector. Haplotype frequencies where the diplotype frequencies of the first iteration are started with. If a haplotype is not specified here, it will get a frequency equal to half of the minimum supplied frequency
#' 
#' @return A matrix with the initial values for the EM-algorithm. Each column represents a possible diplotype (that has been observed for at least one of the individuals) and each row represents a individual.
#' 
# #' @examples 
# #' lst = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' \dontrun{
# #' EM_init_mat(lst)
# #' }
#' 
EM_init_mat = function(lst, init_equal = TRUE, start_HTFs = NULL){
  
  N = length(lst)
  DTs = sort(unique(unlist(lst))); len_DTs = length(DTs)
  
  ##
  if(!is.null(start_HTFs)){
    HTs = names(start_HTFs); len_HTs = length(HTs)
    min_val = min(start_HTFs) / 2
    
    DTs_splt = strsplit(DTs, "\\+")
    start_DTFs = matrix(0, nrow = 1, ncol = len_DTs, dimnames = list(NULL, DTs))[1, ]
    for(i in 1:len_DTs){
      cname = DTs_splt[[i]]
      
      ##
      if(is.na(start_HTFs[cname[1]]) | is.na(start_HTFs[cname[2]])){
        cname[which(is.na(start_HTFs["TEST"]), is.na(start_HTFs[cname[2]]))] = min_val
      }
      ##

      if(cname[1] == cname[2]){
        start_DTFs[i] = start_HTFs[cname[1]] * start_HTFs[cname[2]]
      } else {
        start_DTFs[i] = start_HTFs[cname[1]] * start_HTFs[cname[2]] * 2
      }
    }
      
  } else {
    start_DTFs = matrix(1, nrow = 1, ncol = len_DTs, dimnames = list(NULL, DTs))[1, ]
  }
  ##
  
  matval = matrix(0, nrow = N, ncol = len_DTs, dimnames = list(NULL, DTs))
  for(i in 1:N){
    geno_donor = lst[[i]]
    len_geno_donor = length(geno_donor)
    
    for(j in 1:len_geno_donor){
      cname = geno_donor[j]
      index = which(DTs == cname)
      
      if(init_equal == TRUE){
        matval[i, index] = 1
      } else {
        matval[i, index] = matval[i, index] + 1
      }
    }
  }
  
  matval = t(apply(matval, 1, function(x){x * start_DTFs}))
  
  
  
  # if(!is.null(external_info)){
  #   if(min_val == TRUE){
  #     min_value = min(external_info[, "Frequency"]) / 2
  #   } else {
  #     min_value = 1
  #   }
  # 
  #   diplo_freqs = vector(mode = "numeric", length = len_DTs)
  #   names(diplo_freqs) = DTs
  #   
  #   for(i in 1:len_DTs){
  #     cname = unlist(strsplit(DTs[i], "\\+"))
  #     
  #     if(cname[1] %in% external_info[, "Haplotypes"]){
  #       value1 = external_info[which(external_info[, "Haplotypes"] == cname[1]), 2]
  #     } else {
  #       value1 = min_value
  #     }
  #     
  #     if(cname[2] %in% external_info[, "Haplotypes"]){
  #       value2 = external_info[which(external_info[, "Haplotypes"] == cname[2]), 2]
  #     } else {
  #       value2 = min_value
  #     }
  #     
  #     diplo_freqs[i] = value1 * value2
  #   }
  #   
  #   matval = t(apply(matval, 1, function(x){x * diplo_freqs}))
  # }

  mat = matval / rowSums(matval)
  
  return(mat)
}
