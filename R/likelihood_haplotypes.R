#' @title Calculates the log-likelihood for the EM-algorithm
#' 
#' @description Calculates the log-likelihood based on the data analyzed with \code{\link{EM_algorithm}}.
#'
#' @param mat A matrix. The probability vectors for each individual, can be obtained via \code{\link{EM_algorithm}}.
#' @param haplos_i A vector. The haplotype frequencies for the haplotypes in \code{mat}, can be obtained via \code{\link{EM_algorithm}}.
#' @param indicator An optinal matrix. How much each haplotype (row) contributes to a diplotype (columns). If \code{NULL}, an indicator matrix will be made.
# #' @param normalize An logical scalar. \lars{Not sure what it does}.
#' 
#' @return The (log-)likelihood of the EM-algorithm iteration and the indicator matrix made (or supplied as argument).
#' 
# #' @examples
# #' lst = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' lst_out = EM_algorithm(lst)
# #' 
# #' mat = lst_out$Matrix
# #' haplos_i = lst_out$Frequencies[nrow(lst_out$Frequencies), ]
# #' \dontrun{
# #' likelihood_haplotypes(mat, haplos_i)
# #' }
#' 
likelihood_haplotypes = function(mat, haplos_i, indicator = NULL){  # , normalize = FALSE){
  
  nr_donors = nrow(mat)
  
  haplos = names(haplos_i); nr_haplos = length(haplos)
  diplotypes = colnames(mat); nr_diplotypes = length(diplotypes); diplotypes_splt = strsplit(diplotypes, "\\+")
  
  if(is.null(indicator)){
    indicator = matrix(NA, nrow = nr_haplos, ncol = nr_diplotypes, dimnames = list(haplos, diplotypes))
    for(h in 1:nr_haplos){
      haplo = haplos[h]
      
      haplos_i[h] = ifelse(haplos_i[h] == 0, .Machine$double.xmin, haplos_i[h])
      
      for(d in 1:nr_diplotypes){
        diplo = diplotypes_splt[[d]]
        
        indicator[h, d] = sum(diplo %in% haplo)
      }
    }
  }

  
  # if(normalize){
  #   save_mat = matrix(0, nrow = nr_haplos, ncol = (nr_diplotypes + 1), dimnames = list(haplos, c(diplotypes, "SUM")))
  #   save_mat_lst = as.list(rep(0, nr_donors))
  #   
  #   for(i in 1:nr_donors){
  #     save_mat_lst[[i]] = save_mat
  #   }
  # }

  
  Likeli = 0
  for(h in 1:nr_haplos){
    h_indicator = indicator[h, ]
    
    hi = 0
    for(i in 1:nr_donors){
      diplo_part = mat[i, ] * h_indicator
      
      hi = hi + sum(diplo_part)
 
      if(normalize){
        save_mat_lst[[i]][h, ] = c(diplo_part, sum(diplo_part))
      }

    }

    Likeli = Likeli + log(haplos_i[h]) * hi
  }
  
  # if(normalize){
  #   check = all.equal(unlist(lapply(save_mat_lst, function(x){sum(x[, "SUM"])})), rep(2, nr_donors))  # , tolerance = 1e-100)
  #   
  #   if(!check){
  #     cat("something went wrong, some donors do not equal 2 but have a", check)
  #   }
  # }
  
  return(list("Likelihood" = Likeli, "Indicator" = indicator))
}
