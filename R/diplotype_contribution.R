#' @title Gives all diplotype configurations for a number of haplotypes 
#' 
#' @description Give all theoretically possible diplotype configurations for a specified number of haplotypes. Support function for \code{\link{EM_algorithm}} and \code{\link{information_loss}}.
#' 
#' @param nr_haplo An integer. Number of unique haplotypes for the matrix.
#'
#' @return A matrix with all diplotypes as rows and all haplotypes as columns. The values indicate how often the haplotype contributes to the diplotype.
#' 
# #' @examples 
# #' \dontrun{ 
# #' diplotype_contribution(nr_haplo = 3)
# #' }
diplotype_contribution = function(nr_haplo){
  
  nr_diplo = (nr_haplo^2 + nr_haplo) / 2
  temp = rep(0, nr_haplo)
  
  diplo_mat = matrix(0, nrow = 2, ncol = nr_diplo)
  
  for(i in 0:(nr_diplo - 1)){
    mx = floor(-0.5 + sqrt(0.25 + 2 * i))
    diplo_mat[, (i + 1)] = c(i - (mx^2 + mx) / 2, mx)
  }
  
  dt2ht = matrix(0, nrow = nr_haplo, ncol = nr_diplo)
  for(i in 1:nr_diplo){
    v = temp
    idcs = diplo_mat[, i] + 1
    e = 1; na.rm = 0
    
    v[idcs] = e
    
    dt2ht[, i] = v
  }
  
  dt2ht2 = t(dt2ht) * 2 / colSums(dt2ht)
  return(dt2ht2)
}
