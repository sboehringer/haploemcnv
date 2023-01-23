#' @title Reducing a diplotype probability matrix to the diplotypes of a single gene 
#' 
#' @description Reduces a multi-locus updated diplotype probability matrix to a single-locus diplotype probability matrix. Support function of \code{\link{InformationLoss}}.
#' 
#' @param EM_mat A matrix. The matrix updated via \code{\link{EM_algorithm}} for a combination of genes.
#' @param nr_gene An integer or a vector. To which genes \code{EM_mat} needs to be reduced.
#' @param haplo A logical scalar. Whether or not the matrix should be reduced to a haplotype probability matrix (\code{FALSE} is default). 
#' @param noCNV A logical scalar, whether or not we want to split the CNVs into separate alleles (\code{FALSE} is default).
#' 
#' @return A EM-algorithm matrix for the assigned gene(s).
#' 
# #' @examples
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' EM_out = EM_algorithm(combining_genes(x, y))
# #' EM_mat = EM_out$Matrix
# #' 
# #' reducing_mat(EM_mat, nr_gene = 1)
#'
reducing_mat = function(EM_mat, nr_gene, haplo = FALSE, noCNV = FALSE){
 
  cnames = colnames(EM_mat); len_mat = nrow(EM_mat)
  splt = strsplit(cnames, "\\+"); len_splt = length(splt)
  
  if(identical(1:length(unlist(strsplit(splt[[1]][[1]], "\\-"))), nr_gene)){
    stop("Matrix is already in the format of nr_gene...")
  }
  
  for(i in 1:len_splt){
    splttd = strsplit(splt[[i]], "\\-")
    cnames[i] = paste(unlist(lapply(splttd, function(x){paste(x[nr_gene], collapse = "-")})), collapse = "+")
  }
  un_cnames = sort(unique(cnames)); len_un_cnames = length(un_cnames)

  
  if(TRUE){
    unordered_mat = sapply(un_cnames, function(x){rowSums(as.matrix(EM_mat[, x == cnames]))})
  } else {
    unordered_mat = mat; colnames(unordered_mat) = cnames
  }
  
  cnames_unordered = colnames(unordered_mat); cnames_ordered = unlist(lapply(strsplit(cnames_unordered, "\\+"), function(x){paste(sort(x), collapse = "+")})) 
  un_cnames_ordered = sort(unique(cnames_ordered)); len_un_cnames_ordered = length(un_cnames_ordered)
  
  Repeat_mat = matrix(0, nrow = length(cnames_unordered), ncol = len_un_cnames_ordered, dimnames = list(cnames_unordered, un_cnames_ordered))
  for(j in 1:len_un_cnames_ordered){
    Repeat_mat[, j] = un_cnames_ordered[j] == cnames_ordered
  }
  ordered_mat = unordered_mat %*% Repeat_mat
  
  
  
  if(haplo){  
    un_cnames_ordered_splt = strsplit(un_cnames_ordered, "\\+")
    alleles = sort(unique(unlist(un_cnames_ordered_splt))); len_alleles = length(alleles)
    
    Repeat_mat = matrix(0, nrow = len_un_cnames_ordered, ncol = len_alleles, dimnames = list(un_cnames_ordered, alleles))
    for(i in 1:len_alleles){
      Repeat_mat[, i] = unlist(lapply(un_cnames_ordered_splt, function(x){sum(x %in% alleles[i])}))
    }
    
    if(noCNV & length(nr_gene) == 1){
      alleles_splt = strsplit(alleles, "\\^"); noCNV_alleles = sort(unique(unlist(alleles_splt))); len_noCNV_alleles = length(noCNV_alleles)
      
      Repeat_mat2 = matrix(0, nrow = len_alleles, ncol = len_noCNV_alleles, dimnames = list(alleles, noCNV_alleles))
      for(i in 1:len_noCNV_alleles){
        Repeat_mat2[, i] = unlist(lapply(alleles_splt, function(x){sum(x %in% noCNV_alleles[i])}))
      }
      
      Repeat_mat = Repeat_mat %*% Repeat_mat2
    } else{
      if(noCNV){
        stop("Can perform the noCNV code only when there is a single locus\n")
      }
    }
    
    ordered_mat = ordered_mat %*% Repeat_mat
  }
  
  
  return(ordered_mat)
}
