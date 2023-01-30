#' @title Calculate the Information loss for a gene
#' 
#' @description Calculates the amount of information loss within a reconstruction, as discussed in \cite{}.
#' 
#' @param EM_out An optional EM-algorithm analysis of \code{lst}, as provided by \code{\link{EM_algorithm}}.
#' @param HTFs An optional vector. Haplotype frequencies for all haplotypes, as are estimated in \code{\link{EM_algorithm}}. This was haplo_freqs
#' @param EM_mat An optional matrix. Probabilities for each compatible diplotype in the dataset, as can be provided by \code{\link{EM_algorithm}}.
#' @param haplotype_level A logical scalar. Whether or not the supplied EM_matrix is on haplotype level or on diplotype level (\code{FALSE} is default).
#' @param nr_gene An integer. The number of genes that are in the reconstruction.
#'
#' @return An matrix with the calculated I_{Y}, I_{X}, ratio I_{Y} / I_{X} and the effective sample size (I_{Y} $theta (1 - theta)$) for the specified gene
#' 
# #' @examples 
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' EM_out = EM_algorithm(combining_genes(x, y))
# #' information_loss(EM_out)
#' 
information_loss = function(EM_out = NULL, HTFs = NULL, EM_mat = NULL, haplotype_level = FALSE, nr_gene = 1){
  
  if(!is.null(EM_out)){
    EM_mat = EM_out$Matrix
    HTFs = EM_out$Frequencies[nrow(EM_out$Frequencies), ]
  } 
  
  len_mat = nrow(EM_mat)
  
  if(haplotype_level == TRUE){
    mat_long = Geno2EM_SB(dat = EM_mat, Init = FALSE)
    last_freq = HTFs

  } else {
    temp_mat = reducing_mat(EM_mat = EM_mat, nr_gene = nr_gene)
    mat_long = Geno2EM_SB(dat = temp_mat, Init = FALSE)
    last_freq = Haplo2AF(haplotypes = HTFs)[[nr_gene]]
    
  }
  

  haplo_names = names(last_freq)
  len_haplo_names = length(haplo_names)
  
  ordering = diplotype_contribution(len_haplo_names)  # dtContribMat(len_haplo_names)
  
  HomHetOth = matrix(0, nrow = len_mat, ncol = 3, dimnames = list(NULL, c("Homozygous", "Heterozygous", "Others")))
  Info_Loss_mat = matrix(0, nrow = len_haplo_names, ncol = 5, dimnames = list(haplo_names, c("IY", "IX", "Ratio", "Haplo_freq", "Eff_n")))
  for(i in 1:len_haplo_names){
    ord = ordering[, i]
    
    HomHetOth[, "Homozygous"] = apply(mat_long, 1, function(x){sum(x[which(ord == 2)])})
    HomHetOth[, "Heterozygous"] = apply(mat_long, 1, function(x){sum(x[which(ord == 1)])})
    HomHetOth[, "Others"] = apply(mat_long, 1, function(x){sum(x[which(ord == 0)])})
  
    haplo = haplo_names[i]
    haplo_freq = last_freq[haplo]
    
    IY = sum(apply(HomHetOth, 1, function(x){SS(x, haplo_freq)}))
    IX = sum(apply(HomHetOth, 1, function(x){B(x, haplo_freq)}))
    
    Info_Loss_mat[i, 1] = IY
    Info_Loss_mat[i, 2] = IX
    Info_Loss_mat[i, 3] = IY / IX
    
    Info_Loss_mat[i, 4] = haplo_freq
    
    Info_Loss_mat[i, 5] = IY * haplo_freq * (1 - haplo_freq)
  }
  
  return(Info_Loss_mat)
}
