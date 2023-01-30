#' @title The first derivative of the likelihood for the information Loss
#' 
#' @description Calculating the first derivative of the likelihood for the information loss. Support function of \code{\link{information_loss}}.
#' 
#' @param HHOi A vector. The values of an individual indicating how often the allele of interest occurs Homozygous-Heterozygous-Not.
#' @param allele_freq A numceric value. Allele frequencies of the allele of interest
#'
#' @return The product of the outcome of the first derivative for the individual of the allele of interest.
#' 
# #' @examples  
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' EM_out = EM_algorithm(combining_genes(x, y))
# #' 
# #' temp_mat = Reducing_mat(EM_out$Matrix, nr_gene = 1)
# #' mat_long = Geno2EM_SB(temp_mat, Init = FALSE)
# #' 
# #' last_freq = Haplo2AF(EM_out)[[nrgene]]
# #' alleles_names = names(last_freq)
# #' 
# #' ordering = dtContribMat(length(alleles_names))
# #' HomHetOth = matrix(0, nrow = len_mat, ncol = 3, dimnames = list(NULL, c("Homozygous", "Heterozygous", "Others")))
# #' 
# #' ord = ordering[, i]
# #' 
# #' HomHetOth[, "Homozygous"] = apply(mat_long, 1, function(x){sum(x[which(ord == 2)])})
# #' HomHetOth[, "Heterozygous"] = apply(mat_long, 1, function(x){sum(x[which(ord == 1)])})
# #' HomHetOth[, "Others"] = apply(mat_long, 1, function(x){sum(x[which(ord == 0)])})
# #' 
# #' allele = alleles_names[i]
# #' allele_freq = last_freq[allele]
# #' 
# #' sum(apply(HomHetOth, 1, function(x){SS(x, allele_freq)}))
#' 
SS = function(HHOi, allele_freq){
  S = as.numeric(((2 * HHOi[1] + HHOi[2]) / allele_freq) - ((2 * HHOi[3] + HHOi[2]) / (1 - allele_freq)))
  return(S * S)
}
