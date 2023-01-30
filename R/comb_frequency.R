#' @title Calculating the \emph{true} compound haplotype frequency
#' 
#' @description Calculating the \emph{true} compound haplotype frequency based on a previous selection of which genes are grouped and which are retained. 
#' 
#' @details Because of coding each \emph{x}-locus haplotype is required to have \emph{x} separate gene copies. Thus when a \emph{x}-locus hapltoype is combined into the compound haplotype, each gene copy is replaced by \emph{Comb1:x}. Since a haplotype frequency equals the product of its individual alleles. The compound haplotype frequency of a multi-locus haplotype needs to anticipate on that, requiring adjusting via \code{comb_adjust}. 
#' 
#' @param EM_out An optional EM-algorithm analysis of \code{lst}, as provided by \code{\link{EM_algorithm}}. If \code{NULL}, than \code{haplotypes} must be specified.
#' @param haplotypes An optional vector. The haplotypes with their frequency, can be obtained by \code{\link{EM_algorithm}}. If \code{NULL}, than \code{EM_out} must be specified.
#' @param remain An optional vector. Haplotypes that are retained.
#' @param combine An optional vector. Haplotypes that are grouped.
#' @param comb_adjust A logical scalar. Whether or not the estimated \emph{true} compound haplotype frequency must be adjusted for the fact that the haplotype consists of multiple loci.
#' 
#' @return A single value, representing the \emph{true} compound haplotype frequency.
#' 
comb_frequency = function(EM_out = NULL, haplotypes = NULL, remain = NULL, combine = NULL, comb_adjust = FALSE){
  
  if(!is.null(EM_out)){
    haplotypes = EM_out$Frequencies[nrow(EM_out$Frequencies), ]
  }
  
  if(!is.null(remain)){
    comb_val = sum(haplotypes[!(names(haplotypes) %in% remain)])
  } else if(!is.null(combine)){
    comb_val = sum(haplotypes[names(haplotypes) %in% combine])
  }

  nr_genes = length(unlist(strsplit(names(haplotypes[1]), "\\-")))
  if(comb_adjust && nr_genes > 1){
    comb_val = comb_val^(1 / nr_genes)
    # comb_val = comb_val^(1 / (nr_genes * (nr_genes - 1)))
  }
  
  return(comb_val)
}
