#' @title Kullback-Leibler divergence measure
#' 
#' @description Calculation of the Kullback-Leibler divergence measure between two vectors. Two vectors \code{P} and \code{Q} will be standardized, so not necessarily need to sum to 1. Actual KLD measure is executed by \code{\link[philentropy]{KL}}.
#' 
#' @param P A vector. Values to compare with \code{Q}.
#' @param Q A vector. Correct values, need to consists of the same value names as \code{P}. 
#' @param comb_add A logical scalar. Whether or not the cumulative sum of the vectors \code{P} and \code{Q} need to be composed of a compound haplotype, or that both vectors need to be normalized.
#' @param base An integer. Base value for the logaritmic transformation in \code{\link[philentropy]{KL}}.
#'
#' @return The kullback-leibler divergence measure between \code{P} and \code{Q}.
#' 
KL_divergence = function(P, Q, comb_add = FALSE, base = 2){
  
  if(base == 0){
    unit = "log"
  } else if(base == 2){
    unit = "log2"
  } else if(base == 10){
    unit = "log10"
  } else {
    stop("Choose for base 0, 2 or 10")
  }
  
  cnames = names(P)
    
  splt = strsplit(cnames, "\\-")
  lens_splt = lengths(splt)
    
  gene_sizes = unique(lens_splt)
  nr_genes = length(gene_sizes)


  newP = newQ = NULL
  if(comb_add){
    for(i in gene_sizes){
      tempP = P[lens_splt == i]
      
      comb_valP = (1 - sum(tempP))
      names(comb_valP) = paste(rep(paste("Comb", paste(1:i, collapse = ""), sep = ""), i), collapse = "-")
      
      newP = c(newP, tempP, comb_valP)
      
      
      tempQ = Q[lens_splt == i]
      
      comb_valQ = (1 - sum(tempQ))
      names(comb_valQ) = paste(rep(paste("Comb", paste(1:i, collapse = ""), sep = ""), i), collapse = "-")
      
      newQ = c(newQ, tempQ, comb_valQ)
    }
    
  } else {
    for(i in gene_sizes){
      tempP = P[lens_splt == i]
      tempP = tempP / sum(tempP) / nr_genes
      
      newP = c(newP, tempP)
      
      
      tempQ = Q[lens_splt == i]
      tempQ = tempQ / sum(tempQ) / nr_genes
      
      newQ = c(newQ, tempQ)
    }
  }
  
  

  new_cnames = names(newP)
  
  new_splt = strsplit(new_cnames, "\\-")
  lens_new_splt = lengths(new_splt)
    
  KL_measures = NULL
  index = 1
  for(i in gene_sizes){
    
    subP = newP[lens_new_splt == i]
    subQ = newQ[lens_new_splt == i]
      
#    KL_measures[index] = KL_measure(P = subP, Q = subQ, base = base)
    KL_measures[index] = KL(rbind(P = subP, Q = subQ), unit = unit)
    
    index = index + 1
  }   
  
  return(KL_measures)
}
