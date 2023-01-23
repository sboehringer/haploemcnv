#' @title Calculation of Delta_{ij}
#' 
#' @description Support function for \code{\link{KIR_simulation_preparation_real}}.
#' 
#' @param LD_out The linkage disequilibrium matrix for the full K genes, as provided by \code{\link{LD}}. 
#' @param p_i An optional vector. Haplotype frequencies for the K-1 genes.
#' @param q_j An optional vector. Haplotype frequencies for gene K.
#' 
#' @return Vector with Delta_{ij} values.
#' 
Dprime_to_Delta = function(LD_out, p_i = NULL, q_j = NULL){
  
  D_prime = LD_out[, "D_prime"]
  rnames = rownames(LD_out)
  
  nr_genes = lengths(strsplit(rnames, "\\-"))[1]
  nr_genes_m1 = nr_genes - 1
  
  if(is.null(p_i)){
    p_i = LD_out[, "p_i"]
    names(p_i) = unlist(lapply(strsplit(names(p_i), "\\-"), function(x){paste(x[1:nr_genes_m1], collapse = "-")}))
    p_i = p_i[!duplicated(names(p_i))]
  }
  if(is.null(q_j)){
    q_j = LD_out[, "q_j"]
    names(q_j) = unlist(lapply(strsplit(names(q_j), "\\-"), function(x){paste(x[nr_genes], collapse = "-")}))
    q_j = q_j[!duplicated(names(q_j))]
  }
  
  
  nr_haplos = nrow(LD_out)
  Deltaij = rep(0, nr_haplos); names(Deltaij) = rnames
  for(i in 1:nr_haplos){
    Dprime = D_prime[i]
    Dprime_name = names(Dprime)
    
    Dprime_name_unlist = unlist(strsplit(Dprime_name, "\\-"))
    Dprime_name_splt = c(paste(Dprime_name_unlist[1:nr_genes_m1], collapse = "-"), Dprime_name_unlist[nr_genes])
    
    pi = p_i[Dprime_name_splt[1]]
    qj = q_j[Dprime_name_splt[2]]
    
    pqij = as.numeric(pi * qj)
    
    Dmax = 0
    if(Dprime < 0){
      Dmax = -(min(0, (1 - pi - qj)) - pqij)
    } else {
      Dmax = min(pi, qj) - pqij
    }
    
    Deltaij[i] = Dprime * Dmax
  }
  

  return(Deltaij)
}
