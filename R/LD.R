#' @title linkage disequilibrium calculation
#' 
#' @description Linkage disequilibrium is calculated between two genes or between a combination of genes and the latest added gene. 
#'
#' @param EM_out An EM-algorithm analysis of at least two genes, as provided by \code{\link{EM_algorithm}}.
#' @param haplotypes An optional vector. Haplotypes and their probabilities, must be haplotypes consisting of at least two loci.
#' @param three_plus A logical scalar. Whether or not the linkage disequilibrium is calculated between the K-1 genes and gene K, or between gene 1 and gene 2 (\code{FALSE} is default). 
#' 
#' @return An matrix with the D', delta_{ij}, the R2, the observed and the expected haplotype frequency as separate columns and all possible gene combinations as different rows.
#' 
# #' @examples
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' EM_out = EM_algorithm(combining_genes(x, y))
# #' \dontrun{
# #' LD(EM_out)
# #' }
#' 
LD = function(EM_out = NULL, haplotypes = NULL, three_plus = FALSE){
  
  if(!is.null(haplotypes)){
    names_haplo = sort(names(haplotypes))
    haplo_freqs = haplotypes[names_haplo]
    
    allele_freqs = Haplo2AF(haplotypes = haplo_freqs)
  }
  
  if(!is.null(EM_out)){
    haplo_freqs = EM_out$Frequencies[nrow(EM_out$Frequencies), ]
    names_haplo = names(haplo_freqs)
    
    allele_freqs = Haplo2AF(EM_out = EM_out)
  }

  if(three_plus == FALSE){
    
    gene1_freqs = allele_freqs[[1]]
    gene2_freqs = allele_freqs[[2]]
  } else {
    
    red_haplos = reducing_haplo(haplo_freqs)
    
    temp_haplos = red_haplos$freqs_kmin1
    nr_genes = red_haplos$nr_genes
    nr_comb = red_haplos$nr_comb
    
    
    gene1_freqs = temp_haplos
    gene2_freqs = allele_freqs[[nr_genes]]
  }
  
  
  len_haplo = length(haplo_freqs)
  haplo_splt = strsplit(names_haplo, "\\-")
  
  Linkage = matrix(0, nrow = len_haplo, ncol = 7, dimnames = list(names_haplo, c("D_prime", "delta_ij", "R2", "Obs", "Exp", "p_i", "q_j")))
  
  for(i in 1:len_haplo){
    rho_ij = haplo_freqs[[i]]
    
    if(three_plus == FALSE){
      
      p_i = gene1_freqs[haplo_splt[[i]][1]]
      q_j = gene2_freqs[haplo_splt[[i]][2]]
    } else {
      
      p_i = gene1_freqs[paste(haplo_splt[[i]][1:nr_comb], collapse = "-")]
      q_j = gene2_freqs[haplo_splt[[i]][nr_genes]]
    }
    
    pq_ij = p_i * q_j
    
    Linkage[i, "delta_ij"] = rho_ij - pq_ij
    
    if(Linkage[i, "delta_ij"] < 0){
      d_max = -(min(0, (1 - p_i - q_j)) - pq_ij)
      
    } else {
      d_max = min(p_i, q_j) - pq_ij
      
    }
    
    Linkage[i, "D_prime"] = Linkage[i, "delta_ij"] / d_max
    
    
    p_i1minp_i = p_i * (1 - p_i)
    q_j1minq_j = q_j * (1 - q_j)
    Linkage[i, "R2"] = (Linkage[i, "delta_ij"]**2) / (p_i1minp_i * q_j1minq_j)
    
    Linkage[i, "Obs"] = rho_ij
    Linkage[i, "Exp"] = pq_ij
    
    Linkage[i, "p_i"] = p_i
    Linkage[i, "q_j"] = q_j
  }
  
  return(Linkage)
}
