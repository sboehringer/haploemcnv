#' @title Calculate allele frequencies from haplotype frequencies
#' 
#' @description To calculate allele frequencies for all genes based on the haplotype frequencies. Either an EM-algorithm analysis or the haplotype frequencies need to be supplied.
#' 
#' @param EM_out An optional EM-algorithm analysis of \code{lst}, as provided by \code{\link{EM_algorithm}}.
#' @param haplotypes An optional vector. Containing the haplotype frequencies. 
#' 
#' @return An list with in each list element the allele frequencies of the different genes. The first element contains the AFs of the first gene, the second element of the second gene etc.
#'
# #' @examples
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' EM_out = EM_algorithm(combining_genes(x, y))
# #' \dontrun{
# #' Haplo2AF(EM_out)
# #' }
#' 
Haplo2AF = function(EM_out = NULL, haplotypes = NULL){
  
  if(!is.null(haplotypes)){
    haplo_freqs = haplotypes
  } 
  
  if(!is.null(EM_out)){
    haplo_freqs = EM_out$Freq[nrow(EM_out$Freq), ]
  }
  
  len_haplos = length(haplo_freqs)
  
  allele_names = strsplit(names(haplo_freqs), "\\-")
  nrGenes = max(lengths(allele_names))
  
  AFs = as.list(rep(NA, nrGenes))
  for(i in 1:nrGenes){
    alleles = unlist(lapply(allele_names, function(x){x[i]}))
    
    un_alleles = sort(unique(alleles))
    len_un_alleles = length(un_alleles)
    
    AFs[[i]] = rep(0, len_un_alleles)
    names(AFs[[i]]) = un_alleles
    
    for(j in 1:len_haplos){
      name = allele_names[[j]][i]
      value = as.numeric(haplo_freqs[j])
      
      index = which(names(AFs[[i]]) == name)
      
      AFs[[i]][index] = AFs[[i]][index] + value
    }
  }
  
  return(AFs)
}
