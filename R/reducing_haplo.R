#' @title Reducing haplotype frequencies into the K-1 genes and gene K
#' 
#' @description The haplotype frequencies will be split in two groups, the K-1 genes and gene K. Supporting function for \code{\link{LD}}
#' 
#' @param haplo_freqs A vector. The frequencies of the K-gene haplotypes. 
#' 
#' @return The haplotype frequencies of the K-1 gene combinations and the frequencies of gene K, also the total number of genes (K) and the genes in the combination (K-1) are denoted.
#' 
# #' @examples 
# #' x <- list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y <- list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' EM_out <- EM_algorithm(combining_genes(x, y))
# #' haplo_freqs <- EM_out$Frequencies[nrow(EM_out$Frequencies), ]
# #' \dontrun{
# #' reducing_haplo(haplo_freqs)
# #' }
#' 
reducing_haplo <- function(haplo_freqs){
  names_haplo <- names(haplo_freqs)
  
  splt <- strsplit(names_haplo, "\\-")
  nr_genes <- length(splt[[1]])
  nr_comb <- nr_genes - 1
  
  ## K-1 genes
  comb_haplos <- unlist(lapply(splt, function(x){paste(x[1:nr_comb], collapse = "-")}))
  un_comb_haplos <- unique(comb_haplos)
  len_un_comb_haplos <- length(un_comb_haplos)
  
  temp_haplos <- vector(mode = "numeric", len_un_comb_haplos)
  names(temp_haplos) <- un_comb_haplos
  
  for(i in 1:len_un_comb_haplos){
    index <- which(un_comb_haplos[i] == comb_haplos)
    temp_haplos[i] <- sum(haplo_freqs[index])
  }
  
  
  ## K gene
  comb_haplos_k <- unlist(lapply(splt, function(x){paste(x[nr_genes], collapse = "-")}))
  un_comb_haplos_k <- unique(comb_haplos_k)
  len_un_comb_haplos_k <- length(un_comb_haplos_k)
  
  temp_haplos_k <- vector(mode = "numeric", len_un_comb_haplos_k)
  names(temp_haplos_k) <- un_comb_haplos_k
  
  for(i in 1:len_un_comb_haplos_k){
    index <- which(un_comb_haplos_k[i] == comb_haplos_k)
    temp_haplos_k[i] <- sum(haplo_freqs[index])
  }
  
  return(list("freqs_kmin1" = temp_haplos, "freqs_k" = temp_haplos_k, "nr_genes" = nr_genes, "nr_comb" = nr_comb))
}
