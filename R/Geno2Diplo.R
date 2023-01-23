#' @title Convert a genotype list into a diplotype list
#' 
#' @description To convert the genotype information into diplotype information. 
#'
#' @param genotypes A list. Needs to be in the format: GeneA1+GeneA2&GeneB1-GeneB2.
#' 
#' @return A list with all diplotypes that are compatible with the individual's genotype. The output of a non-ambiguous case is as: GeneA1-GeneB1+GeneA2-GeneB2 and GeneA1-GeneB2+GeneA2-GeneB1.
#' 
# #' @examples 
# #' x <- list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y <- list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' \dontrun{
# #' Geno2Haplo(Haplo2Geno(combining_genes(x, y)))
# #' }
#' 
Geno2Haplo <- function(genotypes){
  len_geno <- length(genotypes)
  
  haplos <- list(NULL)
  for(i in 1:len_geno){
    geno_donor <- genotypes[[i]]
    len_geno_donor <- length(geno_donor)
    
    haplo_all <- NULL
    for(j in 1:len_geno_donor){
      xx <- geno_donor[j]
      
      splt <- strsplit(unlist(strsplit(xx, "&")), "\\+")
      
      grid.splt <- as.matrix(expand.grid(splt))
      len_grid <- nrow(grid.splt)
      len_grid_half <- len_grid / 2
      
      for(k in 1:len_grid_half){
        haplo_all <- c(haplo_all, paste(sort(c(paste(grid.splt[k, ], collapse = "-"), 
                                               paste(grid.splt[(len_grid + 1 - k), ], collapse = "-"))), collapse = "+"))
      }
    }
    
    haplos[[i]] <- sort(unique(haplo_all))
  }
  
  return(haplos)
}
