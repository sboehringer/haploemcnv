#' @title Convert a diplotype list into a genotype list
#'
#' @description To convert the diplotype information into genotype information. 
#'
#' @param diplotypes A list. Needs to be in the format: GeneA1-GeneB1+GeneA2-GeneB2, as obtained from \code{\link{combining_genes}}
#' 
#' @return A list with all genotypes that are compatible with the individual's diplotype. The output of a non-ambiguous case is as: GeneA1+GeneA2&GeneB1+GeneB2
#' 
# #' @examples
# #' x <- list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y <- list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' \dontrun{
# #' Diplo2Geno(combining_genes(x, y))
# #' }
#' 
Diplo2Geno <- function(diplotypes){
  len_haplos <- length(diplotypes)
  nr_genes <- length(unlist(strsplit(unlist(strsplit(diplotypes[[1]][1], "\\+"))[1], "\\-")))
  
  temp_genotypes <- rec_genotypes <- as.list(rep(0, len_haplos))
  for(i in 1:len_haplos){
    xx <- diplotypes[[i]]
    len_xx <- length(xx)
    
    for(j in 1:len_xx){
      xxx <- xx[j]
      
      splt <- strsplit(unlist(strsplit(xxx, "\\+")), "\\-")
      lens_splt <- lengths(splt)
      
      if(TRUE %in% (lens_splt != nr_genes)){
        index <- which(lens_splt != nr_genes)
        len_index <- length(index)
        
        for(k in 1:len_index){
          splt[[index[k]]] <- c(rep("Comb", (nr_genes - lens_splt[index[k]])), splt[[index[k]]])
        }
      }
      
      len_splt <- length(splt[[1]])
      
      geno <- geno_temp <- NULL
      for(k in 1:len_splt){
        geno_temp <- paste(sort(unlist(lapply(splt, function(x){x[k]}))), collapse = "+")
        geno <- c(geno, geno_temp)
      }
      
      temp_genotypes[[i]][j] <- paste(geno, collapse = "&")
    }

    rec_genotypes[[i]] <- sort(unique(temp_genotypes[[i]]))  # sort(true_rec)
  }
  
  return(rec_genotypes)
}
