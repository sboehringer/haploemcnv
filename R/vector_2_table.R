#' @title Making a table from a vector
#' 
#' @description Making from a two-gene haplotype vector a table with the genes separately as rows and columns. 
#' 
#' @param vect A vector. Two-gene haplotypes and their frequencies.
#' 
#' @return A table with the alleles of the first gene as rows and the alleles of the second gene as columns.
#' 
# #' @examples 
# #' vect <- c("001-002" = 0.2, "001-003" = 0.1, "001-NEG" = 0.2, "NEG-002" = 0.08, "NEG-003" = 0.02, "NEG-NEG" = 0.4)
# #' 
# #' \dontrun{
# #' vector_2_table(vect)
# #'}
#'
vector_2_table <- function(vect){
  
  len_vect <- length(vect)

  cnames <- names(vect)
  cnames_splt <- strsplit(cnames, "\\-")
  
  nr_genes <- lengths(cnames_splt)[1]
  if(nr_genes > 2){
    cnames_splt <- lapply(cnames_splt, function(x){c(paste(x[1:(nr_genes - 1)], collapse = "-"), x[3])})
  }
  
  gene1 <- unique(unlist(lapply(cnames_splt, function(x){x[1]})))
  len_gene1 <- length(gene1)
  
  gene2 <- unique(unlist(lapply(cnames_splt, function(x){x[2]})))
  len_gene2 <- length(gene2)
  
  
  tab <- matrix(0, nrow = len_gene1, ncol = len_gene2, byrow = TRUE, dimnames = list(gene1, gene2))
  for(i in 1:len_vect){
    cname <- cnames_splt[[i]]
    
    index1 <- which(gene1 == cname[1])
    index2 <- which(gene2 == cname[2])
    
    tab[index1, index2] <- vect[i]
  }
  
  return(tab)
}
