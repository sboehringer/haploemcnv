#' @title Selecting the relevant haplotype frequencies from an external dataset
#'
#' @description To calculate the haplotype frequencies for only a subset of genes. By discaring irrelevant genes, dubble haplotypes will occur and need to aggregate over them. 
#'
#' @param ext_dat A matrix. Allelic content of each gene in a separate column and the frequencies in a separate column.
#' @param genes A vector. Genes for which the haplotype frequency needs to be determined.
#' 
#' @return A matrix with two columns, the first column denotes the unique haplotypes and the second column the corresponding frequencies.
#' 
# #' @examples 
# #' ext_dat <- matrix(c("001", "NEG", "NEG", 0.2, "001", "NEG", "003", 0.35, "NEG", "002", "003", 0.15, "001", "002", "003", 0.3), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("Gene1", "Gene2", "Gene3", "Frequency")))
# #' genes <- c("Gene1", "Gene2")
# #' 
# #' \dontrun{
# #' external_information_haplotypes(ext_dat, genes)
# #' }
#' 
external_information_haplotypes <- function(ext_dat, genes){
  
  len_genes <- length(genes)
  
  ordering <- c(genes, "Frequency")
  subset <- ext_dat[, colnames(ext_dat) %in% ordering][, ordering]

  haplotypes <- apply(subset, 1, function(x){paste(x[1:len_genes], collapse = "-")})
  haplotypes = unlist(lapply(strsplit(haplotypes, "\\+"), function(x){paste(x, collapse = "^")}))

  dat <- data.frame(cbind("Haplotypes" = haplotypes, "Frequency" = subset[, "Frequency"]))
  dat$Frequency <- as.numeric(dat$Frequency)
  
  red_dat <- aggregate(Frequency ~ Haplotypes, data = dat, sum)
  
  red_dat = red_dat[order(red_dat$Frequency, decreasing = TRUE), ]
  
  return(red_dat)
}
