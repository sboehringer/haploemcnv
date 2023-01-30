#' Code for converting my dataset to Stefan's functions, also support function for \code{\link{information_loss}} and \code{\link{likelihood_haplotypes}}.
#'
#' @param dat The dataset to convert.
#' @param Init Logical. Whether or not make an 0/1-initial value matrix. TRUE is default.
#' 
#' @return The output is a matrix with all possible genotype combinations (not only those that are observed) in a different columns
#' 
# #' @examples 
# #' dat <- list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' 
# #' \dontrun{
# #' Geno2EM_SB(dat)
# #' }
#' 
Geno2EM_SB <- function(dat, Init = TRUE){
  
  if(Init == TRUE){
    len_dat <- length(dat)
    mat <- EM_init_mat(dat)
  } else {
    len_dat <- nrow(dat)
    mat <- dat
  }
  
  cnames <- colnames(mat)
  len_cnames <- length(cnames)
  
  hts <- sort(unique(unlist(strsplit(cnames, "\\+"))))
  Nhts <- length(hts)
  
  Mdt2ht <- diplotype_contribution(nr_haplo = Nhts)  # dtContribMat(Nhts)
  colnames(Mdt2ht) = hts
  
  nrow_Mdt2ht <- nrow(Mdt2ht)
  ncol_Mdt2ht <- ncol(Mdt2ht)
  
  all_colnames <- NULL
  for(i in 1:nrow_Mdt2ht){
    comb_colname <- NULL
    
    for(j in 1:ncol_Mdt2ht){
      
      if(Mdt2ht[i, j] == 1){
        comb_colname <- c(comb_colname, names(Mdt2ht[i, j]))
      }
      
      if(Mdt2ht[i, j] == 2){
        comb_colname <- c(names(Mdt2ht[i, j]), names(Mdt2ht[i, j]))
      }
    }
    
    all_colnames <- c(all_colnames, paste(sort(comb_colname), collapse = "+"))
  }
  
  ord_hts <- matrix(0, nrow = len_dat, ncol = nrow_Mdt2ht, dimnames = list(NULL, all_colnames))
  for(i in 1:len_cnames){
    index = which(all_colnames == cnames[i])
    ord_hts[, index] <- mat[, i]
  }
  
  return(ord_hts)
}
