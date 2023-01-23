#' @title Reordering gene names of reconstruction
#' 
#' @description In order to compare the estimates from different reconstructions and/or starting genes, the haplotypes must be in the same order. Vector names will be put in the same common gene order.
#' 
#' @param names_lst A list. Haplotype notation (\emph{i.e., gene1-gene2-gene3}) of the order in which the genes are added into the reconstruction.
#' @param freqs_lst An optional list. Each list element is a separate vector of estimates. If \code{NULL}, than \code{freqs_mat} must be specified. 
#' @param freqs_mat An optional matrix. Each column is a separate reconstruction, of which the column names need to be put in the common gene order. If \code{NULL}, than \code{freqs_lst} must be specified.
#' @param end_order A character. Haplotype notation of the order in which all reconstructions must be denoted.
#'
#' @return 
#' 
# #' @examples 
# #' gene1 <- list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' gene2 <- list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' gene3 <- list(c("007+008", "007+009"), "NEG+NEG", c("008+008", "008+NEG")) 
# #'
# #' GO1 = combining_genes(combining_genes(gene1, gene2), gene3)
# #' GO1_out = EM_algorithm(GO1)
# #' 
# #' GO2 = combining_genes(combining_genes(gene2, gene3), gene1)
# #' GO2_out = EM_algorithm(GO2)
# #'
# #' names_lst = list("Gene1-Gene2-Gene3", "Gene2-Gene3-Gene1")
# #' freqs_lst = list(GO1_out$Freq[nrow(GO1_out$Freq), ], GO2_out$Freq[nrow(GO2_out$Freq), ])
# #' end_order = "Gene3-Gene2-Gene1"
# #'
# #' \dontrun{
# #' reordering_grouping_results(names_lst, freqs_lst, end_order = end_order)
# #' }
reordering_grouping_results <- function(names_lst, freqs_lst = NULL, freqs_mat = NULL, end_order){
  
  N <- length(names_lst)
  
  splt_end_order <- unlist(strsplit(end_order, "\\-"))
  n1 <- length(splt_end_order)
  
  ordering <- as.list(rep(0, N))
  for(i in 1:N){
    splt_name <- unlist(strsplit(names_lst[[i]], "\\-"))
    
    for(j in 1:n1){
      ordering[[i]][j] <- which(splt_end_order[j] == splt_name)
    }
  }
  
  if(!is.null(freqs_lst)){
    for(i in 1:N){
      freqs <- freqs_lst[[i]]
      fnames <- names(freqs)
      n2 <- length(fnames)
      
      for(j in 1:n2){
        fname <- unlist(strsplit(fnames[j], "\\-"))
        
        fnames[j] <- paste(fname[ordering[[i]]], collapse = "-")
      }
      
      names(freqs) <- fnames
      freqs_lst[[i]] <- freqs
    }
    
    out = freqs_lst
  } else if(!is.null(freqs_mat)){
    for(i in 1:N){
      mat = freqs_mat[[i]]
      mnames = colnames(mat); n2 = length(mnames)
      
      for(j in 1:n2){
        mname = unlist(strsplit(mnames[j], "\\+"))
        
        for(k in 1:2){
          mname[k] = paste(unlist(strsplit(mname[k], "\\-"))[ordering[[i]]], collapse = "-")
        }
        
        mnames[j] = paste(mname, collapse = "+")
      }
      
      colnames(mat) = mnames
      freqs_mat[[i]] = mat
    }
    
    out = freqs_mat
  } else {
    stop("Supply either a list or a matrix.")
  }

  
  return(out)
}
