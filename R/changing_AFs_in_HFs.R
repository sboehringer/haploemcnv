#' @title Alter compound allele frequencies in haplotypes
#' 
#' @description For the profile EM-algorithm the compound allele frequency needs to be kept at bay. Support function for \code{\link{EM-algorithm}}.
#' 
#' @param haplotypes A vector. Haplotypes and their probabilities.
#' @param comb_vals A vector. Compound haplotypes and their \emph{true} frequency
#' @param nonneg a logical scalar. Whether or not to make sure that there are no negative probabilities due to the resetting the compound haplotype frequency.
#' 
#' @return A vector of the haplotypes and some altered probabilities.
#' 
changing_AFs_in_HFs <- function(haplotypes, comb_vals, nonneg = TRUE){
  
  ##
  names_comb_vals <- names(comb_vals)
  len_comb_vals <- length(comb_vals)
  
  for(i in 1:len_comb_vals){
    comb_vals[i] <- comb_vals[i]^(1 / length(unlist(strsplit(str_sub(names_comb_vals[i], start = 6), ""))))
  }
  ##
  
  hnames <- names(haplotypes)
  len_hnames <- length(hnames)

  hnames_splt <- strsplit(hnames, "\\-")
  nr_genes <- length(unlist(strsplit(hnames[1], "\\-")))
  
  AFs <- Haplo2AF(haplotypes = haplotypes)
  len_AFs <- length(AFs)
  
  for(i in 1:len_AFs){
    AF_gene <- AFs[[i]]
    len_AF_gene <- length(AF_gene)
    
    for(j in 1:len_AF_gene){
      
      if(names(AF_gene[j]) %in% names_comb_vals){
        index <- which(names(AF_gene[j]) == names_comb_vals)
        AFs[[i]][j] <- comb_vals[index]
        
      }
    }
  } 

  if(len_AFs == 1){
    out <- unlist(AFs)
    
  } else if(len_AFs == 2){
    new_hname_vals <- as.vector(matrix(0, nrow = 1, ncol = len_hnames))
    names(new_hname_vals) <- hnames
    
    for(i in 1:len_hnames){
      hname <- hnames[i]
      hname_splt <- unlist(strsplit(hname, "\\-"))
      
      new_hname_vals[i] <- AFs[[1]][hname_splt[1]] * AFs[[2]][hname_splt[2]]
    }
    
    LD_values <- LD(haplotypes = haplotypes)[, 2]
    
#    index <- unlist(lapply(strsplit(names(LD_values), "\\-"), function(x){TRUE %in% (names(comb_vals[1]) %in% x)}))
    index <- unlist(lapply(strsplit(names(LD_values), "\\-"), function(x){TRUE %in% (names_comb_vals %in% x)}))
    LD_values[index] <- 0
    
    new_hnames <- LD_values + new_hname_vals

    out <- new_hnames[hnames]
    
  } else {
    haplo_subset <- list(NULL)
    haplo_subset_changed <- list(NULL)
    
    i = nr_genes

    haplo_subset[[i]] <- list("freqs_kmin1" = haplotypes)
    while(i > 1){
      haplo_subset[[i - 1]] <- reducing_haplo(haplo_freqs = haplo_subset[[i]]$freqs_kmin1)

      i <- i - 1
    }
    haplo_subset_changed[[i]] <- list("freqs_kmin1" = AFs[[i]], "freqs_k" = NULL)
    
    
    for(i in 1:(len_AFs - 1)){
      hnames_kmin1 <- names(haplo_subset[[i + 1]]$freqs_kmin1)
      len_hnames_kmin1 <- length(hnames_kmin1)
      
      haplo_subset_changed[[i + 1]] <- list("freqs_kmin1" = NULL)
      for(j in 1:len_hnames_kmin1){
        hname_kmin1 <- unlist(strsplit(hnames_kmin1[j], "\\-"))
        hname_kmin1 <- c(paste(hname_kmin1[-length(hname_kmin1)], collapse = "-"), hname_kmin1[length(hname_kmin1)])
        
        haplo_subset_changed[[i + 1]]$freqs_kmin1[j] <- haplo_subset_changed[[i]]$freqs_kmin1[hname_kmin1[1]] * AFs[[i + 1]][hname_kmin1[2]]
      }
      names(haplo_subset_changed[[i + 1]]$freqs_kmin1) <- hnames_kmin1
      
      if(i == 1){
        LD_values <- LD(haplotypes = haplo_subset[[i + 1]]$freqs_kmin1)[, 2]
      } else {
        LD_values <- LD(haplotypes = haplo_subset[[i + 1]]$freqs_kmin1, three_plus = TRUE)[, 2]
      }
      
      index <- unlist(lapply(strsplit(names(LD_values), "\\-"), function(x){TRUE %in% (names_comb_vals %in% x)}))
      LD_values[index] <- 0
      
      haplo_subset_changed[[i + 1]]$freqs_kmin1 <- LD_values + haplo_subset_changed[[i + 1]]$freqs_kmin1
    }
    
#    all.equal(haplo_subset_changed[[i + 1]]$freqs_kmin1, new_haplotypes)  # Of course only the same if there are no changes due to comb_vals
    out <- haplo_subset_changed[[i + 1]]$freqs_kmin1[hnames]
  }

  if(nonneg){
    out = pmax(out, 0)
  }
  
  return(out)
}
