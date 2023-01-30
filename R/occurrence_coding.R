#' @title Making an occurrence matrix of haplotype occurrences 
#' 
#' @description To make an occurrence matrix indicating how often each haplotype occurs in the individual's diplotypes. 
#' 
#' @param lst An optional list with each individual having a separate list elements with all its compatible diplotypes, as provided by \code{\link{all_options}}
#' @param EM_out An optional list of the EM-algorithm output of the \code{lst} compatible diplotype list
#' @param EM_out_mat The optional matrix of the EM-algorithm of the \code{lst} compatible diplotype list
#' @param haplos_keep An optional vector. For which haplotypes the occurrence matrices need to be made. if \code{NULL} the occurrence matrix will be based on all supplied haplotypes
#' @param gene_indication A logical scalar. Whether or not it is indicated from which gene combination the allele or haplotype is (default is \emph{TRUE}).
#' @param noCNV a logical scalar. Whether or not the copy number variation haplotypes should be accommodated in the allelic predictors or should be kept as separate predictors. Only functional when \code{regression} is not \code{NULL}.
#' @param CNV_option A character. Indicating what to do with the copy number variation alleles with the inclusion of the outcome (either "first", "second", or "third"). 
#' 
#' @return A matrix with each individual having a separate row and each haplotype a separate column. In a full-information setting, a 2 is assigned if the haplotype occurs two times in the individual's diplotype, a 1 if it occurs once and a 0 if it does not occur. With ambiguity, values between 0 and 2 can occur, while the rowsum is always equal to 2. 
#' 
# #' @examples
# #' lst = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' occurrence_coding(lst)
#' 
#' @export
occurrence_coding = function(lst = NULL, EM_out = NULL, EM_out_mat = NULL, haplos_keep = NULL, gene_indication = TRUE, noCNV = FALSE, CNV_option = "second"){  # 
  
  if(!is.null(lst)){
    if(is.null(EM_out)){
      EM_out = EM_algorithm(lst)
    }
    
    lst_mat = EM_out$Matrix
    
    all_cnames = strsplit(colnames(lst_mat), "\\+"); cnames = unique(unlist(all_cnames))
  }
  
  if(!is.null(EM_out_mat)){
    all_cnames = strsplit(colnames(EM_out_mat), "\\+"); cnames = unique(unlist(all_cnames))

    rm_index = unlist(lapply(all_cnames, function(x){FALSE %in% (x %in% cnames)}))  # Why was this added to this function?? It seems useless
    lst_mat = EM_out_mat[, !rm_index]
  }
  
  if(gene_indication & !is.null(haplos_keep)){
    haplos_keep = lapply(strsplit(haplos_keep, "\\_"), function(x){x[2]})
  }
  
  
  if(!is.null(haplos_keep)){
    keep = cnames[cnames %in% haplos_keep]
    remove = cnames[!(cnames %in% haplos_keep)]

    diplo_names = colnames(lst_mat)
    N = length(diplo_names)

    for(i in 1:N){
      haplos = unlist(strsplit(diplo_names[i], "\\+"))
      if(TRUE %in% (haplos %in% remove)){
        index = which(haplos %in% remove)

        haplos[index] = "COMB"
      }

      diplo_names[i] = paste(haplos, collapse = "+")
    }

    colnames(lst_mat) = diplo_names
    cnames = keep

    cnames = sort(cnames)
    len_cnames = length(cnames)
    n = nrow(lst_mat)
    
    dummy_code = matrix(0, nrow = n, ncol = len_cnames, dimnames = list(NULL, cnames))
    for(i in 1:n){
      xx = lst_mat[i, ]
      
      genos = xx[which(xx != 0)]
      names_genos = names(genos)
      len_genos = length(genos)
      
      for(j in 1:len_genos){
        splt = unlist(strsplit(names_genos[j], "\\+"))
        len_splt = length(splt)
        
        for(k in 1:len_splt){
          splttd = splt[k]
          
          if(splttd %in% cnames){
            index = which(cnames %in% splttd)
            dummy_code[i, index] = dummy_code[i, index] + genos[j]
          } else {
            splitted = unlist(strsplit(splttd, "\\^"))
            len_splitted = length(splitted)
            
            for(l in 1:len_splitted){
              index = which(cnames %in% splitted[l])
              dummy_code[i, index] = dummy_code[i, index] + genos[j]
            }
          }
        }
      }
    }
    
    
  } else {
    diplos = colnames(lst_mat); D = length(diplos); diplos_splt = strsplit(diplos, "\\+")
    haplos = sort(unique(unlist(diplos_splt))); M = length(haplos)
    
    nr_genes = length(strsplit(haplos[1], "\\-")[[1]])
    
    if(noCNV){
      haplos_splt = lapply(strsplit(haplos, "\\-"), strsplit, "\\^")
    } else {
      haplos_splt = lapply(strsplit(haplos, "\\-"), as.list)
    }
    
    HTs_in_DTs = lapply(diplos_splt, function(x){c(which(haplos == x[1]), which(haplos == x[2]))})  
    
    ##
    options = lapply(haplos_splt, function(x){apply(expand.grid(x), 1, paste, collapse = "-")})
    if(CNV_option == "first"){
      lens_options = 1; weighting = 1
      
    } else if(CNV_option == "second"){
      lens_options = lengths(options); weighting = 1
      
    } else if(CNV_option == "third"){
      lens_options = lengths(options); weighting = unlist(lapply(haplos_splt, function(x){length(unlist(x))})) / nr_genes
      
    }
    HTs = unique(unlist(options)); len_HTs = length(HTs)
    
    Repeat_mat_HTs = matrix(0, nrow = M, ncol = len_HTs, dimnames = list(haplos, HTs))
    for(j in HTs){
      beta_index_val = unlist(lapply(options, function(x){sum(x %in% j)})) / lens_options * weighting
      
      Repeat_mat_HTs[, j] = beta_index_val
    }
    
    Repeat_mat_DTs = matrix(0, nrow = D, ncol = len_HTs, dimnames = list(diplos, HTs))
    HTs_in_DTs_summed = lapply(HTs_in_DTs, function(x){colSums(Repeat_mat_HTs[x, ])})
    for(j in 1:D){
      Repeat_mat_DTs[j, ] = HTs_in_DTs_summed[[j]]
    }
    
    dummy_code = lst_mat %*% Repeat_mat_DTs

        
    # Repeat_mat = matrix(0, nrow = D, ncol = M, dimnames = list(diplos, haplos))
    # for(j in 1:M){
    #   Repeat_mat[, j] = unlist(lapply(diplos_splt, function(x){sum(x == haplos[j])}))
    # }
    # 
    # if(noCNV){
    #   nr_genes = length(strsplit(haplos[1], "\\-")[[1]])
    # 
    #   if(nr_genes == 1){
    #     haplos_splt = strsplit(haplos, "\\^")
    # 
    #   } else {
    #     haplos_splt = lapply(strsplit(haplos, "\\-"), function(x){apply(expand.grid(strsplit(x, "\\^")), 1, paste, collapse = "-")})
    # 
    #   }
    # 
    #   alleles = unique(unlist(haplos_splt)); A = length(alleles)
    # 
    #   Repeat_mat2 = matrix(0, nrow = M, ncol = A, dimnames = list(haplos, alleles))
    #   for(j in 1:A){
    #     Repeat_mat2[, j] = unlist(lapply(haplos_splt, function(x){sum(x == alleles[j])}))
    #   }
    # 
    #   Repeat_mat = Repeat_mat %*% Repeat_mat2
    # }
    # 
    # dummy_code = lst_mat %*% Repeat_mat
  }
  
  
  return(dummy_code)
}
