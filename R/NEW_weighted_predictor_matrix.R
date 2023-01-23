#'
#'
#'
weighted_predictor_matrix = function(mat, type_regression = "all", noCNV = FALSE, CNV_option = "second"){
  
  N = nrow(mat)
  nr_genes = lengths(strsplit(unlist(strsplit(colnames(mat)[1], "\\+"))[1], "\\-"))
  

# Predictor matrix ----
  ind_mat = t(apply(mat, 1, function(x){x != 0}))
  ind_DTs = unlist(apply(ind_mat, 1, function(x){names(x[x == TRUE])}))
  ind_DTs_nr = apply(ind_mat, 1, function(x){sum(x)}); ind_DTs_tot = sum(ind_DTs_nr)
  
  ind_weights = as.numeric(unlist(apply(mat, 1, function(x){x[x != 0]})))
  
  pred_mat = pred_mat2 = matrix(c(rep(1:N, ind_DTs_nr), ind_DTs, ind_weights), 
                    nrow = ind_DTs_tot, ncol = 3, byrow = FALSE, dimnames = list(NULL, c("ID", "DT", "Weight")))



# Alleles and HTs ----- 
  gene_combinations = as.list(1:nr_genes)
  if(type_regression != "alleles"){
    if(nr_genes == 2){
      gene_combinations = append(gene_combinations, list(c(1, 2)))
      
    } else {
      cutted_vector = cut_vector(1:3)
      gene_combinations = append(gene_combinations, cutted_vector)
      
      for(i in 3:nr_genes){
        cutted_vector = cut_vector(1:i, vect_list = cutted_vector)
        
        gene_combinations = append(gene_combinations, cutted_vector)
      }
    }
  }
  
  
  diplos = colnames(mat); D = length(diplos); diplos_splt = strsplit(diplos, "\\+")
  
  
  haplos = sort(unique(unlist(diplos_splt))); M = length(haplos) 
  if(noCNV){
    haplos_splt = lapply(strsplit(haplos, "\\-"), strsplit, "\\^")
  } else {
    haplos_splt = lapply(strsplit(haplos, "\\-"), as.list)
  }
  
  HTs_in_DTs = lapply(diplos_splt, function(x){c(which(haplos == x[1]), which(haplos == x[2]))})  

  for(i in gene_combinations){
    options = lapply(haplos_splt, function(x){apply(expand.grid(x[i]), 1, paste, collapse = "-")})
    if(CNV_option == "first"){
      lens_options = 1; weighting = 1
      
    } else if(CNV_option == "second"){
      lens_options = lengths(options); weighting = 1
      
    } else if(CNV_option == "third"){
      lens_options = lengths(options); weighting = unlist(lapply(haplos_splt, function(x){length(unlist(x[i]))})) / length(i)
      
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
    
    pred_mat_new = Repeat_mat_DTs[ind_DTs, ]; colnames(pred_mat_new) = paste("G", paste(i, collapse = "-"), "_", colnames(pred_mat_new), sep = "")
    pred_mat = cbind(pred_mat, pred_mat_new)
  }
  
  
  # for(i in gene_combinations){
  #   for(j in 1:len_cnames_splt){
  #     splttd = strsplit(cnames_splt[[j]], "\\-")
  #     
  #     DTs[j] = paste(unlist(lapply(splttd, function(x){paste(x[i], collapse = "-")})), collapse = "+")
  #   }
  #   DTs_splt = strsplit(DTs, "\\+"); HTs = sort(unique(unlist(DTs_splt))); len_HTs = length(HTs)
  #   
  #   Repeat_mat = matrix(0, nrow = len_cnames_splt, ncol = len_HTs, dimnames = list(cnames, HTs))
  #   for(j in 1:len_HTs){
  #     Repeat_mat[, j] = unlist(lapply(DTs_splt, function(x){sum(x == HTs[j])}))
  #   }
  #   
  #   if(noCNV){
  #     nr_genes = length(strsplit(HTs[1], "\\-")[[1]])
  #     
  #     if(nr_genes == 1){
  #       HTs_splt = strsplit(HTs, "\\^")
  #     } else {
  #       HTs_splt = lapply(strsplit(HTs, "\\-"), function(x){apply(expand.grid(strsplit(x, "\\^")), 1, paste, collapse = "-")})
  #     }
  #     As = unique(unlist(HTs_splt)); len_As = length(As)
  #     
  #     Repeat_mat2 = matrix(0, nrow = len_HTs, ncol = len_As, dimnames = list(HTs, As))
  #     for(j in 1:len_As){
  #       Repeat_mat2[, j] = unlist(lapply(HTs_splt, function(x){sum(x == As[j])}))
  #     }
  #     
  #     Repeat_mat = Repeat_mat %*% Repeat_mat2
  #     HTs = sort(unique(unlist(HTs_splt))); len_HTs = length(HTs)
  #   }
  #   
  #   pred_mat_new = Repeat_mat[ind_DTs, ]; colnames(pred_mat_new) = paste("G", paste(i, collapse = "-"), "_", colnames(pred_mat_new), sep = "")
  #   pred_mat2 = cbind(pred_mat2, pred_mat_new)
  #   # 
  #   # pred_cols = c(pred_cols, rep(paste(i, collapse = "-"), len_HTs))
  # }
  
  pred_mat = data.frame(pred_mat, check.names = FALSE, row.names = NULL)
  pred_mat[, -2] = apply(pred_mat[, -2], 2, function(x){as.numeric(x)})
  
  pred_mat = pred_mat[, -2]
  
  
  return(list("pred_mat" = pred_mat, "nr_DTs" = ind_DTs_nr))
}
