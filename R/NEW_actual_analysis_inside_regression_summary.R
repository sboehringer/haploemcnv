#'
#'
#'
actual_analysis_inside_regression_summary = function(last_betas, beta_genes, b0, gene_notation = TRUE, full = TRUE, combinations = NULL, mult_factor = 1){
  
  nr_genes = length(beta_genes); gene_names = names(beta_genes)
  beta_genes_full = beta_genes = lapply(beta_genes, function(x){x * mult_factor})
  for(y in 1:nr_genes){
    beta_genes_y = beta_genes_full[[y]][!str_detect(names(beta_genes_full[[y]]), "NEG")]; names_beta_genes_y = names(beta_genes_y)
    
    nr_CNs = max(lengths(strsplit(names_beta_genes_y, "\\^")))

    new_beta_genes_y = apply(expand.grid(rep(list(beta_genes_y), nr_CNs)), 1, function(x){sum(x)})
    names(new_beta_genes_y) = apply(expand.grid(rep(list(names_beta_genes_y), nr_CNs)), 1, function(x){paste(sort(unlist(strsplit(x, "\\^"))), collapse = "^")})

    new_beta_genes_y = new_beta_genes_y[!duplicated(names(new_beta_genes_y))]
    
    beta_genes_full[[y]] = c(beta_genes_full[[y]], new_beta_genes_y)
  }
  ##
  beta_genes_full = lapply(beta_genes_full, function(x){x[!duplicated(names(x))]})
  ##
  names_beta_genes_full = lapply(beta_genes_full, names)
  
  beta_genes_full_new = NULL

  
  # if(!full & !gene_notation){
  #   gene_notation = TRUE
  #   cat("if full is FALSE, gene_notation needs to be TRUE")
  # }
  
  for(y in 2:nr_genes){
    if(y == 2){
      cutted_vector = cut_vector(1:3); len_cutted_vector = length(cutted_vector)
    } else {
      cutted_vector = cut_vector(1:y, vect_list = cutted_vector); len_cutted_vector = length(cutted_vector)
    }
      
    beta_genes_full_new_sub = NULL
    for(z in 1:len_cutted_vector){
      beta_genes_full_new_sub = apply(expand.grid(beta_genes_full[cutted_vector[[z]]]), 1, sum)
      # if(gene_notation){
        names(beta_genes_full_new_sub) = apply(expand.grid(names_beta_genes_full[cutted_vector[[z]]]), 1, function(x){paste(paste("G", paste(cutted_vector[[z]], collapse = "-"), sep = ""), paste(x, collapse = "-"), sep = "_")})
      # } else {
        # names(beta_genes_full_new_sub) = apply(expand.grid(names_beta_genes_full[cutted_vector[[z]]]), 1, function(x){paste(x, collapse = "-")})
      # }
        
      beta_genes_full_new = c(beta_genes_full_new, beta_genes_full_new_sub)
    }
  }
  
    
  
  if(!full){
    beta_genes_full_new = setting_betas_haplotypes(haplotypes = names(beta_genes_full_new), beta_alleles = beta_genes_full, gene_indication = TRUE, combinations = combinations, add_gene_indication = TRUE)
  } else {
    HT_genes = str_sub(unlist(lapply(strsplit(names(beta_genes_full_new), "\\_"), function(x){x[1]})), start = 2); HT_genes_splt = strsplit(HT_genes, "\\-")
    HTs = unlist(lapply(strsplit(names(beta_genes_full_new), "\\_"), function(x){x[2]})); HTs_splt = strsplit(HTs, "\\-"); len_HTs = length(HTs)
    
    if(!is.null(combinations)){
      len_combinations = length(combinations)
      for(i in 1:len_combinations){
        combinations_i = combinations[[i]]; names_comb_i = which(gene_names %in% combinations_i[[1]]); 
        names_comb_iC = paste(names_comb_i, collapse = "-"); names_comb_iC_splt = unlist(strsplit(names_comb_iC, "\\-"))
        
        alleles_comb_i = combinations_i[[2]]; len_alleles_comb_i = length(alleles_comb_i) - 1
        adj_value_i = as.numeric(alleles_comb_i[len_alleles_comb_i + 1]); alleles_comb_i = alleles_comb_i[-(len_alleles_comb_i + 1)]
        alleles_comb_iC = paste(alleles_comb_i, collapse = "-"); alleles_comb_iC_splt = unlist(strsplit(alleles_comb_iC, "\\-"))
        
        extra_effect = 0; check_names = check_alleles = rep(TRUE, len_HTs); names_index = alleles_index = NULL
        for(j in 1:len_alleles_comb_i){
          extra_effect = extra_effect + (beta_genes[[names_comb_i[j]]][alleles_comb_i[j]] * (adj_value_i - 1))
          
          names_index = lapply(HT_genes_splt, function(x){which(x %in% names_comb_iC_splt[j])}); names_index = unlist(lapply(names_index, function(x){ifelse(length(x) == 0, 0, x)}))
          check_names = check_names & names_index != 0 #  unlist(lapply(HT_genes_splt[1:20], function(x){x[names_comb_i[j]] == names_comb_iC_splt[j]}))
          
          alleles_index = sapply(1:length(HTs_splt), function(x){HTs_splt[[x]][names_index[x]] == alleles_comb_i[j]}); alleles_index = unlist(lapply(alleles_index, function(x){ifelse(length(x) == 0, 0, x)}))
          check_alleles = check_alleles &  alleles_index != 0
        }
        
        beta_genes_full_new[check_names & check_alleles] = beta_genes_full_new[check_names & check_alleles] + extra_effect
        # betas[HT_genes == names_comb_iC & HTs == alleles_comb_iC] = extra_effect
      }
    }
  }
  
  
  
  if(gene_notation){
    for(y in 1:nr_genes){
      names(beta_genes[[y]]) = sapply(names(beta_genes[[y]]), function(x){paste("G", y, "_", x, sep = "")})
      names(beta_genes_full[[y]]) = sapply(names(beta_genes_full[[y]]), function(x){paste("G", y, "_", x, sep = "")})
    }
  } else {
    names(beta_genes_full_new) = unlist(lapply(strsplit(names(beta_genes_full_new), "\\_"), function(x){x[2]}))
  }

  
  names(beta_genes_full) = NULL
  if(!full){
    true_betas = c("Mu" = b0, unlist(beta_genes_full), beta_genes_full_new)
  } else {
    true_betas = c("Mu" = b0, beta_genes_full_new[lengths(HT_genes_splt) == nr_genes])
  }

  
  bias1 = sqrt(colMeans(t(apply(last_betas, 1, function(y){(y - true_betas[names(y)])^2})), na.rm = TRUE))
  bias2 = abs(bias1 / true_betas[names(bias1)])
  
  bias_median = abs(apply(t(apply(last_betas, 1, function(y){y - true_betas[names(y)]})), 2, function(y){median(y, na.rm = TRUE)}))
  
  last_betas_index = which(rowSums(t(apply(last_betas, 1, function(y){abs(y[names(y) != "Mu"]) > min(abs(unlist(beta_genes)[unlist(beta_genes) != 0]))})), na.rm = TRUE) != 0)  # 1.5e-8
  len_last_betas_index = length(last_betas_index)
  
  if(len_last_betas_index == 1){
    bias1_subset = (last_betas[last_betas_index, ] - true_betas[names(last_betas[last_betas_index, ])])^2
    bias1_subset = sqrt(bias1_subset[!is.na(bias1_subset)])
    
    bias2_subset = abs(bias1_subset / true_betas[names(bias1_subset)])
    
  } else if(len_last_betas_index > 1){
    bias1_subset = sqrt(colMeans(t(apply(last_betas[last_betas_index, ], 1, function(y){(y - true_betas[names(y)])^2})), na.rm = TRUE))
    
    bias2_subset = abs(bias1_subset / true_betas[names(bias1_subset)])
    
  } else {
    bias_1_subset = bias_2_subset = NULL
    
  }

  
  return(list("bias1" = bias1, "bias2" = bias2, "bias_median" = bias_median, "last_betas_index" = last_betas_index, "bias1_subset" = bias1_subset, "bias2_subset" = bias2_subset))
}
