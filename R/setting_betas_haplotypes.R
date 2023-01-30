#' @title Setting betas for all haplotypes based individual allele effect sizes
#' 
#' @param haplotypes A vector. Haplotypes and their frequencies.
#' @param beta_alleles An list. Overview of the effect sizes of each allele, each list element contains the alleles of a different gene. Can be obtained via \code{KIR_simulation_preparation}.
#' @param gene_indication A logical scalar. Whether or not it is indicated from which gene combination the allele or haplotype is (default is \emph{TRUE}).
#' @param add_gene_indication A logical scalar. Whether or not the indication from which gene combination the allele or haplotype is should be included (default is \emph{TRUE}).
# #' @param adj_scale A numeric value in \{0, 1\}. The degree that the different sizes of haplotypes contribute to the simulated outcome. If \code{adj_scale} is 0, the simulated outcome is fully based on the one-gene haplotypes. In contrast, if \code{adj_scale} is 1, the simulated outcome is fully based on the full-gene haplotypes. 
#' @param combinations An optional list. Which allele combinations have full effect sizes, and what the loss in effect is of alleles not in this combination. Allows for the addition of haplotype effects. # Note that if \code{adj_scale} is 0 this list is ignored.
# #' @param penetrance A numeric value in \{0, 1\}. How much of the effect all haplotypes not in the combination should be of their \emph{actual} effect.
#'
# #' @description Support function of \code{\link{KIR_simulation_one_subiter}}.
#'
#' @return A vector with the estimates beta for each haplotype combination
#' 
# #' @examples
# #' haplotypes = c("A-E", "A-F", "A-NEG", "B-E", "B-F", "B-NEG", "C-E", "C-F", "C-NEG")
# #' betas_genes = list(c("A" = 0.2, "B" = 0, "C" = 0, "NEG" = 0), c("E" = 0, "F" = -0.1, "NEG" = 0))
# #' 
# #' setting_betas_haplotypes(haplotypes, betas_genes)
#'
#' @export 
setting_betas_haplotypes = function(haplotypes, beta_alleles, gene_indication = TRUE, add_gene_indication = FALSE, combinations = NULL){
  
  N = length(haplotypes)
  # if(gene_indication){
  #   haplotypes = unlist(lapply(strsplit(haplotypes, "\\_"), function(x){x[2]}))
  # }
  HT_genes = str_sub(unlist(lapply(strsplit(haplotypes, "\\_"), function(x){x[1]})), start = 2)
  HTs = unlist(lapply(strsplit(haplotypes, "\\_"), function(x){x[2]}))
  
  
  gene_names = names(beta_alleles)
  # if(is.null(nr_genes)){  # For adj_scale it is important to know the number of genes in the full haplotype set (and not only this subset of haplotypes)
  #   nr_genes = length(beta_alleles) 
  # } 
  nr_genes = length(beta_alleles)
  betas = nr_gene_copies = vector(mode = "numeric", length = N); names(betas) = names(nr_gene_copies) = HTs
  
  
  ##
  for(i in 1:nr_genes){
    beta_alleles_i = beta_alleles[[i]]; len_beta_alleles_i = length(beta_alleles_i)
    betas_index = (HT_genes == i)
    
    for(j in 1:len_beta_alleles_i){
      beta_alleles_ij = beta_alleles_i[j]; names_beta_alleles_ij = names(beta_alleles_ij)
      
      if(!is.na(betas[names_beta_alleles_ij])){
        betas[names_beta_alleles_ij] = beta_alleles_ij
      }
    }
  }
  
  
  if(!is.null(combinations)){
    len_combinations = length(combinations)
    for(i in 1:len_combinations){
      combinations_i = combinations[[i]]
      
      
      names_comb_i = which(gene_names %in% combinations_i[[1]]); names_comb_iC = paste(names_comb_i, collapse = "-")
      
      alleles_comb_i = combinations_i[[2]]; len_alleles_comb_i = length(alleles_comb_i) - 1
      adj_value_i = as.numeric(alleles_comb_i[len_alleles_comb_i + 1]); alleles_comb_i = alleles_comb_i[-(len_alleles_comb_i + 1)]; alleles_comb_iC = paste(alleles_comb_i, collapse = "-")
      
      extra_effect = 0
      for(j in 1:len_alleles_comb_i){
        extra_effect = extra_effect + (beta_alleles[[names_comb_i[j]]][alleles_comb_i[j]] * (adj_value_i - 1))
      }
      
      betas[HT_genes == names_comb_iC & HTs == alleles_comb_iC] = extra_effect
    }
  }
  
  
  if(add_gene_indication){
    names(betas) = haplotypes
  }
  ##
  

  
  
  # if(is.null(combinations)){
  #   for(i in 1:N){
  #     splt = unlist(HT_splt[[i]]); nr_gene_copies[i] = length(splt)
  #     
  #     for(j in 1:nr_gene_copies[i]){
  #       betas_ij = beta_alleles[[j]][splt[j]]
  #       
  #       if(!is.na(betas_ij)){
  #         betas[i] = betas[i] + betas_ij
  #       }
  #     }
  #   }
  #   
  # } else {
  #   len_combinations = length(combinations)
  #   adj_value = as.numeric(unlist(lapply(combinations, function(x){x[2][[1]][length(x[2][[1]])]})))
  #   
  #   for(i in 1:N){
  #     splt = unlist(HT_splt[[i]]); nr_gene_copies[i] = length(splt)
  #     
  #     betas_splt = NULL
  #     for(j in 1:nr_gene_copies[i]){
  #       betas_splt[j] = beta_alleles[[j]][splt[j]]
  #     }
  #     
  #     eligible_i = rep(FALSE, len_combinations)
  #     for(j in 1:len_combinations){
  #       combo = combinations[[j]]; len_combo = length(combo[[1]])
  #       
  #       if(len_combo > nr_gene_copies[i]){
  #         next
  #       }
  #       
  #       eligible_j = rep(FALSE, len_combo)
  #       for(k in 1:len_combo){
  #         gene_index = which(gene_names == combo[[1]][k])
  #         
  #         if(length(gene_index) > 0 && nr_gene_copies[i] > gene_index && HT_splt[[i]][gene_index] == combinations[[j]][[2]][gene_index]){
  #           eligible_j[k] = TRUE
  #         }
  #       }
  #       
  #       if(!(FALSE %in% eligible_j)){
  #         eligible_i[j] = TRUE
  #       }
  #     }
  #     
  #     
  #     if(TRUE %in% eligible_i){
  #       if(sum(eligible_i) > 1){
  #         cat("For haplotype ", HTs[i], " there are multiple combinations... code does not allow that... \n")
  #       }
  #       combo = combinations[[which(eligible_i)]]; len_combo = length(combo[[1]])
  #       
  #       gene_indices = NULL
  #       for(k in 1:len_combo){
  #         gene_indices = c(gene_indices, which(gene_names == combo[[1]][k]))
  #       }
  #       
  #       if(length(gene_indices) == nr_gene_copies[i]){
  #         betas[i] = sum(betas_splt[gene_indices]) * adj_value[eligible_i]
  #       } else {
  #         betas[i] = sum(betas_splt[gene_indices]) * adj_value[eligible_i] + betas_splt[-gene_indices] 
  #       }
  #       
  #     } else {
  #       betas[i] = sum(betas_splt)
  #     }
  #   }
  # }
  # 
  # if(add_gene_indication){
  #   names(betas) = sapply(names(betas), function(x){paste(paste("G", paste(1:lengths(strsplit(x, "\\-")), collapse = "-"), sep = ""), x, sep = "_")})
  # }
  # 
  # 
  # 
  # if(adj_scale != 0.5){
  # 
  #   dbeta_norm = function(x, s1, s2){
  #     db_max = dbeta(0.5, s1, s2)
  #     db = dbeta(x, s1, s2)
  # 
  #     return(db / db_max)
  #   }
  # 
  # 
  #   if(adj_scale == 0){
  #     betas[nr_gene_copies != 1] = 0
  # 
  #   } else if (adj_scale == 1){
  #     betas[nr_gene_copies != nr_genes] = 0
  # 
  #   } else {
  #     warning("adj_scale values that are not 0 or 1 do strange things...\n")
  # 
  #     nr_adj = max(1, (nr_genes - 1))
  #     if(as.integer(nr_adj) %% 2 == 0){  # I think this is to choose dbeta form
  #       adj_range = c((3:nr_genes)[1:(nr_adj / 2)], (3:nr_genes)[1:(nr_adj / 2)])
  #     } else {
  #       adj_range = c((3:nr_genes)[1:floor(nr_adj / 2)], (2:nr_genes)[1:ceiling(nr_adj / 2)])
  #     }
  # 
  #     if(adj_scale > 0.5){
  #       temp_adj_scale = abs(adj_scale - 1)
  #     } else {
  #       temp_adj_scale = adj_scale
  #     }
  # 
  # 
  #     adj_values = 1
  #     mirroring = TRUE
  #     for(i in 1:nr_adj){
  #       if(mirroring == TRUE){
  #         adj_values[i + 1] = -dbeta_norm(0.5 - temp_adj_scale, adj_range[i], 1) + 1
  # 
  #       } else {
  #         adj_values[i + 1] = dbeta_norm(temp_adj_scale, adj_range[i], 1)
  # 
  #       }
  # 
  #       if(adj_range[i] == 3){  # Why is this 3?
  #         mirroring = FALSE
  #       }
  #     }
  # 
  #     if(adj_scale > 0.5){
  #       adj_values = rev(adj_values)
  #     }
  # 
  #     for(i in 1:nr_genes){
  #       betas[nr_gene_copies == i] = betas[nr_gene_copies == i] * adj_values[i]
  #     }
  #   }
  # } else {
  #   warning("adj_scale values that are not 0 or 1 do strange things...\n")
  # }
  
  return(betas)
}
