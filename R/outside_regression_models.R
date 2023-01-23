#'
#'
#'
outside_regression_models = function(mat, candidate_list = NULL, outcome, outcome_true = NULL, status = NULL, gene_indication = TRUE, alpha = 0.5,
                                     regression = NULL, full_regression = NULL, type_regression = "all", regression_reference = "NEG", noCNV = TRUE, CNV_option = "second"){
  
  nr_genes = lengths(strsplit(unlist(strsplit(colnames(mat)[1], "\\+"))[1], "\\-"))

  new_occ_mat_other_genes = new_occ_mat = EPE = NULL
  

  if(full_regression){
    mod_out_res = list("Coefs" = rep(list(NULL), 1))  # nr_genes  
    for(j in 1:nr_genes){
      new_occ_mat = occurrence_coding(EM_out_mat = reducing_mat(EM_mat = mat, nr_gene = j), gene_indication = gene_indication, haplos_keep = candidate_list)
      colnames(new_occ_mat) = paste("G", j, "_", colnames(new_occ_mat), sep = "")
      
      new_occ_mat_other_genes = cbind(new_occ_mat_other_genes, new_occ_mat)
    }
    
    
    ##
    gene_combinations = NULL
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
    len_gene_combinations = length(gene_combinations)
    
    
    occ_mat = new_occ_mat = NULL
    for(j in gene_combinations){
      new_occ_mat = occurrence_coding(EM_out_mat = reducing_mat(EM_mat = mat, nr_gene = j), haplos_keep = candidate_list)
      
      if(ncol(new_occ_mat) != 0){
        colnames(new_occ_mat) = paste("G", paste(j, collapse = "-"), "_", colnames(new_occ_mat), sep = "")
        
        occ_mat = cbind(occ_mat, new_occ_mat)
      }
    }
    
    pred_mat = cbind(new_occ_mat_other_genes, occ_mat)
    pred_mat = pred_mat[, !str_detect(colnames(pred_mat), regression_reference)]
    
    
    ##
    if(TRUE %in% apply(pred_mat, 2, function(x){sd(x) == 0})){
      pred_mat = pred_mat[, !apply(pred_mat, 2, function(x){sd(x) == 0})]
    }
    ##
    # mod_out = glm(outcome ~ pred_mat, family = regression); sum_mod_out = summary(mod_out)
    mod_out = cv.glmnet(y = outcome, x = pred_mat, family = regression, alpha = alpha, nlambda = 200)
    
    mod_out_lambda = mod_out$lambda.1se; mod_out_index = mod_out$index["1se", "Lambda"]
    mod_out_coef = coef(mod_out, s = mod_out_lambda)[, 1]; mod_out_coef_names = names(mod_out_coef); mod_out_coef_names[1] = "Intercept"
    
    fitted.values = predict(mod_out, newx = pred_mat, s = mod_out_lambda, type = "response")
    
    if(regression == "gaussian"){
      # EPE = mean((outcome_true - mod_out$fitted.values)^2)
      EPE = mean((outcome_true - fitted.values)^2)
      
    } else if(regression == "binomial"){
      # EPE = mean((outcome - mod_out$fitted.values)^2)
      EPE = mean((outcome - fitted.values)^2)
      
    }
    
    mod_out_res$Coefs[[1]] = matrix(mod_out_coef, nrow = length(mod_out_coef), ncol = 1, dimnames = list(paste("pred_mat", names(mod_out_coef), sep = ""), "Estimate"))  # sum_mod_out$coef
    
    
    
    ##
    # 
    # for(j in 1:nr_genes){
    #   if(j == 1){
    #     pred_mat = new_occ_mat_other_genes
    #     
    #   } else{
    #     if(j < nr_genes){
    #       red_mat = reducing_mat(EM_mat = mat, nr_gene = c(1:j))
    #     } else {
    #       red_mat = mat
    #     }
    #     new_occ_mat = occurrence_coding(EM_out_mat = red_mat, haplos_keep = candidate_list)
    #     
    #     if(ncol(new_occ_mat) != 0){
    #       colnames(new_occ_mat) = paste("G", paste(1:j, collapse = "-"), "_", colnames(new_occ_mat), sep = "")
    #       
    #       pred_mat = cbind(new_occ_mat_other_genes, new_occ_mat)
    #     }
    #     
    #   } 
    #   pred_mat = pred_mat[, !str_detect(colnames(pred_mat), regression_reference)]
    #   
    #   
    #   mod_out = glm(outcome ~ pred_mat, family = regression); sum_mod_out = summary(mod_out)
    #   
    #   if(j == nr_genes){
    #     if(regression == "gaussian"){
    #       EPE = mean((outcome_true - mod_out$fitted.values)^2)
    #       
    #     } else if(regression == "binomial"){
    #       EPE = mean((outcome - mod_out$fitted.values)^2)
    #       
    #     }
    #   }
    #   
    #   mod_out_res$Coefs[[j]] = sum_mod_out$coef
    #   # mod_out_res$adjR2[j] = sum_mod_out$adj.r.squared
    #   
    #   new_occ_mat_other_genes = pred_mat
    # }
    
  } else {
    mod_out_res = list("Coefs" = rep(list(NULL), 1))  
    for(j in 1:nr_genes){
      new_occ_mat = occurrence_coding(EM_out_mat = reducing_mat(EM_mat = mat, nr_gene = j), noCNV = noCNV, CNV_option = CNV_option)
      colnames(new_occ_mat) = paste("G", j, "_", colnames(new_occ_mat), sep = "")
      
      new_occ_mat_other_genes = cbind(new_occ_mat_other_genes, new_occ_mat)
    }
    
    pred_mat = new_occ_mat_other_genes
    if(!is.null(type_regression) && type_regression == "backward"){
      pred_mat = cbind(pred_mat, iterative_modeling(EM_mat = mat, outcome = outcome, regression = regression, noCNV = noCNV)$pred_mat)
    }
    pred_mat = pred_mat[, !str_detect(colnames(pred_mat), regression_reference)]
    pred_mat = pred_mat[, !str_detect(colnames(pred_mat), "Comb")]
    
    
    mod_out = glm(outcome ~ pred_mat, family = regression); sum_mod_out = summary(mod_out)
    
    if(regression == "gaussian"){
      EPE = mean((outcome_true - mod_out$fitted.values)^2)
        
    } else if(regression == "binomial"){
      EPE = mean((outcome - mod_out$fitted.values)^2)
        
    }
    
    mod_out_res$Coefs[[1]] = sum_mod_out$coef
    
    
    
    # coefs = sum_mod_out$coef[, "Estimate"][-1]; names_coefs = names(coefs)
    # 
    # un_HTs = strsplit(unique(unlist(strsplit(colnames(mat), "\\+"))), "\\-")
    # un_HTs_gene = rep(list(NULL), nr_genes)
    # for(i in 1:nr_genes){
    #   un_HTs_gene_i = unique(unlist(lapply(un_HTs, function(x){x[i]})))
    #   
    #   un_HTs_gene_i = un_HTs_gene_i[!str_detect(un_HTs_gene_i, regression_reference)]
    #   un_HTs_gene_i = un_HTs_gene_i[!str_detect(un_HTs_gene_i, "Comb")]
    #   
    #   un_HTs_gene[[i]] = paste(paste("G", i, sep = ""), un_HTs_gene_i, sep = "_") 
    # }
    # 
    # un_HTs = lapply(un_HTs_gene, function(x){unlist(lapply(strsplit(x, "_"), function(y){y[2]}))})
    # un_HTs_splt = lapply(un_HTs, function(x){strsplit(x, "\\^")})
    # 
    # coefs_genes = as.numeric(unlist(lapply(strsplit(str_sub(names_coefs, start = 10), "\\_"), function(x){x[1]})))
    # coefs_alleles = unlist(lapply(strsplit(str_sub(names_coefs, start = 10), "\\_"), function(x){x[2]}))
    # 
    # coefs_recal = rep(list(NULL), nr_genes)
    # for(i in 1:nr_genes){
    #   un_HTs_splt_i = un_HTs_splt[[i]]; un_HTs_i = un_HTs[[i]]; len_un_HTs_i = length(un_HTs_i)
    #   
    #   coefs_recal[[i]] = rep(0, len_un_HTs_i); names(coefs_recal[[i]]) = un_HTs_i
    #   if(CNV_option == "first"){
    #     weighting = rep(1, len_un_HTs_i)
    #       
    #   } else if(CNV_option == "second"){
    #     weighting = 1 / lengths(un_HTs_splt_i)
    #       
    #   } else if(CNV_option == "third"){
    #     weighting = unlist(lapply(un_HTs_splt_i, function(x){(length(unlist(x)) / 1)})) / lengths(un_HTs_splt_i)
    #       
    #   }
    #   
    #   coefs_genes_i = coefs_genes == i
    #   for(j in 1:len_un_HTs_i){
    #     beta_index_val = sapply(coefs_alleles, function(x){sum(x == un_HTs_splt_i[[j]])})
    #     beta_index = coefs_genes_i & ifelse(beta_index_val != 0, TRUE, FALSE)  
    #     
    #     coefs_recal[[i]][j] = sum((coefs * beta_index_val)[beta_index] * weighting[j])
    #   }
    #   
    #   names(coefs_recal[[i]]) = paste(paste("G", i, sep = ""), names(coefs_recal[[i]]), sep = "_")
    # }
    # 
    # intercept = sum_mod_out$coef[, "Estimate"][1]; names(intercept) = "Intercept"
    # mod_out_res$Coefs[[1]] = c(intercept, unlist(coefs_recal))
  }
  
  
  out = NULL
  
  out$EPE = EPE
  out$mod_out_res = mod_out_res
    
  return(out)
}
