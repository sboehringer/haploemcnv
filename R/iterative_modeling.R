#' @title Running regression models with only significant haplotypes
#' 
#' @description Regression models with only the allele effects are appropriate when the effects are only attributed to alleles, but additional effects that occur when specific alleles come together in a haplotype. Therefore certain haplotypes will need to be included as predictors in the regression models. Because of the abundance of haplotypes, only a limited number of haplotypes are allowed. Via backward selection the haplotypes with the lowest p-values are retained. Multi-loci haplotypes are only considered if one of its subset haplotypes are retained in the predictor matrix.
#' 
#' @param yi An optional vector. The outcome of each individual present in \code{lst}. Only functional when \code{regression} is not \code{NULL}.
#' @param EM_mat An optional matrix. Denoting for each individual the probabilities estimated for each possible diplotypes, as can be obtained via \code{EM-algorithm}.
#' @param regression An optional character. Multiple options are possible: \code{NULL} (default), \code{gaussian} and \code{logistic}. If not \code{NULL} forces the EM-algorithm to run a regression analysis each iteration which is used for haplotype frequency estimation. Then requires an outcome for each individual in \code{outcome}.  If \code{NULL}, haplotype frequency estimation will proceed without the regression analyses.
#' @param regression_reference A character. Which predictor, or which predictors, should be used as reference group in the regression models. 
#' @param only_sign_EES A logical scalar. Whether or not multi-locus haplotypes are only allowed to be constructed from significant alleles. 
#' @param pval_thresh A numeric value \in \{0, 1\}. Multi-locus haplotypes with p-values below this threshold value are retained in the regression model as predictors.
#' @param candidate_list A list. Overview of haplotypes that need to be predictors in the regression models.
#' @param noCNV a logical scalar. Whether or not the copy number variation haplotypes should be accomodated in the allelic predictors or should be kept as separate predictors. Only functional when \code{regression} is not \code{NULL}.
#' 
#' @return A list where each list element contains the estimated betas from each regression step.
#' 
iterative_modeling = iterative_modelling = function(outcome, EM_mat = NULL, weighted_mat = NULL, regression = "gaussian", regression_reference = "NEG", only_sign_EES = TRUE, pval_thresh = 0.10, 
                                                    candidate_list = NULL, noCNV = TRUE, CNV_option = "second", gene_indication = FALSE){
  
  
  
  ##
  suppress_warnings <- function(.expr, .f, ...){
    eval.parent(substitute(
      withCallingHandlers(.expr, warning = function(w){
        cm <- conditionMessage(w)
        cond <- if(is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm,...)
        if(cond){
          invokeRestart("muffleWarning")
        }})))}         
  ##
  
  
  
  
  
  sum_mod_one = list(NULL, NULL); names(sum_mod_one) = c("betas", "pvals")
  sum_mod = list(sum_mod_one)
  
  if(is.null(weighted_mat)){
    N = nrow(EM_mat)
    DTs = colnames(EM_mat); HTs = unique(unlist(strsplit(DTs, "\\+"))); Alleles = unique(unlist(strsplit(HTs, "\\-")))
    nr_genes = length(unlist(strsplit(HTs[1], "\\-")))
    
    
    ## Analysis with only alleles
    occ_mat = occ_mat_new = NULL
    for(i in 1:nr_genes){
      red_mat = reducing_mat(EM_mat = EM_mat, nr_gene = i)
      
      occ_mat_new = occurrence_coding(EM_out_mat = red_mat, haplos_keep = candidate_list, noCNV = noCNV, CNV_option = CNV_option, gene_indication = gene_indication); colnames(occ_mat_new) = paste("G", i, "_", colnames(occ_mat_new), sep = "")
      occ_mat = cbind(occ_mat, occ_mat_new)
    }
    pred_mat = pred_mat_alleles = occ_mat[, !str_detect(colnames(occ_mat), regression_reference)]
    # if(regression == "binomial"){
    pred_mat = pred_mat_alleles = pred_mat[, !str_detect(colnames(pred_mat), "Comb")]
    # } else {
    #   comb_index = which(str_detect(colnames(pred_mat), "Comb"))
    #   pred_mat = pred_mat[, -comb_index[duplicated(unlist(lapply(strsplit(colnames(pred_mat)[comb_index], "\\_"), function(x){x[2]})))]]
    # }
    
    mod_out = glm(outcome ~ pred_mat, family = regression, control = list(maxit = 500))
    
    sum_mod_out = summary(mod_out); sum_mod_out_table = sum_mod_out$coefficients
    sum_mod[[1]]$betas = sum_mod_out_table[, "Estimate"]; sum_mod[[1]]$pvals = sum_mod_out_table[, str_detect(colnames(sum_mod_out_table), "Pr")]  # sum_mod_out_table[, "Pr(>|t|)"]
    
    
    
  } else {
    ID = weighted_mat$ID; N = length(unique(ID)); weights = weighted_mat$Weight; weight_mat = as.matrix(weighted_mat[, !(colnames(weighted_mat) %in% c("ID", "Weight"))])
    
    cnames = colnames(weight_mat)
    
    cnames_which_genes = unlist(lapply(strsplit(cnames, "\\_"), function(x){str_sub(x[1], start = 2)})); nr_genes = max(as.numeric(unique(unlist(strsplit(cnames_which_genes, "\\-")))))
    cnames_which_HTs = unlist(lapply(strsplit(cnames, "\\_"), function(x){x[2]}))
    
    cnames_nr_genes = lengths(strsplit(unlist(lapply(strsplit(cnames, "\\_"), function(x){x[2]})), "\\-"))
    pred_mat = weight_mat[, cnames_nr_genes == 1]
    
    pred_mat = pred_mat_alleles = pred_mat[, !str_detect(colnames(pred_mat), regression_reference)]
    pred_mat = pred_mat_alleles = pred_mat[, !str_detect(colnames(pred_mat), "Comb")]

    mod_out = suppress_warnings(glm(outcome ~ pred_mat, family = regression, weights = weights, control = list(maxit = 500)), "non-integer #successes in a binomial glm!")
    sum_mod_out = summary(mod_out); sum_mod_out_table = sum_mod_out$coefficients
    sum_mod[[1]]$betas = sum_mod_out_table[, "Estimate"]; sum_mod[[1]]$pvals = sum_mod_out_table[, str_detect(colnames(sum_mod_out_table), "Pr")]  # sum_mod_out_table[, "Pr(>|t|)"]
  }
  
  previous_mod_out = mod_out
  
  
  new_cols = NULL; warning_stop = FALSE
  if(TRUE %in% (sum_mod[[1]]$pvals[-1] < pval_thresh)){
    ## Analysis with also HTs
    for(i in 2:nr_genes){
      
      sum_mod[[i]] = sum_mod_one
      if(warning_stop){
        sum_mod[[i]] = sum_mod[[i - 1]]
        next
      }
      
      pvals = sum_mod[[i - 1]]$pvals[!str_detect(names(sum_mod[[i - 1]]$pvals), "Intercept")]
        
      if(only_sign_EES == TRUE){
        pvals = pvals[pvals < pval_thresh]
      }
        
      pvals_splt = NULL
      pvals_splt = lapply(strsplit(names(pvals), "pred_matG"), function(x){unlist(strsplit(x[2], "\\_"))})
        
      if(i == 2){
        pvals_lst = rep(list(NULL), nr_genes)
        for(j in 1:length(pvals_splt)){
          pvals_lst[[as.numeric(pvals_splt[[j]][1])]] = sort(c(pvals_lst[[as.numeric(pvals_splt[[j]][1])]], pvals_splt[[j]][2])) 
        }
        
        cutted_vector = cut_vector(1:3); len_cutted_vector = length(cutted_vector)
        HTs_4_mod = list(NULL); occ_mats = red_mats = NULL
        
        if(is.null(weighted_mat)){
          for(j in 1:len_cutted_vector){
            HTs_4_mod[[j]] = apply(expand.grid(pvals_lst[cutted_vector[[j]]]), 1, function(x){paste(x, collapse = "-")})
            
            if(length(HTs_4_mod[[j]]) != 0){
              red_mats = reducing_mat(EM_mat = EM_mat, nr_gene = cutted_vector[[j]])
              
              occ_mats_new = occurrence_coding(EM_out_mat = red_mats, haplos_keep = HTs_4_mod[[j]], gene_indication = FALSE)
              HTs_4_mod[[j]] = colnames(occ_mats_new)
              
              colnames(occ_mats_new) = paste("G", paste(cutted_vector[[j]], collapse = "-"), "_", colnames(occ_mats_new), sep = "")
              occ_mats = cbind(occ_mats, occ_mats_new)
            }
          }
          
          
        } else {
          for(j in 1:len_cutted_vector){
            HTs_4_mod[[j]] = apply(expand.grid(pvals_lst[cutted_vector[[j]]]), 1, function(x){paste(x, collapse = "-")})
            occ_mats_index = paste(cutted_vector[[j]], collapse = "-") == cnames_which_genes & cnames_which_HTs %in% HTs_4_mod[[j]]
            
            if(sum(occ_mats_index) != 0){
              occ_mats_new = as.matrix(weight_mat[, occ_mats_index]); colnames(occ_mats_new) = colnames(weight_mat)[occ_mats_index]
              occ_mats = cbind(occ_mats, occ_mats_new)
            }
          }
        }
        HTs_4_mod_all = unlist(HTs_4_mod)
        
        
      } else {
        # pvals_splt_index = which(is.na(pvals_splt)); len_pvals_splt_index = length(pvals_splt_index)
        
        pvals_splt_index = which(unlist(lapply(pvals_splt, function(x){str_detect(x[1], "\\-")}))); len_pvals_splt_index = length(pvals_splt_index)
        if(len_pvals_splt_index == 0){
          next
          
        } else {
          cutted_vector_OLD = cutted_vector; HTs_4_mod_OLD = HTs_4_mod
          
          cutted_vector = cut_vector(vect = 1:i, vect_list = cutted_vector_OLD); un_cutted_vector = unlist(lapply(cutted_vector, function(x){paste(x, collapse = "-")})); len_cutted_vector = length(cutted_vector)        
          HTs_4_mod = list(NULL); occ_mats = red_mats = NULL
          
          for(j in 1:len_pvals_splt_index){
            # pvals_splt[[pvals_splt_index[j]]] = c(paste(unlist(cutted_vector_OLD[unlist(lapply(HTs_4_mod_OLD, function(x){HTs_4_mod_all[j] %in% x}))]), collapse = "-"), HTs_4_mod_all[j])
            
            for(k in 1:len_cutted_vector){
              pvals_splt_NEW_k_index = cutted_vector[[k]][!(cutted_vector[[k]] %in% as.numeric(unlist(strsplit(pvals_splt[[pvals_splt_index[j]]][1], "\\-"))))]
              pvals_splt_opts = which(unlist(lapply(pvals_splt, function(x){x[1] == pvals_splt_NEW_k_index}))); len_pvals_splt_opts = length(pvals_splt_opts)
              
              if(len_pvals_splt_opts == 0){
                next
              }
              
              new_pred_k = unlist(strsplit(pvals_splt[[pvals_splt_index[j]]][2], "\\-")); names(new_pred_k) = unlist(strsplit(pvals_splt[[pvals_splt_index[j]]][1], "\\-"))
              for(l in 1:len_pvals_splt_opts){
                new_pred_k_extra = pvals_splt[[pvals_splt_opts[l]]][2]; names(new_pred_k_extra) = pvals_splt_NEW_k_index
                
                # new_pred_kl = c(new_pred_k, new_pred_k_extra); new_pred_kl = paste("G", paste(cutted_vector[[k]], collapse = "-"), "_", paste(new_pred_kl[order(names(new_pred_kl))], collapse = "-"), sep = "")
                # HTs_4_mod[[k]] = unique(c(HTs_4_mod[[k]], new_pred_kl))
                
                new_pred_kl = c(new_pred_k, new_pred_k_extra); new_pred_kl = paste(new_pred_kl[order(names(new_pred_kl))], collapse = "-")
                HTs_4_mod[[k]] = unique(c(HTs_4_mod[[k]], new_pred_kl))
                
              }
            }
          }
          
          ##
          no_allele_index = which((lengths(lapply(pvals_splt, function(x){unlist(strsplit(x, "\\-"))})) / 2) == (i - 1)); len_no_allele_index = length(no_allele_index)
          HTs_4_mod_new = HTs_4_mod_extra = NULL
          
          if(len_no_allele_index > 1){
            splt1_pvals_splt = lapply(pvals_splt[no_allele_index], function(x){unlist(strsplit(x[1], "\\-"))})
            splt2_pvals_splt = lapply(pvals_splt[no_allele_index], function(x){unlist(strsplit(x[2], "\\-"))})
            for(j in 1:(len_no_allele_index - 1)){
              for(k in j:len_no_allele_index){
                no_allele_index_match = which(!(splt1_pvals_splt[[k]] %in% splt1_pvals_splt[[j]]))
                
                if(length(no_allele_index_match) == 1){
                  no_allele_index_match_order = order(c(splt1_pvals_splt[[j]], splt1_pvals_splt[[k]][no_allele_index_match]))
                  
                  HTs_4_mod_new = paste(c(splt2_pvals_splt[[j]], splt2_pvals_splt[[k]][no_allele_index_match])[no_allele_index_match_order], collapse = "-")
                  names(HTs_4_mod_new) = paste(no_allele_index_match_order[no_allele_index_match_order], collapse = "-")
                  
                  HTs_4_mod_extra = c(HTs_4_mod_extra, HTs_4_mod_new)
                }
              }
            }
            
            if(!is.null(HTs_4_mod_extra)){
              names_HTs_4_mod_extra = names(HTs_4_mod_extra); len_HTs_4_mod_extra = length(HTs_4_mod_extra)
              for(j in 1:len_HTs_4_mod_extra){
                no_allele_index = which(names_HTs_4_mod_extra[j] %in% un_cutted_vector)
                HTs_4_mod[[no_allele_index]] = unique(c(HTs_4_mod[[no_allele_index]], HTs_4_mod_extra[j]))
              }
            }
          }
          ##
          
          
          occ_mats = NULL
          if(is.null(weighted_mat)){
            for(k in 1:len_cutted_vector){
              if(length(HTs_4_mod[[k]]) != 0){
                
                if(all.equal(1:nr_genes, cutted_vector[[k]])){
                  red_mats = EM_mat
                } else {
                  red_mats = reducing_mat(EM_mat = EM_mat, nr_gene = cutted_vector[[k]])
                }
                
                occ_mats_new = occurrence_coding(EM_out_mat = red_mats, haplos_keep = HTs_4_mod[[k]], gene_indication = FALSE)
                
                ##          
                HTs_4_mod[[k]] = colnames(occ_mats_new)
                ##
                
                colnames(occ_mats_new) = paste("G", paste(cutted_vector[[k]], collapse = "-"), "_", colnames(occ_mats_new), sep = "")
                occ_mats = cbind(occ_mats, occ_mats_new)
                
              } 
            }
            
          } else {
            for(k in 1:len_cutted_vector){
              occ_mats_index = paste(cutted_vector[[k]], collapse = "-") == cnames_which_genes & cnames_which_HTs %in% HTs_4_mod[[k]]
              
              if(sum(occ_mats_index) != 0){
                occ_mats_new = as.matrix(weight_mat[, occ_mats_index]); colnames(occ_mats_new) = colnames(weight_mat)[occ_mats_index]
                occ_mats = cbind(occ_mats, occ_mats_new)
              }
            }
          }
          
          HTs_4_mod_all = unlist(HTs_4_mod)
        }
      }
      
      
      if(length(HTs_4_mod_all) == 0){
        sum_mod[[i]] = sum_mod[[i - 1]]
        
        warning_stop = TRUE
        next
      }
      
        
      if(!is.null(HTs_4_mod_all)){
        ### HTs_4_mod_all contains all theoretical HTs, but not all HTs are present in the dataset (e.g., collapsing), so need to accomodate that
        if(FALSE %in% (HTs_4_mod_all %in% unlist(lapply(strsplit(colnames(occ_mats), "\\_"), function(x){x[2]})))){
          missing_HTs_4_mod_all = HTs_4_mod_all[!(HTs_4_mod_all %in% unlist(lapply(strsplit(colnames(occ_mats), "\\_"), function(x){x[2]})))]
          
          # HTs_4_mod_all = HTs_4_mod_all[!(HTs_4_mod_all %in% missing_HTs_4_mod_all)]
          HTs_4_mod_all = unlist(lapply(strsplit(colnames(occ_mats), "\\_"), function(x){x[2]}))
          
          cat("HTs_4_mod_all and the column names of occ_mats are different. Some HTs are not present in occ_mats:", missing_HTs_4_mod_all, "so they have been discarded\n")
          
          # stop("HTs_4_mod_all and the column names of occ_mats are different...")
        } else {
          HTs_4_mod_all = unlist(lapply(strsplit(colnames(occ_mats), "\\_"), function(x){x[2]}))
        }
        ###
        
      
      
      
        cnames_start = colnames(occ_mats)
        repeat_model = TRUE
        
        cat("Dropped: ")
        while(repeat_model & !warning_stop){
  
          if(ncol(occ_mats) == 0){
            pred_mat = pred_mat_alleles
          } else {
            pred_mat = cbind(pred_mat_alleles, occ_mats)
          }
          
          
          previous_mod_out = mod_out
          if(is.null(weighted_mat)){
            # mod_out = glm(outcome ~ pred_mat, family = regression, control = list(maxit = 500))
            mod_out = tryCatch(glm(outcome ~ pred_mat, family = regression, control = list(maxit = 500)), warning = function(w){TRUE})
          } else {
            # mod_out = glm(outcome ~ pred_mat, family = regression, weights = weights, control = list(maxit = 500))
            # mod_out = tryCatch(glm(outcome ~ pred_mat, family = regression, weights = weights, control = list(maxit = 500)), warning = function(w){TRUE})
            mod_out = tryCatch(suppress_warnings(glm(outcome ~ pred_mat, family = regression, weights = weights, control = list(maxit = 500)), "non-integer #successes in a binomial glm!"),
                               warning = function(w){TRUE})
          }
  
          if(is.logical(mod_out) && mod_out == TRUE){
            cat("A warning occurred in the glm model, we fall back to iteration ", i - 1, "\n")
            warning_stop = TRUE
            next
          }
          
          sum_mod_out = summary(mod_out); sum_mod_out_table = sum_mod_out$coefficients
          
          mod_out_coefs = sum_mod_out_table[, "Estimate"]; mod_out_pvals = sum_mod_out_table[, str_detect(colnames(sum_mod_out_table), "Pr")]  # sum_mod_out_table[, "Pr(>|t|)"]
          names_mod_out_coefs = names(mod_out_coefs); len_mod_out_coefs = length(mod_out_coefs)
          
          
          cnames = colnames(pred_mat); len_cnames = length(cnames); cnames_splt = strsplit(cnames, "\\_")
          save_vector_coefs = rep(0, len_cnames + 1); save_vector_pvals = rep(1, len_cnames + 1); names(save_vector_coefs) = names(save_vector_pvals) = c("(Intercept)", paste("pred_mat", cnames, sep = ""))
          
          for(j in 1:len_mod_out_coefs){
            mod_out_j = names_mod_out_coefs[j]
            
            save_vector_coefs[mod_out_j] = mod_out_coefs[j]
            save_vector_pvals[mod_out_j] = mod_out_pvals[j]
          }
          
          sum_mod[[i]]$betas = save_vector_coefs 
          sum_mod[[i]]$pvals = save_vector_pvals
          
          
          if(is.matrix(occ_mats)){
            # pvals = sum_mod[[i]]$pvals[sapply(names(sum_mod[[i]]$pvals), function(x){TRUE %in% str_detect(x, HTs_4_mod_all)})]  
            pvals = sum_mod[[i]]$pvals[unlist(lapply(strsplit(names(sum_mod[[i]]$pvals), "\\_"), function(x){x[3] %in% HTs_4_mod_all}))]
            
          } else {
            pvals = sum_mod[[i]]$pvals[which(str_detect(names(sum_mod[[i]]$pvals), "occ_mats"))]
            
          }
          
          if(TRUE %in% (pvals > pval_thresh)){
            if(is.matrix(occ_mats)){
              occ_mats_index = -which(pvals == max(pvals)); occ_mats_index = occ_mats_index[length(occ_mats_index)]
              
              cat(names(occ_mats_index), " ")
              # all.equal(str_sub(names(occ_mats_index), start = 9), colnames(occ_mats)[-occ_mats_index])
              
              occ_mats = occ_mats[, occ_mats_index]
              
              HTs_4_mod_all = HTs_4_mod_all[occ_mats_index]
              cnames_start = cnames_start[occ_mats_index]
              
              if(!is.matrix(occ_mats)){
                
                occ_mats = matrix(occ_mats, ncol = 1, dimnames = list(NULL, cnames_start[str_detect(cnames_start, HTs_4_mod_all)]))
              }
              
            } else {
              occ_mats = NULL
              
            }
            
          } else {
            repeat_model = FALSE
          }
        }
        
        if(!warning_stop){
          cat("as predictors\n")
  
          new_cols = cbind(new_cols, occ_mats)
          pred_mat_alleles = cbind(pred_mat_alleles, occ_mats)
        }
        
      } else {
        cat("There were no possible", i, "level HTs possible, so stop with searching.\n")
        
        warning_stop = TRUE
      }
      
      
    }
  }
  
  if(warning_stop){
    mod_out = previous_mod_out
    sum_mod[[i]] = sum_mod[[i - 1]]
  }

  
  out = NULL
  
  out$sum_mod = sum_mod
  
  out$new_cols = new_cols
  out$pred_mat = pred_mat_alleles
  
  if(regression == "gaussian"){
    out$sigma = sqrt(sum((predict(mod_out) - outcome)^2) / (N - nrow(sum_mod_out_table)))  # out$sigma = mod_out$sigma
  }
  
  return(out)
}
