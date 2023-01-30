#'
#'
#'
actual_analysis_HT_summary = function(all_HTFs = NULL, last_HTFs, all_true_HTFs, HTF_selection_sim, HTF_selection_max_val, selection_index = 1, 
                                      nr_config = 1, paper1 = FALSE, selection_threshold = 0.01){
  
  N = nrow(last_HTFs)
  nr_genes = length(unlist(strsplit(colnames(last_HTFs)[1], "\\-")))
  
  if(paper1){
    all_HTFs_index = all_HTFs[, lengths(strsplit(colnames(all_HTFs), "\\-")) != nr_genes]
    last_HTFs = cbind(all_HTFs_index, last_HTFs)
      
    HTF_selection_sim[[selection_index]] = unique(c(HTF_selection_sim[[selection_index]], colnames(all_HTFs_index)))
  } else {
    HTF_selection_sim[[selection_index]] = unique(c(HTF_selection_sim[[selection_index]], colnames(last_HTFs)[apply(last_HTFs, 2, function(x){TRUE %in% (x > selection_threshold)})]))
      
      
    temp_cnames = colnames(last_HTFs); len_temp_cnames = length(temp_cnames)
    if(FALSE %in% (temp_cnames %in% names(HTF_selection_max_val[[selection_index - 1]]))){
      max_val_index = !(temp_cnames %in% names(HTF_selection_max_val[[selection_index - 1]]))
        
      new_HTF_selection_max_val = rep(0, sum(max_val_index)); names(new_HTF_selection_max_val) = temp_cnames[max_val_index]
      HTF_selection_max_val[[selection_index - 1]] = c(HTF_selection_max_val[[selection_index - 1]], new_HTF_selection_max_val)
    }
    
    for(y in 1:len_temp_cnames){
      temp_cname = temp_cnames[y]
        
      HTF_selection_max_val[[selection_index - 1]][temp_cname] = max(HTF_selection_max_val[[selection_index - 1]][[temp_cname]], last_HTFs[!is.na(last_HTFs[, temp_cname]), temp_cname])
    }
  } 
    
    
  cnames_last_HTFs = colnames(last_HTFs); len_cnames_last_HTFs = length(cnames_last_HTFs)
  true_HTFs = rep(0, len_cnames_last_HTFs); names(true_HTFs) = cnames_last_HTFs
    
  name_index = cnames_last_HTFs[cnames_last_HTFs %in% names(all_true_HTFs)]
  true_HTFs[name_index] = all_true_HTFs[name_index]

    
  SE = t(apply(last_HTFs, 1, function(x){(x - true_HTFs[names(x)])^2}))
  SE = ifelse(!is.na(SE), SE, 0)
    
  bias1 = sqrt(colMeans(SE))
  bias2 = bias1 / true_HTFs[names(bias1)]
    
  MCSE_HT = sqrt(colSums(t(((apply(last_HTFs, 1, function(x){(x - true_HTFs[names(x)])^2}) - (bias1^2))^2))) * (1 / (N * (N - 1))))
    
  bias_median = abs(apply(t(apply(last_HTFs, 1, function(x){x - true_HTFs[names(x)]})), 2, function(x){median(x, na.rm = TRUE)}))
  
  return(list("bias1" = bias1, "bias2" = bias2, "bias_median" = bias_median, "MCSE" = MCSE_HT,
              "HTF_selection_sim" = HTF_selection_sim, "HTF_selection_max_val" = HTF_selection_max_val))
}
