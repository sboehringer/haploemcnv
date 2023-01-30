#'
#'
#'
actual_analysis_outside_regression_summary = function(all_betas, true_betas){
 
  N = nrow(all_betas)
  
  cmeans = sqrt(colMeans(t(apply(all_betas, 1, function(x){(x - true_betas[names(x)])^2})), na.rm = TRUE)) 
  cmeans_median = abs(apply(t(apply(all_betas, 1, function(x){x - true_betas[names(x)]})), 2, function(y){median(y, na.rm = TRUE)}))
  
  beta_empSE = colSums(t(apply(all_betas, 1, function(x){(x - true_betas[names(x)])^2})), na.rm = TRUE) / (N - 1)
  
  MCSE_beta = sqrt(colSums(t(((apply(all_betas, 1, function(x){(x - true_betas[names(x)])^2}) - (cmeans^2))^2))) * (1 / (N * (N - 1))))
  
   
  return(list("cmeans" = cmeans, "cmeans_median" = cmeans_median, "beta_empSE" = beta_empSE, "MCSE_beta" = MCSE_beta))
}
  