#' @title P-value calculation for penalized regression
#' 
#' @description Due to the bias introduced in the estimates due to the penalty, p-values are not (entirely) valid. However, they can be estimated, as described in \cite{Cule et al., 2011}. Support function for \code{\link{iterative_modelling}}.
#' 
#' @param pred_mat A matrix. Predictor matrix used in the penalized regression model.
#' @param betas A numeric value. Betas estimated in the penalized regression model with the chosen \code{lambda}.
#' @param lambda A numeric value. Lambda chosen for the penalized regression model.
#' 
#' @return A vector with estimated p-values for the predictors.
#' 
penalized_significance = function(pred_mat, betas, lambda){
  m = length(betas) - 1
  
  tXX = t(pred_mat) %*% pred_mat
  tXXlIm1 = (tXX + diag(lambda, nrow = m))^-1
  
  H = pred_mat %*% tXXlIm1 %*% t(pred_mat)
  df = n - sum(diag((2 * H) - (H %*% t(H))))
#  df = n - m
  
  regression_sigma = (t(yi - (pred_mat %*% betas[-1])) %*% (yi - (pred_mat %*% betas[-1]))) / df
  regression_SE = sqrt(diag(as.numeric(regression_sigma) * tXXlIm1 %*% tXX %*% tXXlIm1))
  
  T_lam = betas[-1] / regression_SE; T_lam = ifelse(is.na(T_lam), 0, T_lam)
  HT_pvals = 2 * pt(-abs(T_lam), df = df)
  names(HT_pvals) = c("(Intercept)", sapply(names(HT_pvals)[!str_detect(names(HT_pvals), "Intercept")], function(x){paste("pred_mat", x, sep = "")}))
  
  
  return(HT_pvals)
}
