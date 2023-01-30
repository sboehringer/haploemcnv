#' @title The Expectation-Maximization (EM-)algorithm
#' 
#' @import stats
#' @import Matrix 
#' @import stringr
#' @import glmnet
# #' @import ggplot2
#' @import philentropy
#' @import MASS
#' @import dplyr
#' @import gridExtra
#' @import gtools
#' @import rlang 
#' 
#' @description The EM-algorithm that underlies the complete profile EM-algorithm from \code{\link{HT_reconstruction}}. Haplotype frequencies are estimated in each iteration until convergence. Haplotype frequencies can be estimated 
#' 
# #' @details
#' 
#' @param lst A list. Data containing the possible diplotypes of each donor as a separate list element, can be obtained via \code{\link{all_options}}.
#' @param nl An integer. Total number of iterations till the EM-algorithm is terminated, default = 200.
#' @param freq_epsilon An numeric value. The minimal difference the biggest changing haplotype frequencies between two subsequent iterations must have to keep the EM-algorithm going. If the biggest difference is smaller, the EM-algorithm will be terminated.
#' @param likelihood_estimation A logical scalar. Whether or not the likelihood is calculated each iteration (\code{FALSE} is default). This is quite computational intensive.
#' @param likelihood_epsilon A numeric value. The minimal difference the changing likelihood between two subsequent iterations must have to keep the EM-algorithm going. Only relevant if \code{likelihood_estimation} is \code{TRUE}.
#' @param init_equal A logical scalar. Whether or not the probabilities in the initial value matrix are equal (a 0 when diplotype is not compatible with individuals genotype call, and a 1 when the diplotype is compatible) or that the number of times a compatible diplotype is stated in \code{lst} equals the value the initial value matrix should have (\code{TRUE} is default).
#' @param equal_first A logical scalar. Whether or not the estimated haplotype frequencies of the first iteration are equal (\code{FALSE} is default).
#' @param comb_vals An optional vector. The \emph{true} frequency of the combinational haplotypes is, thus at which values the \code{comb_vals} should be fixed in each iteration.
#' @param regression An optional character. Multiple options are possible: \code{NULL} (default), \code{linear}, \code{lin_ridge}, \code{logistic}. If not \code{NULL} forces the EM-algorithm to run a regression analysis each iteration which is used for haplotype frequency estimation. Then requires an outcome for each individual in \code{outcome}.  If \code{NULL}, haplotype frequency estimation will proceed without the regression analyses.
#' @param full_regression A logical scalar. Whether or not the regression models should contain all haplotype predictors next to the allelic predictors. Only functional when \code{regression} is not \code{NULL}. NOTETOSELF: Dit is nu waarschijnlijk al overbodig.
#' @param type_regression An optional character. If the EM-algorithm with outcome is requested, which inside EM-algorithm regression model needs to be run (either "alleles", "backward", "penalized"). One an also run all three by supplying "all".
#' @param outcome An optional vector. The outcome of each individual present in \code{lst}. Only functional when \code{regression} is not \code{NULL}.
#' @param status An optional vector. The status of each individual present in \code{lst}, only applicable for survival analysis. Only functional when \code{regression} is not \code{NULL}.
#' @param regression_reference A character. Which predictor, or which predictors, should be used as reference group in the regression models. 
#' @param noCNV A logical scalar. Whether or not the copy number variation haplotypes should be accomodated in the allelic predictors or should be kept as separate predictors. Only functional when \code{regression} is not \code{NULL}.
# #' @param true_betas An optional vector. The true betas of the haplotypes, if not \code{NULL} than the regression models will not be run but instead the betas here are used. Only functional when \code{regression} is not \code{NULL}. NOTETOSELF: Dit kan ook weg.
#' @param nfolds An optional integer. The group size of the cross-validation in \code{\link[glmnet]{cv.glmnet}} used for the penalized regression models. Only functional when \code{regression} is not \code{NULL}.
#' @param alpha An optional numeric value in \{0, 1\}. Which penalized regression should be performed. If \code{alpha} is 0 the ridge penalty is used, when \code{alpha} is 1 the LASSO penalty is applied and when \code{alpha} (0, 1) the elastic-net is used. Only functional when \code{regression} is not \code{NULL}.
#' @param lambda An optional integer. The number of standard errors that should be added to the optimal lambda, chosen by \code{\link[glmnet]{cv.glmnet}}. Literature suggest that lambda is 1 is the best choice, if \code{NULL}, this will be chosen. Only functional when \code{regression} is not \code{NULL}.
#' @param split_evenly An optional logical scalar. Whether or not for the calculation of the probability of the outcome given the diplotype (probYH) the probabilities are equally divided over the combinational and the non-combinatonal haplotypes, or that the combinational haplotype gets its \emph{true} frequency. Note that this \emph{true} frequency is only the true frequency for the frequency, not for probYH.
#' @param start_HTFs An optional vector. Haplotype frequencies where the diplotype frequencies of the first iteration are started with. If a haplotype is not specified here, it will get a frequency equal to half of the minimum supplied frequency
#' @param weighted_predictor A logical scalar. Whether or not each observed diplotype of an individual should get its own row where its frequency is added as a weight in the regression models (\code{TRUE} is default).
#' @param CNV_option A character. Indicating what to do with the copy number variation alleles with the inclusion of the outcome (either "first", "second", or "third"). 
#' @param probYH_epsilon A numeric value. Threshold value below which it is assumed that the outcome information is purely noise and should be discarded for a given individual.
#' 
#' @return Returns both the fully updated probability matrix with the probability vector for each donor, as the calculated haplotype frequencies for each iteration for each haplotype
#' 
# #' @examples
# #' lst = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' EM_out = EM_algorithm(lst)
#' 
#' @export EM_algorithm
EM_algorithm = function(lst, nl = 200, freq_epsilon = 1e-5, likelihood_estimation = FALSE, likelihood_epsilon = 0, init_equal = TRUE, equal_first = FALSE, start_HTFs = NULL,
                        regression = NULL, full_regression = FALSE, type_regression = "all", outcome = NULL, status = NULL, regression_reference = "NEG", comb_vals = NULL, 
                        noCNV = TRUE, nfolds = 50, alpha = 0.5, lambda = NULL, split_evenly = TRUE, weighted_predictor = TRUE, CNV_option = "second", probYH_epsilon = 1e-10){

  
  if(!is.null(start_HTFs)){
    init_equal = TRUE 
    equal_first = FALSE
  }
    
  mat = init_mat = EM_init_mat(lst = lst, init_equal = init_equal, start_HTFs = start_HTFs) 
  init_mat = Matrix(init_mat, sparse = TRUE)
  
  N = length(lst)
  nr_genes = length(unlist(strsplit(unlist(strsplit(lst[[1]][1], "\\+"))[1], "\\-")))
  
  
  diplos = colnames(mat);D = length(diplos)
  diplos_splt = sapply(diplos, function(x){strsplit(x, "\\+")})
  
  if(TRUE %in% (lengths(diplos_splt) != 2)){
    stop("A diplotype does not consists of 2 copies")
  }    
    
  haplos = sort(unique(unlist(diplos_splt))); M = length(haplos)
  haplos_splt = strsplit(haplos, "\\-")
    
  
  Repeat_mat = matrix(0, nrow = D, ncol = M, dimnames = list(diplos, haplos))
  for(j in 1:M){
    Repeat_mat[, j] = unlist(lapply(diplos_splt, function(x){sum(x == haplos[j])}))
  }
  Repeat_mat_Sp = Matrix(Repeat_mat, sparse = TRUE)
  
  
  ## HERE SOME CHECKS
  if(!is.null(regression)){
    outcome_EM = TRUE
  } else {
    outcome_EM = FALSE
  }
  
  if(outcome_EM){
    if(!(type_regression %in% c("alleles", "backward", "penalized"))){
      stop("Type of regression models within the EM-algorithm is not well specified, only have: 'alleles', 'backward' or 'penalized'")
    }

    if(is.null(full_regression)){  # Probably want to do this already at an earlier step
      if(type_regression %in% c("alleles", "backward")){
        full_regression = FALSE
      } else if(type_regression %in% c("penalized")){
        full_regression = TRUE
      }
    }
  }
  ## HERE SOME CHECKS
  
  
  only_comb_i = no_comb_i = NULL
  l = 1; e = NULL; e[l] = 1; likelihood_diff = 1; nl = nl + 1; Yi_sigmas = NULL
  
# Iterative ----
  cat("This EM-algorithm took: ")
  while(e[l] > freq_epsilon & l != nl){  #& likelihood_diff > likelihood_epsilon){)

# Haplotype frequency estimation ----
    if(l == 1){
      haplo_freqs = matrix(0, nrow = 1, ncol = M, dimnames = list(NULL, haplos))
      
      adj_index = list(NULL)
    } else {
      haplo_freqs = rbind(haplo_freqs, 0)
      
      adj_index = append(adj_index, list(NULL))
    }

    
    if(!outcome_EM){
      haplo_freqs[l, ] = colSums(Matrix(mat, sparse = TRUE) %*% Repeat_mat_Sp)
      print(haplo_freqs[l, ]);
      if(l == 1 & equal_first){
        haplo_freqs[l, ] = 1 / M
      }
      haplo_freqs[l, ] = haplo_freqs[l, ] / sum(haplo_freqs[l, ])
      
      
      
    } else {
# Predictor matrix ----
      marginal_regress.mat = full_regress.mat = NULL; marginal_genes_vect = 0  # maringal_genes_vect is redundant right?
          
      if(weighted_predictor){
        weight_pred_mat = weighted_predictor_matrix(mat = mat, type_regression = type_regression, noCNV = noCNV, CNV_option = CNV_option); pred_mat = weight_pred_mat$pred_mat
        ID = pred_mat$ID; weights = pred_mat$Weight; nr_DTs = weight_pred_mat$nr_DTs
      
        outcome_weight = rep(outcome, weight_pred_mat$nr_DTs)
        
        if(type_regression == "backward"){
          backward_selection_columns = iterative_modeling(outcome = outcome_weight, regression = regression, weighted_mat = weight_pred_mat$pred_mat, regression_reference = regression_reference, noCNV = noCNV)
          
          regress.mat = backward_selection_columns$pred_mat
          
        } else {
          regress.mat = as.matrix(pred_mat[, !(colnames(pred_mat) %in% c("ID", "Weight"))])
          
        }
        
        
        ####
        # new_mat = mat
        # for(i in 1:nrow(mat)){
        #   for(j in 1:ncol(mat)){
        #     if(new_mat[i, j] != 0){
        #       new_mat[i, j] = max(new_mat[i, j], 5e-01)
        #     }
        #   }
        # }
        # new_mat = new_mat / rowSums(new_mat)
        # 
        # weight_pred_mat = weighted_predictor_matrix(mat = new_mat, type_regression = type_regression, noCNV = noCNV, CNV_option = CNV_option); pred_mat = weight_pred_mat$pred_mat
        # regress.mat = as.matrix(pred_mat[, !(colnames(pred_mat) %in% c("ID", "Weight"))]); weights = pred_mat$Weight
        # 
        # mod_out = glm(outcome_weight ~ regress.mat, family = regression, weights = weights, control = list(maxit = 10000))
        # table(tapply(weight_pred_mat$pred_mat$Weight, weight_pred_mat$pred_mat$ID, sum))
        ####
        

      } else {
        if(nr_genes != 1){
          for(j in 1:nr_genes){
            red_mat = occurrence_coding(EM_out_mat = reducing_mat(EM_mat = mat, nr_gene = c(j)), noCNV = FALSE, gene_indication = FALSE, CNV_option = CNV_option); colnames(red_mat) = paste("G", j, "_", colnames(red_mat), sep = "")
            
            marginal_genes_vect = c(marginal_genes_vect, rep(j, ncol(red_mat)))
            marginal_regress.mat = cbind(marginal_regress.mat, red_mat)
          }
          
          if(full_regression){
            if(nr_genes > 2){
              cutted_vector_recent = cutted_vector = NULL
              for(j in 3:nr_genes){
                cutted_vector_recent = cut_vector(1:j, vect_list = cutted_vector_recent)
                
                cutted_vector = c(cutted_vector, cutted_vector_recent)
              }
              len_cutted_vector = length(cutted_vector)
              
              for(j in 1:len_cutted_vector){
                red_mat = occurrence_coding(EM_out_mat = reducing_mat(EM_mat = mat, nr_gene = cutted_vector[[j]]), noCNV = noCNV, gene_indication = FALSE); colnames(red_mat) = paste("G", paste(cutted_vector[[j]], collapse = "-"), "_", colnames(red_mat), sep = "")
                marginal_regress.mat = cbind(marginal_regress.mat, red_mat)
                
                marginal_genes_vect = c(marginal_genes_vect, rep(max(marginal_genes_vect) + 1, ncol(red_mat)))
              }
            }
            
            full_regress.mat = as.matrix(mat %*% Repeat_mat_Sp); colnames(full_regress.mat) = paste("G", paste(1:nr_genes, collapse = "-"), "_", colnames(full_regress.mat), sep = "")
            marginal_genes_vect = c(marginal_genes_vect, rep(max(marginal_genes_vect) + 1, ncol(full_regress.mat)))
            
            full_regress.mat = full_regress.mat[, colSums(full_regress.mat) > (1 / (2 * N))]
            marginal_regress.mat = marginal_regress.mat[, colSums(marginal_regress.mat) > (1 / (2 * N))]
            
          } else if(type_regression == "backward"){
            backward_selection_columns = iterative_modeling(EM_mat = mat, outcome = outcome, regression = regression, regression_reference = regression_reference, noCNV = noCNV)$new_cols
            marginal_regress.mat = cbind(marginal_regress.mat, backward_selection_columns)
            
          }
          
          regress.mat = cbind(marginal_regress.mat, full_regress.mat)

          # comb_index = which(str_detect(colnames(regress.mat), "Comb"))
          # # if(length(comb_index) != 0){
          #   regress.mat = regress.mat[, -comb_index]
          # # } else {
          # #   cat("In iteration", l, "there was no Comb allele or HT in the predictor matrix...")
          # # }
          # # regress.mat = regress.mat[, -comb_index[duplicated(unlist(lapply(strsplit(colnames(regress.mat)[comb_index], "\\_"), function(x){x[2]})))]]
          
          
        } else {
          regress.mat = as.matrix(mat %*% Repeat_mat_Sp); colnames(regress.mat) = paste("G1_", colnames(regress.mat), sep = "")
          
          marginal_genes_vect = c(marginal_genes_vect, rep(1, ncol(regress.mat)))
        }
      }
      

          
      names_regress.mat = colnames(regress.mat)
      if(l == 1){
        regress.mat_index = which(str_detect(names_regress.mat, regression_reference) | str_detect(names_regress.mat, "Comb"))  # Added COMB as extra regression_reference
        regress.mat_index_haplo = unique(names_regress.mat[regress.mat_index])
        
        regress.mat = regress.mat[, !(names_regress.mat %in% regress.mat_index_haplo)]
        marginal_genes_vect = marginal_genes_vect[c(TRUE, !(names_regress.mat %in% regress.mat_index_haplo))]
          
        names_regress.mat = colnames(regress.mat)
          
          
        haplos_sub = haplos[!(haplos %in% regress.mat_index_haplo)]
        if(noCNV){
          haplos_sub_splt = lapply(strsplit(haplos_sub, "\\-"), strsplit, "\\^")
        } else {
          haplos_sub_splt = lapply(strsplit(haplos_sub, "\\-"), as.list)
        }
        # max_haplos_sub_splt = unlist(lapply(haplos_sub_splt, function(x){max(lengths(x))}))
        max_haplos_sub_splt = nr_genes / unlist(lapply(haplos_sub_splt, function(x){sum(lengths(x))}))
        
        betas = matrix(0, nrow = 1, ncol = (length(haplos_sub) + 1), dimnames = list(NULL, c("Mu", haplos_sub)))
        betas_mod = matrix(0, nrow = 1, ncol = (length(names_regress.mat) + 1), dimnames = list(NULL, c("Mu", names_regress.mat)))
        
        
        lst_haplo = rep(list(NULL), N)
        for(i in 1:N){
          lst_haplo[[i]] = unique(unlist(strsplit(lst[[i]], "\\+")))
        }
          
        
        if(regression == "survival"){
          pred.mat = as.data.frame(Repeat_mat); colnames(pred.mat) = sapply(colnames(pred.mat), function(x){paste(unlist(strsplit(x, "\\-")), collapse = ".")})
        }
        
        foldID = foldID_setup = NULL
        if(type_regression == "penalized"){
          ###
          foldID_setup = rep(nfolds, N)
          
          Ntosample = 1:N
          
          Nsize = N / nfolds
          Nsizeceiling = ceiling(Nsize); Nsizefloor = floor(Nsize)
          
          Nsize = rep(Nsizefloor, (nfolds - 1))
          if(Nsizeceiling != Nsizefloor){
            modulo = N %% nfolds
            Nsize[1:modulo] = Nsize[1:modulo] + 1
          }
          
          for(m in 1:(nfolds - 1)){
            samp = sort(sample(x = Ntosample, size = Nsize[m], replace = FALSE))
            
            foldID_setup[samp] = m
            Ntosample = Ntosample[!(Ntosample %in% samp)]
          }
          
          if(weighted_predictor){
            foldID = rep(foldID_setup, nr_DTs)
          } else {
            foldID = foldID_setup
          }
          ###
          
          lambdas = NULL
        }  
        
          
          
      } else {
        regress.mat = regress.mat[, !(names_regress.mat %in% regress.mat_index_haplo)]
        names_regress.mat = colnames(regress.mat)
        
        betas = rbind(betas, 0)
        betas_mod = rbind(betas_mod, 0)
        

        if(FALSE %in% (names_regress.mat %in% colnames(betas_mod))){
          betas_mod_index = which(!(names_regress.mat %in% colnames(betas_mod)))
          betas_mod_extra = matrix(0, nrow = l, ncol = length(betas_mod_index), dimnames = list(NULL, names_regress.mat[betas_mod_index]))
          
          betas_mod = cbind(betas_mod, betas_mod_extra)
        }  
        
        if(weighted_predictor){
          foldID = rep(foldID_setup, nr_DTs)
        } else {
          foldID = foldID_setup
        }
      }
        

      regress.mat = round(regress.mat, 50)  
      
      

# Inside regression analyses ----
      if(type_regression %in% c("alleles", "backward")){
        if(!weighted_predictor){
          mod_out = glm(outcome ~ regress.mat, family = regression, control = list(maxit = 10000))
        } else {
          mod_out = glm(outcome_weight ~ regress.mat, family = regression, weights = weights, control = list(maxit = 10000))
        }
          
        mod_out_coef = mod_out$coefficients; mod_out_coef_names = str_sub(names(mod_out_coef), start = 12); mod_out_coef_names[1] = "G0_Mu"
          
      } else if(type_regression == "penalized"){
        if(!weighted_predictor){
          mod_out = cv.glmnet(y = outcome, x = regress.mat, family = regression, alpha = alpha, nlambda = 200, foldid = foldID)
        } else {
          mod_out = cv.glmnet(y = outcome_weight, x = regress.mat, weights = weights, family = regression, alpha = alpha, nlambda = 200, foldid = foldID)
        }
  
        if(is.null(lambda) || lambda == 1){
          mod_out_lambda = mod_out$lambda.1se; mod_out_index = mod_out$index["1se", "Lambda"]
            
        } else {
          mod_out_lambda = mod_out$lambda.min + ((mod_out$lambda.1se - mod_out$lambda.min) * lambda); mod_out_lambda_diff = abs(mod_out_lambda - mod_out$lambda)
          mod_out_index = which(mod_out_lambda_diff == min(mod_out_lambda_diff))
            
        }
        
        # cat("Iteration ", l, " has lambda ", mod_out_lambda, "\n")
          
        mod_out_coef = coef(mod_out, s = mod_out_lambda)[, 1]; mod_out_coef_names = names(mod_out_coef); mod_out_coef_names[1] = "G0_Mu"
      }
      
       
      mod_out_coef_names_splt = strsplit(mod_out_coef_names, "\\_")
      mod_out_coef_names_nr_gene = lapply(mod_out_coef_names_splt, function(x){as.numeric(unlist(strsplit(str_sub(x[1], start = 2), "\\-")))})
      mod_out_coef_names_HTs = unlist(lapply(mod_out_coef_names_splt, function(x){x[2]}))
          
      splt_mods_out_coef = strsplit(mod_out_coef_names_HTs, "\\-"); lens_splt_mods_out_coef = lengths(splt_mods_out_coef)
      for(j in 1:length(mod_out_coef)){
        mod_out_coef_j = mod_out_coef[j]; mod_out_coef_names_j = mod_out_coef_names[j]
        mod_out_coef_names_nr_gene_j = mod_out_coef_names_nr_gene[[j]]; mod_out_coef_names_HTs_j = mod_out_coef_names_HTs[j]
            
        if(j == 1){
          betas[l, "Mu"] = mod_out_coef[j]
          betas_mod[l, "Mu"] = mod_out_coef[j]
          
          next
              
        # } else if(lens_splt_mods_out_coef[j] == 1){
        #   beta_index_val = c(FALSE, unlist(lapply(haplos_sub_splt, function(x){sum(unlist(x[mod_out_coef_names_nr_gene_j]) == mod_out_coef_names_HTs_j)})))
        #   beta_index = ifelse(beta_index_val != 0, TRUE, FALSE)
        #   
        #   # if(noCNV){
        #   #   beta_index_val = c(FALSE, unlist(lapply(haplos_sub_splt, function(x){sum(unlist(x[mod_out_coef_names_nr_gene_j]) == mod_out_coef_names_HTs_j)})))
        #   #   beta_index = ifelse(beta_index_val != 0, TRUE, FALSE)
        #   #   
        #   # } else {
        #   #   beta_index = c(FALSE, unlist(lapply(haplos_sub_splt, function(x){x[mod_out_coef_names_nr_gene_j] == mod_out_coef_names_HTs_j})))
        #   #   beta_index_val = as.numeric(beta_index)
        #   # }
        #   
        # ##
        } else {
          options = lapply(haplos_sub_splt, function(x){apply(expand.grid(x[mod_out_coef_names_nr_gene_j]), 1, paste, collapse = "-")})
          if(CNV_option == "first"){
            lens_options = 1; weighting = 1
            
          } else if(CNV_option == "second"){
            lens_options = lengths(options); weighting = 1
            
          } else if(CNV_option == "third"){
            lens_options = lengths(options); weighting = unlist(lapply(haplos_sub_splt, function(x){length(unlist(x[mod_out_coef_names_nr_gene_j]))})) / length(mod_out_coef_names_nr_gene_j)
              
          }
            
          ## When noCNV == TRUE, each HT can be different, with differing lengths of options and weights
          ## When noCNV == FALSE, each HT remains its HT. So the lengths of options and weighting are always 1
          
          beta_index_val = c(FALSE, unlist(lapply(options, function(x){sum(x %in% mod_out_coef_names_HTs_j)})) / lens_options * weighting)
          beta_index = ifelse(beta_index_val != 0, TRUE, FALSE)
        }
          
        betas[l, beta_index] = betas[l, beta_index] + mod_out_coef_j * beta_index_val[beta_index]
        betas_mod[l, mod_out_coef_names_j] = mod_out_coef_j
      }
        
        
      if(TRUE %in% is.na(betas[l, ])){
        temp_betas = betas[l, ]
        temp_betas[is.na(temp_betas)] = 0
          
        betas_diplo = betas[l, "Mu"] + colSums(apply(Repeat_mat, 1, function(x){x * temp_betas[-1]}))  # Repeat_mat[, -regress.mat_index]
      } else {
        betas_diplo = betas[l, "Mu"] + colSums(apply(Repeat_mat, 1, function(x){x * betas[l, -1]}))  # Repeat_mat[, -regress.mat_index]
      }
      
      
        
      if(regression == "gaussian"){
        if(type_regression %in% c("alleles", "backward")){
          if(!weighted_predictor){
            Yi_sigma = sqrt(sum((predict(mod_out) - outcome)^2) / (N - sum(!is.na(mod_out_coef))))  # mod_out$df.residual  # summary(mod_out)$sigma
          } else {
            Yi_sigma = sqrt(sum(weights * (predict(mod_out) - outcome_weight)^2) / (length(weights) - sum(!is.na(mod_out_coef))))
            # Yi_sigma = sqrt(mod_out$deviance / mod_out$df.residual)  # sigma(mod_out)
          }
          
        } else if(type_regression == "penalized"){
          if(!weighted_predictor){
            Yi_sigma = sqrt(sum((outcome - predict(mod_out, regress.mat, s = mod_out_lambda))^2) / (N - mod_out$glmnet.fit$df[mod_out_index] + 1))
          } else {
            Yi_sigma = sqrt(sum((outcome_weight - predict(mod_out, regress.mat, s = mod_out_lambda))^2) / (N - mod_out$glmnet.fit$df[mod_out_index] + 1))
          }
          
        }
        Yi_sigmas[l] = Yi_sigma
          
      } else if(regression == "binomial"){
        prob_diplo = 1 / (1 + exp(betas[l, "Mu"] + colSums(apply(Repeat_mat, 1, function(x){x * betas[l, -1]}), na.rm = TRUE)))
          
      }
        
      
      
# Calculation probYH ----
      probYH = matrix(0, nrow = N, ncol = D, dimnames = list(NULL, diplos))
      for(i in 1:N){
        comp_diplo = lst[[i]]; index_k = which(diplos %in% comp_diplo)
        comp_haplo = lst_haplo[[i]]; index_j = which(haplos %in% comp_haplo)
         
        if(regression == "gaussian"){
          probYH[i, index_k] = dnorm(outcome[i], betas_diplo[index_k], Yi_sigma)
        
        } else if(regression == "binomial"){
          probYH[i, index_k] = dbinom(outcome[i], 1, prob_diplo[index_k])
          
        }
        
        if(sum(probYH[i, index_k]) == 0){
          # cat("For donor", i, "all probYHs were 0 and are thus given the same value\n")
          probYH[i, index_k] = 1
          
        # } else if(sum(probYH[i, index_k]) < probYH_epsilon){
        } else if(sum(probYH[i, index_k[mat[i, index_k] != 0]]) < probYH_epsilon){
#          cat("For donor", i, "the sum of the probYHs was lower than", probYH_epsilon, "and are thus given the same value\n")
          probYH[i, index_k] = 1
          
          adj_index[[l]] = c(adj_index[[l]], i)
        }
          
        # } else if(!(TRUE %in% (which(mat[i, index_k] != 0) %in% which(probYH[i, index_k] != 0)))){
        #   cat("For donor", i, "probYHs and mat had no overlap in non-zero DTs\n")
        #   
        #   probYH_index = which(probYH[i, index_k] == 0); probYH_min = min(probYH[i, index_k[-probYH_index]])
        #   probYH[i, index_k[probYH_index]] = probYH_min * 0.01
        # }
        
        
        # if(i == 1379){
        #   cat("the unnormalized probYH:\n")
        #   print(probYH[i, index_k])
        #   cat("\n")
        # }
        
        
        probYH[i, index_k] = probYH[i, index_k] / sum(probYH[i, index_k])
      }
      
      # if(!is.null(adj_index[[l]])){
      #   cat("For donors", adj_index[[l]], "the sum of the probYHs was lower than 0.001 and are thus given the same value\n")
      # }
    }
    


    
    
# Fixing probYH comb_values ----
    if(!is.null(comb_vals) & outcome_EM){
      comb_vals_index = comb_vals[length(comb_vals)]

      for(i in 1:N){
        comp_diplo = lst[[i]]; index_k = which(diplos %in% comp_diplo)
        comp_haplo = lst_haplo[[i]]; index_j = which(haplos %in% comp_haplo)

        comb_index_probYH = str_detect(names(probYH[i, index_k]), names(comb_vals_index))

        if(TRUE %in% comb_index_probYH){
          
          if(sum(probYH[i, index_k[!comb_index_probYH]]) != 0 & sum(probYH[i, index_k[comb_index_probYH]]) != 0){

            if(split_evenly){
              probYH[i, index_k[!comb_index_probYH]] = probYH[i, index_k[!comb_index_probYH]] / sum(probYH[i, index_k[!comb_index_probYH]]) * 0.5
              probYH[i, index_k[comb_index_probYH]] = probYH[i, index_k[comb_index_probYH]] / sum(probYH[i, index_k[comb_index_probYH]]) * 0.5

            } else {
              probYH[i, index_k[!comb_index_probYH]] = probYH[i, index_k[!comb_index_probYH]] / sum(probYH[i, index_k[!comb_index_probYH]]) * (1 - comb_vals_index)
              probYH[i, index_k[comb_index_probYH]] = probYH[i, index_k[comb_index_probYH]] / sum(probYH[i, index_k[comb_index_probYH]]) * comb_vals_index

            }

          } else if(sum(probYH[i, index_k[!comb_index_probYH]]) == 0){
            only_comb_i = unique(c(only_comb_i, i))

          } else if(sum(probYH[i, index_k[comb_index_probYH]]) == 0){
            no_comb_i = unique(c(no_comb_i, i))

          }
        }
        
        # if(i == 1379){
        #   cat("The standardized probYH:\n")
        #   print(probYH[i, index_k])
        #   cat("\n")
        # }
        
      }
    }
    
    
    
# Calculating x_{ij} ----
    if(outcome_EM){
      xij = matrix(0, nrow = N, ncol = M, dimnames = list(NULL, haplos))
      for(i in 1:N){
        comp_diplo = lst[[i]]; index_k = which(diplos %in% comp_diplo)
        comp_haplo = lst_haplo[[i]]; index_j = which(haplos %in% comp_haplo)
        
        xij.mat_i = probYH[i, index_k] * mat[i, index_k]
        # if(i == 1379){
        #   cat("xij.mat_i:\n")
        #   print(xij.mat_i)
        #   cat("\n")
        #   
        #   cat("mat:\n")
        #   print(mat[i, index_k])
        #   cat("\n")
        # }
        xij.mat_i = xij.mat_i / sum(xij.mat_i)
        
        
        for(j in index_j){
          xij[i, j] = sum(Repeat_mat[index_k, haplos[j]] * xij.mat_i)
        }
      }
      
      haplo_freqs[l, ] = colSums(xij) / (2 * N)
    }
    
    

# Fixing comb_values ----
    if(!is.null(comb_vals)){
      if(l == 1){
        comb_names = names(comb_vals)
        haplos_splt = strsplit(haplos, "\\-")
        
        index_haplos_w_comb = unlist(lapply(haplos_splt, function(x){TRUE %in% (comb_names %in% x)}))
        index_comb_w_haplos = apply(sapply(comb_names, function(x){lapply(haplos_splt, function(y){x %in% y})}), 2, function(z){TRUE %in% z})
        len_comb_index = sum(index_comb_w_haplos)
      }
      
      haplo_freqs[l, ] = changing_AFs_in_HFs(haplotypes = haplo_freqs[l, ], comb_vals = comb_vals)

      if(len_comb_index == 0){
        haplo_freqs[l, ] = haplo_freqs[l, ] / sum(haplo_freqs[l, ])
      
      } else if(len_comb_index == 1){
        haplo_freqs[l, index_haplos_w_comb] = haplo_freqs[l, index_haplos_w_comb] / sum(haplo_freqs[l, index_haplos_w_comb]) * comb_vals[index_comb_w_haplos]
        haplo_freqs[l, !index_haplos_w_comb] = haplo_freqs[l, !index_haplos_w_comb] / sum(haplo_freqs[l, !index_haplos_w_comb]) * (1 - comb_vals[index_comb_w_haplos])
        
      } else {
        comb_indeces = as.list(rep(0, len_comb_index))
        for(j in 1:len_comb_index){
          comb_indeces[[j]] = which(index_haplos_w_comb)[j]
            
          haplo_freqs[l, comb_indeces[[j]]] = haplo_freqs[l, comb_indeces[[j]]] / sum(haplo_freqs[l, comb_indeces[[j]]]) * comb_vals[which(index_comb_w_haplos)[j]]
        }
          
        comb_index = unique(unlist(comb_indeces))
        haplo_freqs[l, -comb_index] = haplo_freqs[l, -comb_index] / sum(haplo_freqs[l, -comb_index]) * (1 - sum(comb_vals[index_comb_w_haplos]))
      }
      
      if(is.finite(FALSE %in% haplo_freqs[l, ])){
        haplo_freqs[l, which(!is.finite(haplo_freqs[l, ]))] = 0
      }
    }
    

    
# Diplotype frequency estimation ----
    diplo_freqs = vector(mode = "numeric", length = D); names(diplo_freqs) = diplos
    for(d in 1:D){
      cname = diplos_splt[[d]]
          
      if(cname[1] == cname[2]){
        diplo_freqs[d] = haplo_freqs[l, cname[1]] * haplo_freqs[l, cname[2]]
      } else {
        diplo_freqs[d] = haplo_freqs[l, cname[1]] * haplo_freqs[l, cname[2]] * 2
      }
    }
      
    index_no0 = which(diplo_freqs != 0)
    diplo_freqs[index_no0] = diplo_freqs[index_no0] / sum(diplo_freqs)


  
# Updating matrix ----
    if(outcome_EM){
      mat = t(apply(init_mat, 1, function(x){x * diplo_freqs})) * probYH
    } else {
      mat = t(apply(init_mat, 1, function(x){x * diplo_freqs}))
    }
    mat = mat / rowSums(mat) 
    
    
    
    
# Likelihood ----
    if(l == 1){
      obs_likelihoods = NULL
    }
#    obs_likelihoods[l] = sum(apply(init_mat, 1, function(x){log(sum(diplo_freqs * x))}))
#    likelihood_diff = abs(obs_likelihoods[i - 1] - obs_likelihoods[l])
      
    if(likelihood_estimation){
      if(l == 1){
        est_likelihoods = NULL
        
        likeli_iter = likelihood_haplotypes(mat, haplos_i = haplo_freqs[l, ])
        indicator = likeli_iter$Indicator
        
        est_likelihoods[l] = likeli_iter$Likelihood
      } else {
        est_likelihoods[l] = likelihood_haplotypes(mat, haplos_i = haplo_freqs[l, ], indicator = indicator)$Likelihood
      }
    }
    
    if(l == 1){
      e[l + 1] = 1
    } else {
      e[l + 1] = max(abs(haplo_freqs[l, ] - haplo_freqs[l - 1, ]))
    }
    
    cat(l, " ")
    l = l + 1
  }
  cat("iterations \n")
  
  
  e = e[-1]
  
  if(l == nl){
    cat("The number of iterations equals the maximum allowed iterations, epsilon has finished at", e[l - 1], "\n")
  }
  
  
  
  
  
  out = NULL
  
  out$Matrix = as.matrix(mat)
  out$Frequencies = haplo_freqs
  out$Epsilon = e
  out$Obs_likelihood = obs_likelihoods
  out$Yi_sigmas = Yi_sigmas
  
  if(likelihood_estimation == TRUE){
    out$Est_likelihood = est_likelihoods
  }
  if(outcome_EM){
    out$adj_index = adj_index
    out$betas = betas; out$betas_mod = betas_mod
    
    if(type_regression == "penalized"){
      out$lambdas = lambdas
    }
    
    out$no_comb_i = sort(no_comb_i)
    out$only_comb_i = sort(only_comb_i)
  }
  
  return(out)
}
