#' @title Profile EM-algorithm for a single starting gene
#' 
#' @description The best order in which the genes are added into the haplotype reconstruction is based via some heuristic strategy. The gene in \code{gene1} is the starting gene, and all other genes that need to be added into the reconstruction are denoted in \code{other_genes}. The profile EM-algorithm is described in \cite{}.
#' 
#' @param gene1 A list. Diplotype list of the starting gene, can be obtained via \code{\link{all_options}}.
#' @param other_genes An optional list. Diplotype list of the other genes that need to be added into the reconstruction, can be obtained via \code{\link{all_options}}. 
#' @param ... The diplotype lists of the other genes can also be supplied without specification of \code{other_genes}.
#' @param gene_names An optional vector. Names of genes, need to be in the order in which they are supplied.
#' @param orderings An optional vector. Numerical order of the genes.
#' @param comb_vals An optional vector. The \emph{true} frequency of the combinational haplotypes is, thus at which values the \code{comb_vals} should be fixed in each iteration.
#' @param naive A logical scalar. Whether or not the naive analyses also should be run with the data
#' @param fixed_order A logical scalar. Whether or not the order in which the genes are grouped is based on the order in which they are supplied, or that the order is based on some heuristic strategy.
#' @param candidate_list An optional vector. Which haplotypes are never allowed to be grouped into the compound haplotype.
#' @param comb_gene1 A logical scalar. Whether or not the alleles of the starting gene need to be grouped before the reconstruction starts (\code{TRUE} is default).
#' @param combine_haplos A logical scalar. Whether or not low-frequency haplotypes in the haplotype reconstruction are grouped (\code{TRUE} is default). 
#' @param CO_thresh An optional numeric value in \{0, 1\}. The cut-off threshold, only haplotypes lower than this threshold are be grouped.
#' @param CO_perc An optional numerical value in \{0, 100\}. The cut-off percentage, only the haplotypes belonging to the percentage of haplotypes with the lowest frequency are eligible for grouping.
#' @param CO_min An optional integer. The cut-off minimum, the minimum number of haplotypes which need to be grouped to actually group haplotypes.
#' @param all_previous A logical scalar. Whether or not all previous compound haplotypes should be grouped in the newly formed compound haplotype. The name of all previous compound haplotypes should start with \emph{Comb}.
#' @param outcome_EM A logical scalar. Whether or not the EM-algorithm with outcome must be run (\code{TRUE} is default).
#' @param regression An optional character. Multiple options are possible: \code{NULL} (default), \code{linear}, \code{lin_ridge}, \code{logistic}. If not \code{NULL} forces the EM-algorithm to run a regression analysis each iteration which is used for haplotype frequency estimation. Then requires an outcome for each individual in \code{outcome}.  If \code{NULL}, haplotype frequency estimation will proceed without the regression analyses.
#' @param type_regression An optional character. If the EM-algorithm with outcome is requested, which inside EM-algorithm regression model needs to be run (either "alleles", "backward", "penalized"). One an also run all three by supplying "all".
#' @param outcome An optional vector. The outcome of each individual present in \code{lst}. Only functional when \code{regression} is not \code{NULL}.
#' @param status An optional vector. The status of each individual present in \code{lst}, only applicable for survival analysis. Only functional when \code{regression} is not \code{NULL}.
#' @param full_regression A logical scalar. Whether or not the regression models should contain all haplotype predictors next to the allelic predictors. Only functional when \code{regression} is not \code{NULL}.
#' @param regression_reference A character. Which predictor, or which predictors, should be used as reference group in the regression models. 
#' @param noCNV a logical scalar. Whether or not the copy number variation haplotypes should be accomodated in the allelic predictors or should be kept as separate predictors. Only functional when \code{regression} is not \code{NULL}.
#' @param nfolds An optional integer. The group size of the cross-validation in \code{\link[glmnet]{cv.glmnet}} used for the penalized regression models. Only functional when \code{regression} is not \code{NULL}.
#' @param alpha An optional numeric value in \{0, 1\}. Which penalized regression should be performed. If \code{alpha} is 0 the ridge penalty is used, when \code{alpha} is 1 the LASSO penalty is applied and when \code{alpha} (0, 1) the elastic-net is used. Only functional when \code{regression} is not \code{NULL}.
#' @param lambda An optional integer. The number of standard errors that should be added to the optimal lambda, chosen by \code{\link[glmnet]{cv.glmnet}}. Literature suggest that lambda = 1 is the best choice, if \code{NULL}, this will be chosen. Only functional when \code{regression} is not \code{NULL}.
#' @param CNV_option A character. Indicating what to do with the copy number variation alleles with the inclusion of the outcome (either "first", "second", or "third"). 
#' @param weighted_predictor A logical scalar. Whether or not each observed diplotype of an individual should get its own row where its frequency is added as a weight in the regression models (\code{TRUE} is default).
#' 
#' @return While running the different combinations and iterations, the order in which the genes are combined are concatenated. After each inclusion, the included gene is removed from the list and cannot be included again. The output is all grouping information of the last grouping and the retained haplotypes, additionally, EM-algorithm output for the intermediate and last reconstruction is supplied.
#' 
# #' @examples 
# #' gene1 = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' gene2 = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' gene3 = list(c("007+008", "007+009"), "NEG+NEG", c("008+008", "008+NEG"))
# #' 
# #' \dontrun{
# #' grouping_order_all(gene1, gene2, gene3)
# #' grouping_order_all(gene1, other_genes = list(gene2, gene3))
# #' }
#' 
grouping_order_all = function(gene1, ..., other_genes = NULL, gene_names = NULL, orderings = NULL, comb_vals = NULL, naive = TRUE, fixed_order = FALSE, candidate_list = NULL, 
                              comb_gene1 = TRUE, combine_haplos = TRUE, CO_thresh = 1e-7, CO_perc = 100, CO_min = 1, all_previous = TRUE, outcome_EM = FALSE,
                              regression = "gaussian", full_regression = FALSE, type_regression = "all", regression_reference = "NEG", outcome = NULL, status = NULL, noCNV = TRUE, CNV_option = "second",
                              nfolds = 50, alpha = 0.5, lambda = NULL, weighted_predictor = TRUE){
  
  
  if(is.null(other_genes)){
    other_genes = list(...)
  }
  nr_other_genes = length(other_genes); nr_genes = nr_other_genes + 1
  
  if(is.null(gene_names)){
    gene_names = 1:(nr_other_genes + 1)
  }
  if(is.null(orderings)){
    orderings = 1:length(gene_names)
  }
  
  out_gene1 = EM_algorithm(lst = gene1, comb_vals = comb_vals, regression = NULL)
  all_EMs_before = list(out_gene1)
  
  if(!is.null(candidate_list)){
    candidate_list_splt = strsplit(candidate_list, "\\_"); CL_genes = str_sub(unlist(lapply(candidate_list_splt, function(x){x[1]})), start = 2); CL_HTs = unlist(lapply(candidate_list_splt, function(x){x[2]}))
    first_gene = combining_AF(lst = gene1, EM_out = out_gene1, CO_thresh = CO_thresh, CO_perc = CO_perc, CO_min = CO_min, comb_name = "CombG1", excluding_haplo = CL_HTs[CL_genes == 1])  # what if there is a different starting gene...
  } else {
    first_gene = combining_AF(lst = gene1, EM_out = out_gene1, CO_thresh = CO_thresh, CO_perc = CO_perc, CO_min = CO_min, comb_name = "CombG1", excluding_haplo = candidate_list)
  }

  if(length(first_gene$comb) != 0){
    if("CombG1" %in% names(comb_vals)){
      comb_vals["CombG1"] = comb_frequency(haplotypes = first_gene$haplotypes, combine = first_gene$comb)
    } else {
      comb_vals = c(comb_vals, "CombG1" = comb_frequency(haplotypes = first_gene$haplotypes, combine = first_gene$comb))
    }
  }
  
  if(comb_gene1 == TRUE & length(first_gene$comb) != 0){
    gene1_comb = first_gene$comb  # [[1]]
    gene1 = first_gene$lst
    
  } else if(length(first_gene$comb) == 0){
    comb_gene1 = FALSE
    
  }
  
  out_gene1 = EM_algorithm(lst = gene1, comb_vals = comb_vals, regression = NULL)

  reference_freqs = out_gene1$Frequencies[nrow(out_gene1$Frequencies), ]
  if(comb_gene1 == TRUE){
    reference_info = list("freqs" = reference_freqs,
                           "comb" = gene1_comb)  
  } else {
    reference_info = list("freqs" = reference_freqs,
                           "comb" = NULL)
  }
  
  if(comb_gene1 == TRUE){
    final_freqs = reference_freqs
  } else {
    final_freqs = reference_freqs[first_gene$remain]
  }
  
  
  if(naive){
    all_naive_out = rep(list(NULL), (length(other_genes) + 1))
    all_naive_out[[1]] = Naive_EMbased_analysis(gene1, haplo = TRUE)  #, adj_comb_vals = comb_vals)
  }
  

  all_lists = list(gene1)
  all_EMs_after = list(out_gene1)
  
  GOs = list(NULL)
  for(i in 1:nr_other_genes){
    GO_i = grouping_order_iter(gene1 = gene1, other_genes, reference_info = reference_info, combine_haplos = combine_haplos, regression = NULL,
                               all_previous = all_previous, orderings = orderings, fixed_order = fixed_order, gene_names = gene_names, candidate_list = candidate_list, 
                               comb_vals = comb_vals, CO_perc = CO_perc, CO_thresh = CO_thresh, CO_min = CO_min)  # comb_gene1 = FALSE
                               # outcome = outcome, lambda = lambda, alpha = alpha, nfolds = nfolds)
    
    final_freqs = c(final_freqs, GO_i$final_freqs)
    
    if(combine_haplos){
      gene1 = GO_i$lst$lst
    } else {
      gene1 = GO_i$lst
    }
    
    all_lists[[i + 1]] = gene1
    all_EMs_before[[i + 1]] = GO_i$EM_out_before
    all_EMs_after[[i + 1]] = GO_i$EM_out_after
    
    GOs[[i]] = GO_i$GOs[[1]]
    
    if(naive){
      all_naive_out[[(i + 1)]] = Naive_EMbased_analysis(lst = gene1, haplo = TRUE)  #, adj_comb_vals = GO_i$comb_vals)
    }
    
    
    other_genes[[GO_i$max]] = NULL
    
    gene_names = GO_i$gene_names
    orderings = GO_i$orderings
    comb_vals = GO_i$comb_vals
  }
  
  
  


  if(outcome_EM){
    start_HTFs = all_EMs_after[[nr_genes]]$Frequencies[nrow(all_EMs_after[[nr_genes]]$Frequencies), ]
    
    if(!(type_regression %in% c("all", "alleles", "backward", "penalized"))){
      cat("Type of regression models within the EM-algorithm is not well specified, only have: 'alleles', 'backward' or 'penalized'. Therefore we will run all three.\n")
      
      type_regression = "all"
    }
    
    if(type_regression == "all"){

      # lst = gene1; nl = 200; freq_epsilon = 1e-5; likelihood_estimation = FALSE; likelihood_epsilon = 0; init_equal = TRUE; equal_first = FALSE; status = NULL; nfolds = 50; split_evenly = TRUE; probYH_epsilon = 1e-10
      EM_out_Y1 = EM_algorithm(lst = gene1, regression = regression, full_regression = FALSE, type_regression = "alleles", regression_reference = regression_reference, start_HTFs = start_HTFs,
                               outcome = outcome, status = status, comb_vals = comb_vals, alpha = alpha, lambda = lambda, nfolds = nfolds, noCNV = TRUE, CNV_option = CNV_option, weighted_predictor = weighted_predictor)

      EM_out_Y2 = EM_algorithm(lst = gene1, regression = regression, full_regression = FALSE, type_regression = "backward", regression_reference = regression_reference, start_HTFs = start_HTFs,
                               outcome = outcome, status = status, comb_vals = comb_vals, alpha = alpha, lambda = lambda, nfolds = nfolds, noCNV = TRUE, CNV_option = CNV_option, weighted_predictor = weighted_predictor)

      EM_out_Y3 = EM_algorithm(lst = gene1, regression = regression, full_regression = TRUE, type_regression = "penalized", regression_reference = regression_reference, start_HTFs = start_HTFs,
                               outcome = outcome, status = status, comb_vals = comb_vals, alpha = alpha, lambda = lambda, nfolds = nfolds, noCNV = FALSE, CNV_option = CNV_option, weighted_predictor = weighted_predictor)
      
      EM_out_Y = list("analysis" = list("alleles" = EM_out_Y1, "backward" = EM_out_Y2, "penalized" = EM_out_Y3))
      
    } else {
      outcome_EM_one = EM_algorithm(lst = gene1, regression = regression, full_regression = full_regression, type_regression = type_regression, regression_reference = regression_reference,
                                    outcome = outcome, status = status, comb_vals = comb_vals, alpha = alpha, lambda = lambda, nfolds = nfolds, noCNV = noCNV, CNV_option = CNV_option, weighted_predictor = weighted_predictor)
      outcome_EM_one_lst = list(outcome_EM_one); names(outcome_EM_one_lst) = type_regression
      
      
      EM_out_Y = list("analysis" = outcome_EM_one_lst)
    }
    
    
    if(naive){
      naive_out_Y = Naive_EMbased_analysis(lst = gene1, haplo = TRUE)  #, adj_comb_vals = comb_vals)
    }
        
  } else {
    EM_out_Y = naive_out_Y = NULL
  }
  
  
  
  
  
  out = NULL
  
  out$EM_out = all_EMs_after
  out$EM_out_before = all_EMs_before
  out$comb_vals = comb_vals
  
  out$final_freqs = final_freqs
  out$Lists = all_lists
  out$GO_deci = GOs
  out$gene_names = gene_names
  out$orderings = orderings
  
  if(outcome_EM){
    out$EM_out_Y = EM_out_Y
    
    if(naive){
      out$Naive_out_Y = naive_out_Y
    }
  }
  if(naive){
    out$Naive_out = all_naive_out
  }
  
  
  return(out)
}
