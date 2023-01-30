#' @title Individual haplotype reconstruction measure
#' 
#' @description Multiple estimates of the individual haplotype reconstruction measure, where the individual's probability vector is compared to the truth. This haplotype reconstruction measure gives a better indication on how well the haplotype reconstruction methods perform on an individual level.
#' 
#' @details There are different haplotype reconstruction measures possible, depending on which arguments are specified. \itemize{
#' \item Compare the probability estimated for the correct diplotype with 1. Executed when \code{only_corr} is \code{TRUE}.
#' \item Compare the probabilities estimated for all compatible diplotypes, diplotypes that have alleles that are partly correct are considered to be correct according to ratio. Executed when \code{only_corr} is \code{FALSE} and \code{DT_error} is \code{TRUE}.
#' \item Compare the probabilities estimated for all compatible diplotypes, haplotypes that have alleles that are partly correct are considered to be correct according to ratio. Executed when \code{only_corr} is \code{FALSE} and \code{DT_error} is \code{FALSE}.
#' }
#'   The second option is considered to be the best measure to summarize the individual haplotype frequencies when no haplotype effects are considered. Diplotypes consisting of different haplotype configurations (which the algorithm cannot distinguish between) are not punished this way. 
#' 
#' @param EM_mat,naive_mat An optional matrix. Denoting for each individual the probabilities estimated for each possible diplotypes, as can be obtained via \code{EM-algorithm} or \code{Naive_analysis}. One of the two must be specified.
#' @param correct_mat A matrix. The matrix with the correct probabilities for each diplotype, thus each individual has one 1 for the actual simulated diplotype and a 0 for all other diplotypes.
#' @param measure1 A logical scalar, whether or not we are interested only in the completely correct DT or that we are interested in all diplotypes and consider the diplotypes that have alleles that are partly correct to be correct according to ratio (\code{TRUE} is default).
#' @param noCNV A logical scalar, whether or not we want to split the CNVs into separate alleles (\code{TRUE} is default).
#' 
# #' @param only_corr A logical scalar. Whether or not only the probability of the correct diplotype is compared (\code{FALSE} is default).
# #' @param DT_error A logical scalar. Whether or not an allelic mistake in the diplotype is wrong, or that an allelic mistake in the haplotype is wrong (\code{TRUE} is default). Only when \code{only_corr} if \code{FALSE}.
# #' @param HTFs An optional vector. 
#' 
#' @return The sum of squares measure of the chosen haplotype reconstruction measure.
#' 
HT_reconstruction_criteria = function(EM_mat = NULL, naive_mat = NULL, correct_mat, measure1 = TRUE, noCNV = TRUE){  
  # , only_corr = FALSE, DT_error = TRUE, HTFs = NULL){
  
  if(!is.null(EM_mat) & !is.null(naive_mat) | is.null(EM_mat) & is.null(naive_mat)){
    stop("Please supply only an EM matrix or only a naive analysis matrix")
  }
  if(!is.null(EM_mat)){
    test_mat = EM_mat
  } else if(!is.null(naive_mat)){
    test_mat = naive_mat
  }
  
  N = nrow(correct_mat)
  
  ##
  cnames = colnames(test_mat); cnames_corr = colnames(correct_mat)
  
  cnames_missing = cnames_corr[!(cnames_corr %in% cnames)]; len_cnames_missing = length(cnames_missing)
  test_mat = cbind(test_mat, matrix(0, nrow = N, ncol = len_cnames_missing, dimnames = list(NULL, cnames_missing)))
  
  cnames_missing = cnames[!(cnames %in% cnames_corr)]; len_cnames_missing = length(cnames_missing)
  correct_mat = cbind(correct_mat, matrix(0, nrow = N, ncol = len_cnames_missing, dimnames = list(NULL, cnames_missing)))
  
  cnames = colnames(test_mat); len_cnames = length(cnames); cnames = sort(cnames)
  test_mat = test_mat[, cnames]; correct_mat = correct_mat[, cnames]
  

  reconstruction_error = HT_criteria = NULL
  if(measure1){
    diff_mat = (test_mat - correct_mat)^2
    reconstruction_error = rowSums(diff_mat)
    
    
  } else {
    cnames_splt = lapply(strsplit(cnames, "\\+"), function(x){strsplit(x, "\\-")}); nr_genes = length(cnames_splt[[1]][[1]])
    
    reconstruction_error_sub = rep(list(NULL), nr_genes); reconstruction_error = NULL
    for(i in 1:nr_genes){
      correct_mat_i = reducing_mat(EM_mat = correct_mat, nr_gene = i, haplo = TRUE, noCNV = noCNV)
      test_mat_i = reducing_mat(EM_mat = test_mat, nr_gene = i, haplo = TRUE, noCNV = noCNV)
      
      reconstruction_error_sub[[i]] = rowSums((test_mat_i - correct_mat_i)^2)
    }
    
    for(i in 1:N){
      reconstruction_error[i] = mean(unlist(lapply(reconstruction_error_sub, function(x){x[i]})))
    }
  }
  
  
  HT_criteria = mean(reconstruction_error)
  
  #   cnames_splt = lapply(strsplit(cnames, "\\+"), function(x){strsplit(x, "\\-")})
  #   
  #   nr_genes = length(cnames_splt[[1]][[1]])
  #   reconstruction_error_sub = rep(list(rep(0, nr_genes)), N)
  #   
  #   for(i in 1:N){
  #     cnames_corr_splt = cnames_splt[[which(cnames == names(which(correct_mat[i, ] == 1)))]]
  #     
  #     test_index = which(cnames %in% names(which(test_mat[i, ] != 0))); len_test_index = length(test_index)
  #     cnames_test_splt = cnames_splt[test_index]; test_mat_sub = test_mat[i, test_index]  # ; diff_mat_sub = diff_mat[i, test_index]
  #     
  #     for(j in 1:nr_genes){
  #       corr_ij = unlist(lapply(cnames_corr_splt, function(x){x[[j]]}))
  #       
  #       for(k in 1:len_test_index){
  #         test_ijk = unlist(lapply(cnames_test_splt[[k]], function(x){x[[j]]}))
  #         
  #         if(noCNV){
  #           test_ijk = unlist(strsplit(test_ijk, "\\^"))
  #           corr_ijk = unlist(strsplit(corr_ij, "\\^"))
  #           
  #         } else {
  #           corr_ijk = corr_ij
  #         }
  #         
  #         len_test_ijk = length(test_ijk); len_corr_ijk = length(corr_ijk)
  #         temp = rep(TRUE, len_test_ijk); matching = 0
  #         for(l in 1:len_test_ijk){
  #           if(corr_ijk[l] %in% test_ijk[temp]){
  #             
  #             temp[which(temp & test_ijk == corr_ijk[l])[1]] = FALSE
  #             matching = matching + 1
  #           }
  #         }
  #         
  #         reconstruction_error_sub[[i]][j] = reconstruction_error_sub[[i]][j] + (test_mat_sub[k] * (matching / max(len_test_ijk, len_corr_ijk)))  # (diff_mat_sub[k] * (matching / len_test_ijk))
  #       }
  #     }
  #   }
  #   reconstruction_error = (unlist(lapply(reconstruction_error_sub, function(x){sum(x) / nr_genes})) - 1)^2
  #   
  #   ##
  #   reconstruction_error_sub = unlist(lapply(reconstruction_error_sub, function(x){sum(x) / nr_genes}))
  #   reconstruction_error = ((reconstruction_error_sub - 1)^2 + (1 - reconstruction_error_sub)^2)
  #   ##
  #   
  # }
  
  
  # collapsed_index = which(!(cnames_corr %in% cnames)); len_collapsed_index = length(collapsed_index)
  # if(!is.null(HTFs)){
  #   HTFs_names = sort(names(HTFs))
  #   
  #   HTFs = HTFs[HTFs_names]
  # }
  # 
  # if(only_corr){
  #   # if(is.null(HTFs)){
  #     if(len_collapsed_index == 0){
  #       for(i in 1:N){
  #         corr_index = names(which(correct_mat[i, ] == 1))
  #         
  #         reconstruction_error[i] = (test_mat[i, corr_index] - 1)^2
  #       }
  #       
  #     } else {
  #       for(i in 1:N){
  #         corr_index = names(which(correct_mat[i, ] == 1))
  #         
  #         if(corr_index %in% cnames_test){
  #           reconstruction_error[i] = (test_mat[i, corr_index] - 1)^2
  #         } else {
  #           reconstruction_error[i] = 1
  #         }
  #       }
  #       
  #       #    stop("There are DTs that have been (partly) collapsed, what to do...")  These are 
  #     }
  #     
## Weighted score, but that is not considered good/best so is discarded for now
#     } else {
#       DTs = t(outer(HTFs_names, HTFs_names, paste, sep = "+"))
#       DTs = DTs[lower.tri(DTs, diag = TRUE)] #  DTs[!upper.tri(DTs)]
#         
#       DTFs = t(outer(HTFs, HTFs, "*"))
#       DTFs[lower.tri(DTFs)] = DTFs[lower.tri(DTFs)] * 2
#       DTFs = DTFs[lower.tri(DTFs, diag = TRUE)] #  DTFs[!upper.tri(DTFs)]
#         
#       names(DTFs) = DTs
#       
#         
#       if(len_collapsed_index == 0){
#         for(i in 1:N){
# #          corr_index = names(which(correct_mat[i, ] == 1))
#           corr_index = cnames_corr[correct_mat[i, ] == 1]
# 
#           reconstruction_error[i] = (test_mat[i, corr_index] - 1)^2 * DTFs[DTs == corr_index]
#         }
#         
#       } else {
#         for(i in 1:N){
# #          corr_index = names(which(correct_mat[i, ] == 1))
#           corr_index = cnames_corr[correct_mat[i, ] == 1]
# 
# 
#           if(0 %in% cnames_test){
#             cat(i, "", test_mat[i, corr_index], "\n")
#             reconstruction_error[i] = (test_mat[i, corr_index] - 1)^2 * DTFs[DTs == corr_index]
# 
#           } else {
#             reconstruction_error[i] = 0 # Because the DT was collapsed, its DTFs is assumed to be 0 (?)
#           }
#         }
#         
#       }
#     }
  #   
  # } else {
  #   cnames_test_splt = lapply(strsplit(cnames_test, "\\+"), function(x){strsplit(x, "\\-")})
  #   cnames_corr_splt = lapply(strsplit(cnames_corr, "\\+"), function(x){strsplit(x, "\\-")})
  #   
  #   nr_genes = length(cnames_corr_splt[[1]][[1]])
  #   reconstruction_error_sub = rep(list(rep(0, nr_genes)), N)
  #   
  #   
  #   if(DT_error){  # Only when an allele in the DT is incorrect
  #     for(i in 1:N){
  #       corr_DTi = names(which(correct_mat[i, ] == 1)); corr_index = which(cnames_corr == corr_DTi)
  #       cnames_corr_splt_sub = cnames_corr_splt[[corr_index]]
  #       
  #       test_DTi = names(which(test_mat[i, ] != 0)); test_index = which(cnames_test %in% test_DTi); len_test_index = length(test_index)
  #       cnames_test_splt_sub = cnames_test_splt[test_index]
  #       test_mat_sub = test_mat[i, test_index]
  #       
  #       for(j in 1:nr_genes){
  #         corr_ij = unlist(lapply(cnames_corr_splt_sub, function(x){x[[j]]}))
  #         
  #         for(k in 1:len_test_index){
  #           test_ijk = unlist(lapply(cnames_test_splt_sub[[k]], function(x){x[[j]]}))
  #           test_ijk = unlist(strsplit(test_ijk, "\\^"))
  #           
  #           ##
  #           # temp = c(TRUE, TRUE); matching = 0
  #           # for(l in 1:2){
  #           #   if(test_ijk[l] %in% corr_ij[temp]){
  #           # 
  #           #     temp[which(corr_ij == test_ijk[l])[1]] = FALSE
  #           #     matching = matching + 1
  #           #   }
  #           # }
  #           # 
  #           # reconstruction_error_sub[[i]][j] = reconstruction_error_sub[[i]][j] + (test_mat_sub[k] * (matching / 2))   # (sum(test_ijk %in% corr_ij) / 2))
  #           ##
  #           ### OR
  #           ##
  #           temp = rep(TRUE, length(test_ijk)); matching = 0
  #           for(l in 1:2){
  #             if(corr_ij[l] %in% test_ijk[temp]){
  #               
  #               temp[which(temp & test_ijk == corr_ij[l])[1]] = FALSE
  #               matching = matching + 1
  #             }
  #           }
  #           
  #           reconstruction_error_sub[[i]][j] = reconstruction_error_sub[[i]][j] + (test_mat_sub[k] * (matching / length(test_ijk)))  
  #           ##
  #         }
  #       }
  #     }
  #     reconstruction_error = (unlist(lapply(reconstruction_error_sub, function(x){sum(x) / nr_genes})) - 1)^2
  #     
  #   } else {  # When an allele in the HT is incorrect
  #     
  #     for(i in 1:N){
  #       corr_DTi = names(which(correct_mat[i, ] == 1)); corr_index = which(cnames_corr == corr_DTi)
  #       cnames_corr_splt_sub = cnames_corr_splt[[corr_index]]
  #       
  #       test_DTi = names(which(test_mat[i, ] != 0)); test_index = which(cnames_test %in% test_DTi); len_test_index = length(test_index)
  #       cnames_test_splt_sub = cnames_test_splt[test_index]
  #       test_mat_sub = test_mat[i, test_index]
  #       
  #       for(j in 1:nr_genes){
  #         corr_ij = unlist(lapply(cnames_corr_splt_sub, function(x){x[[j]]}))
  #         
  #         for(k in 1:len_test_index){
  #           test_ijk = unlist(lapply(cnames_test_splt_sub[[k]], function(x){x[[j]]}))
  # 
  #           reconstruction_error_sub[[i]][j] = reconstruction_error_sub[[i]][j] + (test_mat_sub[k] * (sum(test_ijk == corr_ij) / 2))
  #         }
  #       }
  #     }
  #     reconstruction_error = (unlist(lapply(reconstruction_error_sub, function(x){sum(x) / nr_genes})) - 1)^2
  #     
  #   }
  # }
  # 
  # HT_criteria = sum(reconstruction_error) / N
  
  
  return(list("HT_criteria" = HT_criteria, "reconstruction_error" = reconstruction_error))
}

