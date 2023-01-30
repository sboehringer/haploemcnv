#' @title Single reconstruction step for the profile EM-algorithm. 
#' 
#' @description Determines which gene is best suitable to be added into the reconstruction next. Analyses are run with this gene added into the reconstruction. Support function of \code{\link{grouping_order_all}},
#' 
#' @param gene1 A list. Diplotype list of the starting gene, can be obtained via \code{\link{all_options}}.
# #' @param other_genes An optional list. Diplotype list of the other genes that need to be added into the reconstruction, can be obtained via \code{\link{all_options}}. 
#' @param ... The diplotype lists of the other genes can also be supplied without specification of \code{other_genes}.
#' @param reference_info An optional list. Contains the frequencies and diplotypes in the combined group, required for determining the grouping order. If no grouping is required for the reference group, the `comb` needs to be specified as NULL. If no vector is supplied, the analysis of gene1 alone will be used.
#' @param gene_names An optional vector. Names of genes, need to be in the order in which they are supplied.
#' @param orderings An optional vector. Numerical order of the genes.
#' @param comb_vals An optional vector. The \emph{true} frequency of the combinational haplotypes is, thus at which values the \code{comb_vals} should be fixed in each iteration.
#' @param all_previous A logical scalar. Whether or not all previous compound haplotypes should be grouped in the newly formed compound haplotype. The name of all previous compound haplotypes should start with \emph{Comb}.
#' @param fixed_order A logical scalar. Whether or not the order in which the genes are grouped is based on the order in which they are supplied, or that the order is based on some heuristic strategy.
#' @param candidate_list An optional vector. Which haplotypes are never allowed to be grouped into the compound haplotype.
#' @param combine_haplos A logical scalar. Whether or not low-frequency haplotypes in the haplotype reconstruction are grouped (\code{TRUE} is default). 
#' @param CO_thresh An optional numeric value in \{0, 1\}. The cut-off threshold, only haplotypes lower than this threshold are be grouped.
#' @param CO_perc An optional numerical value in \{0, 100\}. The cut-off percentage, only the haplotypes belonging to the percentage of haplotypes with the lowest frequency are eligible for grouping.
#' @param CO_min An optional integer. The cut-off minimum, the minimum number of haplotypes which need to be grouped to actually group haplotypes.
#' 
#' @return The output of \code{\link{grouping_order}} for each of the grouping configurations with the starting gene and the gene added into the reconstruction (with analysis).
#'
# #' @examples 
# #' gene1 = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' gene2 = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' gene3 = list(c("007+008", "007+009"), "NEG+NEG", c("008+008", "008+NEG"))
# #' 
# #' \dontrun{
# #' grouping_order_iter(gene1, gene2, gene3)
# #' grouping_order_iter(gene1, list(gene2, gene3))
# #' }
#' 
grouping_order_iter = function(gene1, ..., reference_info = NULL, gene_names = NULL, orderings = NULL, comb_vals = NULL, all_previous = FALSE, fixed_order = FALSE, candidate_list = NULL, 
                                combine_haplos = TRUE, CO_thresh = 1e-7, CO_perc = 100, CO_min = 2){  # comb_gene1 = FALSE
  
  if(is.list(c(...))){
    other_genes = c(...)
  } else {
    other_genes = list(...)
  }
  
  
  if(fixed_order == FALSE){
    nr_other_genes = length(other_genes)
  } else {
    nr_other_genes = 1
  }
  
  if(is.null(gene_names)){
    gene_names = 1:(nr_other_genes + 1)
  }
  
  if(is.null(orderings)){
    orderings = 1:length(gene_names)
  }
  
  
  reference_comb = reference_info$comb
  reference_freqs = reference_info$freqs

  
  gene1_secC = as.list(rep(0, nr_other_genes))
  out_gene1_sec = out_gene1_secC = as.list(rep(0, nr_other_genes))
  
  GOs = as.list(rep(0, nr_other_genes))
  all_GO_diffs = NULL
  
  for(i in 1:nr_other_genes){
    sec_gene = other_genes[[i]]
    
    gene1_sec = combining_genes(gene1, sec_gene, haplo = TRUE)

    if(combine_haplos == TRUE){
      out_gene1_sec[[i]] = EM_algorithm(lst = gene1_sec, comb_vals = comb_vals, regression = NULL)  # Before

      comb_name = paste("CombG", orderings[1], orderings[1 + i], sep = "")
      gene1_secC[[i]] = combining_AF(lst = gene1_sec, EM_out = out_gene1_sec[[i]], haplo_grouping = TRUE, all_previous = all_previous,
                                      CO_thresh = CO_thresh, CO_perc = CO_perc, CO_min = CO_min, comb_name = comb_name, excluding_haplo = candidate_list)
      
      if(length(gene1_secC[[i]]$comb) != 0){
        comb_vals = c(comb_vals, comb_frequency(haplotypes = gene1_secC[[i]]$haplotypes, remain = gene1_secC[[i]]$remain))
        names(comb_vals) = c(names(comb_vals)[-length(comb_vals)], comb_name)
      }
      

      un_haplos = unique(unlist(strsplit(unique(unlist(gene1_secC[[i]]$lst)), "\\+")))
      comb_haplos = un_haplos[str_detect(un_haplos, "Comb")]

            
      if(length(comb_haplos) == 0){
        out_gene1_secC[[i]] = EM_algorithm(lst = gene1_secC[[i]]$lst, regression = NULL)

      } else if(length(comb_haplos) == 1){
        out_gene1_secC[[i]] = EM_algorithm(lst = gene1_secC[[i]]$lst, comb_vals = comb_vals, regression = NULL)  # After

      } else {
        cat("In the combination of ", gene_names[1], ", with ", gene_names[i + 1], " there are multiple comb haplotypes:", comb_haplos, "\n")
        out_gene1_secC[[i]] = EM_algorithm(lst = gene1_secC[[i]]$lst, comb_vals = comb_vals, regression = NULL)
      }


    } else {
      gene1_secC[[i]] = gene1_sec
      
      out_gene1_secC[[i]] = EM_algorithm(gene1_secC[[i]], regression = NULL)
    }
    
    if(!is.null(reference_comb)){
      GOs[[i]] = grouping_order(freqs1 = reference_freqs, freqs2 = Haplo2AF(out_gene1_secC[[i]])[[1]], comb_group = reference_comb)
    } else {
      GOs[[i]] = grouping_order(freqs1 = reference_freqs, freqs2 = Haplo2AF(out_gene1_secC[[i]])[[1]], freqs1_no_comb = TRUE)
    }
    
    all_GO_diffs[i] = GOs[[i]]$Diff
  }

  max_index = which(all_GO_diffs == max(all_GO_diffs))
  
  if(length(max_index) != 1){
    cat("\n", "Both genes", gene_names[max_index + 1], "gives the highest difference in combination with", gene_names[1], "continue with gene", gene_names[max_index[1] + 1], "\n")
    
    max_index = max_index[1]
  } else {
    cat("\n", "The combination of", gene_names[1], "with gene", gene_names[max_index + 1], "gives the highest difference", "\n")
  }
  
  if(combine_haplos == TRUE){
    final_freqs = sort(out_gene1_secC[[max_index]]$Frequencies[nrow(out_gene1_secC[[max_index]]$Frequencies), ], decreasing = TRUE)
    
  } else {
    final_freqs = NULL
    
  }
  
  
  gene_names[1] = paste(gene_names[c(1, (max_index + 1))], collapse = "-")
  gene_names = gene_names[-(max_index + 1)]
  
  orderings[1] = paste(orderings[1], orderings[max_index + 1], sep = "")
  orderings = orderings[-(max_index + 1)]
  
  
  out = list("lst" = gene1_secC[[max_index]], "EM_out_before" = out_gene1_sec[[max_index]], "EM_out_after" = out_gene1_secC[[max_index]], "GOs" = GOs, "max" = max_index, 
              "final_freqs" = final_freqs, "gene_names" = gene_names, "orderings" = orderings, "comb_vals" = comb_vals)
  
  return(out)
}
