#' @title Haplotype reconstruction via the profile EM-algorithm
#' 
#' @description General function to run the profile EM-algorithm as described in \cite{}.
#' 
#' @param gene_GTs An optional data frame. Each column contains the genotype calls of a separate gene. If \code{NULL} then \code{gene_lsts} must be supplied.
#' @param gene_lsts A optional list. Each list element has the compatible DTs, as can be provided by \code{\link{all_options}}. If \code{NULL} then \code{gene_GTs} must be supplied, if both are supplied, \code{gene_lsts} is overwritten by input of \code{gene_GTs}.
#' @param nr_start_genes An optional integer. The number of genes that should be used as a starting gene. The reconstruction is executed with the indicated number of starting genes denoted in the order of \code{gene_lsts}. As default each gene will be used as starting gene once.
#' @param gene_order A optional vector indicating the order in which the results of each reconstruction should be reordered.
#' @param CO_thresh A numerical value. Threshold value that determines which alleles or haplotypes are rare and will be grouped in the compound HT (\emph{$1e^{-07}$} is default).
#' @param all_previous A logical scalar. Whether or not all old combinational haplotypes should be forced to be grouped in the new combinational haplotype in a new reconstruction step if their probability is higher than \code{CO_thresh} (\code{TRUE} is default).
#'
#' @return \code{HT_reconstruction} returns a list with two elements: 
#' \describe{
#' \item{Matrices}{a matrix containing the estimated diplotype frequencies for each donor. The diplotypes of each starting gene are reordered into the common gene order}
#' \item{Frequencies}{the haplotype frequencies estimated for the final reconstructions. The haplotypes of each starting gene are reordered into the common gene order}
#' }
#' 
# #' @references 
#'
#' @examples
#' gene1 = list(c("A+A", "A+B", "B+B"), "C+D", c("A+D", "C+D"))
#' gene2 = list(c("E+E", "E+F"), c("E+E", "E+F", "F+F"), "E+F")
#' gene3 = list("G+H", "G+H", c("G+G", "G+H", "H+H"))
#'
#' gene_lsts = list("Gene1" = gene1, "Gene2" = gene2, "Gene3" = gene3)
#'
#' \dontrun{
#' HT_reconstruction(gene_lsts = gene_lsts, gene_order = names(gene_lsts))
#' HT_reconstruction(gene_lsts = gene_lsts, nr_start_genes = 1, gene_order = names(gene_lsts))
#' }
#' 
#' @export HT_reconstruction
HT_reconstruction = function(gene_GTs = NULL, gene_lsts = NULL, nr_start_genes = NULL, gene_order = NULL, CO_thresh = 1e-7, all_previous = TRUE,
                             outcome_EM = FALSE, regression = NULL, outcome = NULL){  #, CO_perc = 100, CO_min = 1, reduce_singles = FALSE, CO_singles = 1e-10){
  
  if(!is.null(gene_GTs)){
    nr_genes = ncol(gene_GTs); gene_lsts = rep(list(NULL), nr_genes); names(gene_lsts) = colnames(gene_GTs)
    for(i in 1:nr_genes){
      gene_lsts[[i]] = all_options(dat = gene_GTs[, 1])
    }
  } else if(!is.null(gene_lsts)){
    nr_genes = length(gene_lsts)
    
  } else {
    stop("Please supply either gene_GTs or gene_lsts, can't continue without data")
  }
  
  
  if(is.null(nr_start_genes)){
    nr_start_genes = nr_genes
    
  } else if(nr_start_genes > nr_genes){
    nr_start_genes = nr_genes
    
    warning("The number of starting genes was too large, it has been reduced to match the number of genes supplied in with the data.")
  }
  
  gene_names = names(gene_lsts)
  if(is.null(gene_order)){
    gene_order = gene_names
  }


  multiple_out = rep(list(NULL), nr_start_genes)
  names(multiple_out) = gene_names[1:nr_start_genes]
  
  Matrices = Rec_gene_names = Frequencies = rep(list(NULL), nr_start_genes)
  names(Matrices) = names(Rec_gene_names) = names(Frequencies) = gene_names[1:nr_start_genes]

  
  start_time = Sys.time()
  if(!outcome_EM){
    for(i in 1:nr_start_genes){
      cat("\n This iteration starts with ", gene_names[i], "\n")
      
      temp_gene_names = c(gene_names[i], gene_names[-i])
      one_reconstruction = grouping_order_all(gene1 = gene_lsts[[i]], other_genes = gene_lsts[-i], comb_gene1 = TRUE, gene_names = temp_gene_names, 
                                              CO_thresh = CO_thresh, all_previous = all_previous)
      
      
      multiple_out[[i]] = one_reconstruction
      
      Rec_gene_names[[i]] = one_reconstruction$gene_names
      
      Matrices[[i]] = one_reconstruction$EM_out[[nr_genes]]$Matrix
      Frequencies[[i]] = sort(one_reconstruction$EM_out[[nr_genes]]$Frequencies[nrow(one_reconstruction$EM_out[[nr_genes]]$Frequencies), ], decreasing = TRUE); names_Frequencies = names(Frequencies[[i]])
    }
    
    one_rec = NULL
  
  } else {
    
    analysis_lst = rep(list(NULL), nr_start_genes); names(analysis_lst) = gene_names[1:nr_start_genes]
    
    len_analysis = 3
    Betas_mod = Betas = Frequencies_Y = rep(list(analysis_lst), len_analysis); names(Betas_mod) = c("alleles", "backward", "penalized")
    
    for(i in 1:nr_start_genes){
      cat("\n This iteration starts with ", gene_names[i], "\n")
      
      temp_gene_names = c(gene_names[i], gene_names[-i])
      one_reconstruction = grouping_order_all(gene1 = gene_lsts[[i]], other_genes = gene_lsts[-i], comb_gene1 = TRUE, gene_names = temp_gene_names, 
                                              CO_thresh = CO_thresh, all_previous = all_previous,
                                              outcome_EM = outcome_EM, regression = regression, outcome = outcome)
      
      
      multiple_out[[i]] = one_reconstruction
      
      Rec_gene_names[[i]] = one_reconstruction$gene_names
      
      Matrices[[i]] = one_reconstruction$EM_out[[nr_genes]]$Matrix
      Frequencies[[i]] = sort(one_reconstruction$EM_out[[nr_genes]]$Frequencies[nrow(one_reconstruction$EM_out[[nr_genes]]$Frequencies), ], decreasing = TRUE); names_Frequencies = names(Frequencies[[i]])
      
      for(y in 1:len_analysis){
        nr_iter = nrow(one_reconstruction$EM_out_Y$analysis[[y]]$Frequencies)
        
        Betas_mod[[y]][[i]] = one_reconstruction$EM_out_Y$analysis[[y]]$betas_mod[nr_iter, ]
        Betas[[y]][[i]] = one_reconstruction$EM_out_Y$analysis[[y]]$betas[nr_iter, ]
        
        Frequencies_Y[[y]][[i]] = one_reconstruction$EM_out_Y$analysis[[y]]$Frequencies[nr_iter, ]
      }
    }
  }
  cat("\nThe complete analysis took", difftime(Sys.time(), start_time, units = "hours"), "hours (i.e.", difftime(Sys.time(), start_time, units = "mins"), "minutes or", difftime(Sys.time(), start_time, units = "secs"), "seconds) \n")
  
  Matrices = reordering_grouping_results(names_lst = Rec_gene_names, freqs_mat = Matrices, end_order = gene_order)
  Frequencies = reordering_grouping_results(names_lst = Rec_gene_names, freqs_lst = Frequencies, end_order = gene_order)
  
  
  if(outcome_EM){
    for(y in 1:len_analysis){
      Betas[[y]] = reordering_grouping_results(names_lst = Rec_gene_names, freqs_lst = Betas[[y]], end_order = gene_order)
      
      Frequencies_Y[[y]] = reordering_grouping_results(names_lst = Rec_gene_names, freqs_lst = Frequencies_Y[[y]], end_order = gene_order)
    }
  } else {
    Betas_mod = Betas = Frequencies_Y = NULL
  }

  
  
  return(list("Matrices" = Matrices, "Frequencies" = Frequencies, "All" = multiple_out, "Betas_mod" = Betas_mod, "Betas" = Betas, "Frequencies_Y" = Frequencies_Y))
}
