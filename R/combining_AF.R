#' @title Grouping of the haplotypes based on their estimated frequencies
#' 
#' @description Grouping of haplotypes into a common compound haplotype. There are multiple grouping criteria possible that can be chosen \enumerate{
#'   \item Based on a minimum frequency the haplotype should have (\code{CO_thresh}).
#'   \item Based on a number of haplotypes that should be kept (\code{cumulative_n}).
#'   \item Based on the cumulative frequency that should be kept (\code{cumulative_threshold}).
#' }. It is strongly advised to make use of the first criteria (\code{CO_thresh}) since this is for general use the most logical one. Further alterations are possible with other arguments.
#' 
#' @details With \code{CO_perc} and \code{CO_min} one can add extra limitations to the \code{CO_thresh} grouping criteria, but in practice their gain is limited and they primarily complicate the understanding of the grouping 
#' 
#' @param lst A List. The data containing the possible diplotypes of each donor as a separate list element, can be obtained via \code{\link{all_options}}.
#' @param EM_out An optional EM-algorithm analysis of \code{lst}, as provided by \code{\link{EM_algorithm}}.
#' @param haplotypes An optional vector. The haplotypes with their frequency, can be obtained by \code{\link{EM_algorithm}}.
#' @param haplo_grouping A logical scalar. Whether or not grouping occurs on the haplotype level (\code{TRUE} is default) or on the allele level. In the case of 1 gene, haplotype and allele level are the same. 
#' @param comb_name An optional character. Indicating the name of the combined group. If \code{NULL} the name is either `Comb` or `Comb` followed by a number (if previous combining step have been performed for \code{lst}).
#' @param cumulative_threshold An optional numeric value. The cumulative threshold value above which you want to group all alleles. 
#' @param cumulative_n An optional integer. The number of haplotypes you want to retain, the \code{cumulative_n} haplotypes with the highest frequencies are retained.
#' @param CO_thresh An optional numeric value in \{0, 1\}. The cut-off threshold, only haplotypes lower than this threshold are be grouped.
#' @param CO_perc An optional numerical value in \{0, 100\}. The cut-off percentage, only the haplotypes belonging to the percentage of haplotypes with the lowest frequency are eligible for grouping.
#' @param CO_min An optional integer. The cut-off minimum, the minimum number of haplotypes which need to be grouped to actually group haplotypes.
#' @param excluding_haplo An optional vector. Which haplotypes never should be grouped, even if they do not meet the chosen criteria. NOTETOSELF: WAS excluding_geno, CAN BE WRONG IN OTHER FUNCTIONS
#' @param all_previous A logical scalar. Whether or not all previous compound haplotypes should be grouped in the newly formed compound haplotype. The name of all previous compound haplotypes should start with \emph{Comb}.
#' @param gene_subset An optional vector. Which genes need to be grouped if grouping_haplo = FALSE. NOTETOSELF: Volgens mij kan ik alles wat haplo_grouping == FALSE is wel weghalen, wordt denk ik niet meer gebruikt en weet ook niet wat ik er doe...}
#' 
#' @return A list with all possible genotypes of each donor, but now with some alleles replaced by `Comb`. Additionally are also the alleles which are combined into the `Comb` group provided.
#' 
# #' @examples 
# #' lst = list(c("001+001", "002+002"), "001+NEG", "003+NEG", c("001+NEG", "003+NEG"), "001+003", c("001+003", "004+004"), "001+001", "003+003")
# #' \dontrun{
# #' combining_AF(lst, haplo_grouping = TRUE, CO_thresh = 1e-7)
# #' }
#' 
combining_AF = function(lst, EM_out = NULL, haplotypes = NULL, haplo_grouping = TRUE, comb_name = NULL, cumulative_threshold = NULL, cumulative_n = NULL, CO_thresh = 1e-7, CO_perc = 100, CO_min = 2, 
                        excluding_haplo = NULL, all_previous = TRUE, gene_subset = NULL){
  
  len_lst = length(lst)
  
  if(is.null(haplotypes)){
    if(is.null(EM_out)){
      EM_out = EM_algorithm(lst = lst)
    }
    haplotypes = EM_out$Freq[nrow(EM_out$Freq), ]
  }

  
  nr_alleles = length(haplotypes)
  nr_genes = unique(lengths(strsplit(names(haplotypes), "\\-")))  
  if(length(nr_genes) != 1){
   cat("Something has gone wrong, not all haplotypes consist of the same amount of alleles \n")
  }
  
  
  un_alleles = unique(unlist(strsplit(unique(unlist(lst)), "[+-]+")))
  if(is.null(comb_name)){
    if(TRUE %in% str_detect(un_alleles, "Comb")){
      index = which(str_detect(un_alleles, "Comb"))
      
      if(TRUE %in% (un_alleles == "Comb")){
        comb_name = "Comb" 
      } else{    
        nr_combines = max(as.numeric(str_sub(un_alleles[index], start = 5, end = 5)))
        comb_name = paste(rep(paste("Comb", nr_combines + 1, sep = ""), nr_genes), collapse = "-")
      }
    } else {
      comb_name = paste(rep("Comb1", nr_genes), collapse = "-")
    }
  } else {
    comb_name = paste(rep(comb_name, nr_genes), collapse = "-")
  }
  
  
  
  if(haplo_grouping){
    sorted_haplos = sort(haplotypes, decreasing = TRUE); nr_haplos = length(sorted_haplos); sorted_names = names(sorted_haplos)
    
    cumsum_thresh = cumsum(sorted_haplos)
    if(!is.null(cumulative_threshold)){
      cat("Note that this does not work for very small probabilities, because the cumulative sum is 1 before those small probabilities are included \n")
      cum_thresh = cumsum_thresh < cumulative_threshold
      comb = sorted_names[which(cum_thresh == FALSE)[-1]]
      
    } else if(!is.null(cumulative_n)){
      if(nr_haplos > cumulative_n){
        cum_thresh = cumsum_thresh[1:cumulative_n]
        comb = sorted_names[!(sorted_names %in% names(cum_thresh))]
      } else {
        stop("You already have the same or less genotypes remaining")
      }
      
    } else {
      passing_threshold = sorted_haplos[sorted_haplos < CO_thresh]
      names_passing_threshold = names(passing_threshold); nr_passing_threshold = length(names_passing_threshold)
      
      nr_passing_percentage = floor(nr_haplos * CO_perc / 100)
      passing_percentage = tail(sorted_haplos, nr_passing_percentage); names_passing_percentage = names(passing_percentage)
      
      comb = names_passing_threshold[which(names_passing_threshold %in% names_passing_percentage)]; len_comb = length(comb)
      
      if(len_comb < CO_min){
        cat("There are ", len_comb, " haplotypes to combine, there should have been more than ", CO_min, ", so continue without combining \n")

        comb = NULL
        len_comb = 0
      }
    }
    
    
    if(!is.null(excluding_haplo)){
      
      len_excluding_haplo = length(excluding_haplo)
      for(i in 1:len_excluding_haplo){
        if(excluding_haplo[i] %in% comb){
          index = which(excluding_haplo %in% comb)
          
          comb = comb[!(comb %in% excluding_haplo[index])]
        }
      }
    }
    remain = sorted_names[!(sorted_names %in% comb)]
  
    if(all_previous == TRUE){
      comb = c(comb, remain[str_detect(remain, "Comb")])
      remain = remain[!str_detect(remain, "Comb")]
    }

  
    
  } else {
    freq_names = as.character(names(haplotypes))
    freq_name_splt = strsplit(freq_names, "\\-")
    
    if(!is.null(gene_subset)){
      thresh_subset = gene_subset
      len_thresh_subset = length(thresh_subset)
    } else {
      thresh_subset = 1:nr_genes
      len_thresh_subset = length(thresh_subset)
    }
    
    
    all_combs = as.list(rep(0, nr_genes))
    all_remains = as.list(rep(0, nr_genes))
    allel_freqs = Haplo2AF(haplotypes = haplotypes)
    
    for(i in thresh_subset){
      sorted_alleles = NULL
      sorted_alleles = sort(allel_freqs[[i]], decreasing = TRUE)
      names_alleles = names(sorted_alleles)
      len_names_alleles = length(names_alleles)
      
      if(TRUE %in% str_detect(names_alleles, "Comb")){
        cat("Gene ", i, " has previously already been grouped...")
      }
      
      cumsum_thresh = cumsum(sorted_alleles)
      
      if(!is.null(cumulative_threshold)){
        cat("Note that this doesnot work for very small probabilities, because the cumulative sum is 1 before those small probabilities are included \n")
        cum_thresh = cumsum_thresh < cumulative_threshold
        all_combs[[i]] = names_alleles[which(cum_thresh == FALSE)[-1]]
        
      } else if(!is.null(cumulative_n)){
        if(len_names_alleles > cumulative_n){
          cum_thresh = cumsum_thresh[1:cumulative_n]
          all_combs[[i]] = names_alleles[!(names_alleles %in% names(cum_thresh))]
        } else {
          cat("You already have the same or less genotypes remaining \n")
          all_combs[[i]] = 0
        }
        
      } else {
        passing_threshold = sorted_alleles[sorted_alleles < CO_thresh]
        names_passing_threshold = names(passing_threshold)
        nr_passing_threshold = length(names_passing_threshold)
        
        nr_passing_percentage = floor(len_names_alleles * CO_perc / 100)
        passing_percentage = tail(sorted_alleles, nr_passing_percentage)
        names_passing_percentage = names(passing_percentage)
        
        all_combs[[i]] = names_passing_threshold[which(names_passing_threshold %in% names_passing_percentage)]
        
        if(length(all_combs[[i]]) < CO_min){
          cat("There are less than ", CO_min, " haplotypes to group for gene ", i, " so grouping isnt executed \n")
          all_combs[[i]] = 0
        }
      }
      
      if(!is.null(excluding_haplo)){
        if(TRUE %in% (excluding_haplo %in% all_combs[[i]])){
          index = which(excluding_haplo %in% all_combs[[i]])
          
          all_combs[[i]] = all_combs[[i]][!(all_combs[[i]] == excluding_haplo[index])]
        }
      }
      
      all_remains[[i]] = names_alleles[!(names_alleles %in% all_combs[[i]])]
      
    }
  }
  


  
  if(haplo_grouping == TRUE){
    for(i in 1:len_lst){
      splt = strsplit(lst[[i]], "\\+")
      len_splt = length(splt)
      
      for(j in 1:len_splt){
        splttd = unlist(splt[[j]])
        len_splttd = length(splttd)
        
        for(k in 1:len_splttd){
          if(splttd[k] %in% comb){
            splttd[k] = comb_name
          }
        }  
        
        splt[[j]] = paste(unlist(sort(splttd)), collapse = "+")
      }
      
      lst[[i]] = unique(unlist(splt))  # Added unique(...) here so that each option only occurs once
    }
    
    
  } else {
    for(i in 1:len_lst){
      splt = strsplit(lst[[i]], "\\+")
      len_splt = length(splt)
      
      for(j in thresh_subset){
        for(k in 1:len_splt){
          splttd = strsplit(splt[[k]], "\\-")
          len_splttd = length(splttd)
          
          for(l in 1:len_splttd){
            
            if(splttd[[l]][j] %in% all_combs[[j]]){
              splttd[[l]][j] = comb_name
            }
          }
          
          splt[[k]] = unlist(lapply(splttd, function(x){paste(sort(x), collapse = "-")}))  # sort added
        }
      }
      
      lst[[i]] = unlist(lapply(splt, function(x){paste(sort(x), collapse = "+")}))
    }
  }

  
  if(haplo_grouping == TRUE){
    out = list("lst" = lst, "comb" = comb, "remain" = remain, "haplotypes" = haplotypes)
  } else {
    out = list("lst" = lst, "comb" = all_combs, "remain" = all_remains, "haplotypes" = haplotypes)
  }
  
  return(out)
}
