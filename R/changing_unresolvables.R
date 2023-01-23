#' @title Combining alleles that cannot be distinguished from each other in the genotyping process
#' 
#' @description To remove allelic ambiguities that originates from the fact that certain exons are not included in the sequencing process.
#' 
#' @details For consistency, if only a part of the unresolvable alleles is present in the individual's genotype call these alleles are still replaced by the indicated code.
#' 
#' @param dat A data frame. The genotype calls of each individual as a separate row and each gene as a separate column.
#' @param unresolvable A three-column data frame. The first column must contain the corresponding gene; the second column the allele combination that cannot be distinguished and the third column the code to which that combination needs to be changed
#'
#' @return A vector with the unresolvable allele combinations combined into the indicated code.
#' 
# #' @examples
# #' dat = matrix(c("001", "001/002", "003", "003", "005", "003/004"), nrow = 3, byrow = FALSE, dimnames = list(NULL, c("Gene1", "Gene2")))
# #' unresolvable = matrix(c("Gene1",	"001/002", "001c", "Gene2", "003/004",	"003c"), nrow = 2, byrow = TRUE)
# #' 
# #' changing_unresolvables(dat, unresolvable)
#' 
changing_unresolvables = function(dat, unresolvable){
  len_dat = nrow(dat); cnames = colnames(dat)
  len_unresolvable = nrow(unresolvable)
  
  dat_out = dat
  for(i in 1:len_unresolvable){
    
## Selec the gene, combination and code from the i-th unresolvable allele combination.
    gene = as.character(unresolvable[i, 1])
    combination = as.character(unresolvable[i, 2])
    code = as.character(unresolvable[i, 3])
    
    
## Because there is a 2DS4 and a 2DS4N, need the 2DS4.
    col = which(str_detect(cnames, gene))[1]  
    comb = unlist(strsplit(combination, "\\/"))
    
    for(j in 1:len_dat){
      dat_j = dat_out[j, col]
      
## If individual allele is part of the unresolvable combination, it is replaced.      
      if(!is.na(dat_j)){
        for(k in 1:length(comb)){
          dat_j = str_replace_all(dat_j, comb[k], code)
        }
        
        if(!str_detect(dat_j, "\\|")){
          if(!str_detect(dat_j, "\\/")){
            dat_out[j, col] = dat_j
            next
          }
          
          if(str_detect(dat_j, "\\/")){
            if(!str_detect(dat_j, "\\+")){
              splt = unlist(strsplit(dat_j, "\\/"))
              
              dat_out[j, col] = paste(sort(unique(splt)), collapse = "/")
              next
            }
            if(str_detect(dat_j, "\\+")){
              splt = unlist(strsplit(dat_j, "\\+"))
              
              for(k in 1:length(splt)){
                splttd = unlist(strsplit(splt[k], "\\/"))
                splt[k] = paste(sort(unique(splttd)), collapse = "/")
              }
              
              dat_out[j, col] = paste(sort(splt), collapse = "+")
              next
            }
          }
        }
        
        if(str_detect(dat_j, "\\|")){
          if(!str_detect(dat_j, "\\/")){
            splt = unlist(strsplit(dat_j, "\\|"))
              
            dat_out[j, col] = paste(sort(unique(splt)), collapse = "|")  
            next
          }
          
          if(str_detect(dat_j, "\\/")){
            if(!str_detect(dat_j, "\\+")){
              splt = unlist(strsplit(dat_j, "\\|"))
              
              for(k in 1:length(splt)){
                splttd = unlist(strsplit(splt[k], "\\/"))
                splt[k] = paste(sort(unique(splttd)), collapse = "/")
              }
              
              dat_out[j, col] = paste(sort(unique(splt)), collapse = "|")
              next
            }
            if(str_detect(dat_j, "\\+")){
              splt = unlist(strsplit(dat_j, "\\|"))
              
              for(k in 1:length(splt)){
                splttd = unlist(strsplit(splt[k], "\\+"))
                
                for(l in 1:length(splttd)){
                  splitted = unlist(strsplit(splttd[l], "\\/"))
                  splttd[l] = paste(sort(unique(splitted)), collapse = "/")
                }
                
                splt[k] = paste(sort(splttd), collapse = "+")
              }
              
              dat_out[j, col] = paste(sort(unique(splt)), collapse = "|")
              next
            }
          }
        }
      }
    }
  }    
  
  
  return(dat_out) 
}
