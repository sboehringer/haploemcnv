#' @title Combining the diplotypes of two genes 
#'
#' @description To combine the information of two gene(s) combination(s), which is required for haplotype reconstruction based on multiple loci. Support function for \code{\link{HT_reconstruction}}.
#'
#' @details If the genes are combined on haplotype level, the diplotypes notation is as follows: \emph{x1-y1+x2-y2} and \emph{x1-y2+x2-y1}. If the genes are combined on genotypic level, the notation is: \emph{x1+x2&y1+y2}.
#'
#' @param x,y A list. The diplotype or haplotype list of a gene or of a combination of genes, can be obtained via \code{\link{all_options}}.
#' @param haplo a logical scalar. Whether or not \code{x} and \code{y} are combined on haplotype level, or are combined on genotype level (\code{TRUE} is default).
#' @param un_diplos a logical scalar. Whether or not only the unique diplotypes should be retained, or that each diplotype should occur the number of times it is created (\code{TRUE} is default).
#' 
#' @return A list with the combined diplotypes of the two specified genes. The output in a non-ambiguous case is like this: GeneA1+GeneA2&GeneB1+GeneB2
#' 
# #' @examples 
# #' x = list(c("001+001", "001+002", "002+002"), "003+NEG", c("001+NEG", "003+NEG"))
# #' y = list("006+NEG", c("004+NEG", "006+NEG"), c("004+004", "004+005", "005+005"))
# #' 
# #' \dontrun{
# #' combining_genes(x, y)
# #' }
combining_genes = function(x, y, haplo = TRUE, un_diplos = TRUE){
  len_comb = length(x)
  
  comb_dat = list(NULL)
  if(haplo == TRUE){
    for(i in 1:len_comb){

      xx = x[[i]]
      len_xx = length(xx)
      
      yy = y[[i]]
      len_yy = length(yy)
      
      haplo_all = NULL
      for(j in 1:len_xx){
        xxx = unlist(strsplit(xx[j], "\\+"))
        
        for(k in 1:len_yy){
          yyy = unlist(strsplit(yy[k], "\\+"))
          
          dat = list(xxx, yyy)
          grid.dat = as.matrix(expand.grid(dat))
          
          len_grid = nrow(grid.dat)
          len_grid_half = len_grid / 2
          for(l in 1:len_grid_half){
            haplo_all = c(haplo_all, paste(sort(c(paste(grid.dat[l, ], collapse = "-"), 
                                                  paste(grid.dat[(len_grid + 1 - l), ], collapse = "-"))), collapse = "+"))
          }
        }
      }
      
      if(un_diplos){
        comb_dat[[i]] = unique(sort(haplo_all))
      } else {
        comb_dat[[i]] = sort(haplo_all)
      }
    }
    
  } else {
    for(i in 1:len_comb){
      xx = x[[i]]
      yy = y[[i]]
      
      dat = list(xx, yy)
      grid.dat = as.matrix(expand.grid(dat))
      
      len_grid = nrow(grid.dat)
      haplo_all = NULL
      
      for(j in 1:len_grid){
        haplo_all = c(haplo_all, paste(grid.dat[j, ], collapse = "&"))
      }
      
      if(un_diplos){
        comb_dat[[i]] = unique(sort(haplo_all))
      } else {
        comb_dat[[i]] = sort(haplo_all)
      }
    }
  }

  
  return(comb_dat)
}
