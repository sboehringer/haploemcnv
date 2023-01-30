#' @title Removal of specific strange symbols
#'
#' @description Removing certain strange symbols (`*`, `?`, `:`) that originate in the genotype calls from altering the initial \emph{KIR} dataset.
#' 
#' @details Function is tailored to the occurrences in the obtained \emph{KIR} dataset. For extensions to other dataset always look at the patterns observed in that data and check if they match. \itemize{
#'   \item Some gene copies are preceded by a `*`.
#'   \item Some single gene genotype calls are followed by `:X?`.
#'   \item Some multi-gene genotype calls are followed by `:X?`.
#' }
#' 
#' @param dat A vector. The genotype calls of a single gene. 
#' 
#' @return A vector where all the strange symbols are discarded from the genotypes.
#'
# #' @examples
# #' \dontrun{remove_strange_symbols("004+001:X?")}
#' 
remove_strange_symbols = function(dat){
  len_dat = length(dat)
  
  out = dat
  for(i in 1:len_dat){
    xx = dat[i]
    
    if(is.na(xx) == TRUE){
      cat(i)
      next
    }
    
    
## There are three cases which have unique strange symbols, special lines to adjust those 
    if(str_detect(xx, "\\*")){
      xx = unlist(strsplit(xx, "\\*"))[2]
    }
    
    if(str_detect(xx, "^\\:X?+")){
      out[i] = unlist(strsplit(xx, "\\+"))[2]
      next
    }
    
    if(str_detect(xx, "copy")){
      test = unlist(strsplit(xx, ""))
      
      out[i] = paste(test[1], test[2], test[3], collapse = "", sep = "")
      next
    }
    
    
## Select donors with ":" and with/without ambiguity, index only on the first, correct, argument
    if(str_detect(xx, ":")){
      if(!str_detect(xx, "\\/")){
        if(!str_detect(xx, "\\+")){
          test = strsplit(xx, ":")[[1]][1] 
          
          out[i] = test
          next
        }
      }
      
      
## the ":" only occurs after an allele with a "+"
      if(str_detect(xx, "\\+")){
        test = strsplit(xx, "\\+")[[1]]
        test = strsplit(test, ":")
        
        lntest = length(test)
        

        spltd = as.vector(rep(NA, length = lntest))
        for(j in 1:lntest){
          spltd[j] = test[[j]][1]
        }
        
        out[i] = paste(spltd, collapse = "+")
        next
      }
    }
    
    
## Select donors with a "?" and index only on the first, correct, argument
    if(str_detect(xx, "\\?")){
      test = unlist(strsplit(xx, "X"))
      
      out[i] = test[1]
      next
    }
  }
  
  return(out)
}