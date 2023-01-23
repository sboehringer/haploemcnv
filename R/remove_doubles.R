#' @title Removes redundant genotypic ambiguities
#' 
#' @description Remove genotypic ambiguities are double and/or redundant and places the alleles of the remaining ambiguities in the correct order for further processing.
#' 
#' @param dat A vector. The genotyping data of a single gene.
#'
#' @return a vector with the unique genotypic ambiguities in the correct order.
#' 
# #' @examples 
# #' \dontrun{remove_doubles("002+006|006+002")}
#' 
remove_doubles <- function(dat){
  len_dat <- length(dat)
  
  dat_out <- NULL
  for(i in 1:length(dat)){
    xx <- dat[i]
    
## Skip all donors which are "NEG"
    if(is.na(xx)){
      dat_out[i] <- "NEG"
      next
    }
    
## Select donors with a "|" and with/without a "+"
    if(str_detect(xx, "\\|")){
      if(!str_detect(xx, "\\+")){
        splt <- unlist(strsplit(xx, "\\|"))
        
        dat_out[i] <- paste(sort(unique(splt)), collapse = "|")
        next
      }
      
      if(str_detect(xx, "\\+")){
        splt <- unlist(strsplit(xx, "\\|"))
        
        for(j in 1:length(splt)){
          splt[j] <- paste(sort(unlist(strsplit(splt[j], "\\+"))), collapse = "+")
        }
        
        dat_out[i] <- paste(sort(unique(splt)), collapse = "|")
        next
      }
    }
    
## Select donors without a "|" and with a "+
    if(!str_detect(xx, "\\|")){
      if(str_detect(xx, "\\+")){
        splt <- unlist(strsplit(xx, "\\+"))
        
        dat_out[i] <- paste(sort(splt), collapse = "+")
        next
      }
    }
    
    dat_out[i] <- xx
  }
  
  return(dat_out)
}