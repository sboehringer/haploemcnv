#' @title Combining the KIR2DS4 and KIR2DS4N KIR genes into one gene
#' 
#' @description To combine the two KIR2DS4 genes into one. The gene copies of the two genes are simply added together.
#' 
#' @details If any of the two genes has a `POS` genotype call, than the combined gene also gets the `POS` genotype.
#' 
#' @param KIR2DS4 A vector. The genotyping data of KIR2DS4.
#' @param KIR2DS4N A vector. The genotyping data of KIR2DS4N.
#'
#' @return A vector with the combined genotyping result of the two genes. The combination made for each individual is also printed in the Rstudio console. 
#' 
# #' @examples 
# #' KIR2DS4 = c("NEG", "001", "NEG")
# #' KIR2DS4N = c("004", "003", "NEG")
# #' 
# #' \dontrun{
# #' KIR2DS4_combining(KIR2DS4, KIR2DS4N)
# #' }
#' 
KIR2DS4_combining = function(KIR2DS4, KIR2DS4N){
  
  N = length(KIR2DS4N)
  
  KIR2DS4_comb = NULL
  for(i in 1:N){
    x1 = KIR2DS4[i]; x2 = KIR2DS4N[i]
    
    if(x1 == "NEG" & x2 == "NEG"){
      KIR2DS4_comb[i] = "NEG"
      next
    } else if(x1 == "NEG"){
      KIR2DS4_comb[i] = x2
      next
    } else if(x2 == "NEG"){
      KIR2DS4_comb[i] = x1
      next
    }
    
    if(x1 == "POS" | x2 == "POS"){  # If any of the two are POS, then both are POS?
      KIR2DS4_comb[i] = "POS"
      next
    }
    
    x1_plus = str_detect(x1, "\\+")
    x1_slash = str_detect(x1, "\\/")
    x1_bar = str_detect(x1, "\\|")
    
    x2_plus = str_detect(x2, "\\+")
    x2_slash = str_detect(x2, "\\/")
    x2_bar = str_detect(x2, "\\|")
    
    if(x1_bar == FALSE & x2_bar == FALSE){
      KIR2DS4_comb[i] = paste(sort(c(x1, x2)), collapse = "+")
      next
    }

    if(x1_bar == TRUE | x2_bar == TRUE){
      splt1 = unlist(strsplit(x1, "\\|"))
      splt2 = unlist(strsplit(x2, "\\|"))
      
      KIR2DS4_comb[i] = paste(sort(c(apply(expand.grid(splt1, splt2), 1, function(x){paste(sort(c(x[1], x[2])), collapse = "+")}))), collapse = "|")
    }
    
    cat("For individual", i, " ", x1, " and ", x2, " become ", KIR2DS4_comb[i], "\n")
  }
  
  return(KIR2DS4_comb)
}
