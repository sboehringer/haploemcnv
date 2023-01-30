#' @title Denoting all compatible diplotype combinations per individual
#' 
#' @description \code{all_options} is be used to transform genotype data of a single gene into a list containing all compatible diplotypes for that gene.
#' 
#' @details Different forms of ambiguity within the genotype data are recognized. A `+` indicates that there are multiple copies of a single gene in the individual's genotype; a `/` indicates allelic ambiguity, different allelic possibilities for the gene copy, and `|` indicates genotypic ambiguity, different genotypic possibilities for the complete genotype.
#' 
#' @param dat A vector. The genotyped data of the individuals, stating the whole genotype of a individual in one argument.
#' @param NEGs An optional character. How the `NEG` allele should be called, if \code{NULL} the `NEG` alleles will just be called `NEG`.
#' @param POS.rep A logical scalar. Whether or not all `POS` diplotype should be replaced by all other observed diplotypes. Individuals with \code{NA} as genotype will then also obtain all observed diplotypes.
#' 
#' @return A list with each element stating all unique diplotypes of a individual that are compatible with the individual's genotype.
#' 
#' @examples
#' dat = c("001/002+001/002", "003+004+005")
#' \dontrun{
#' all_options(dat = dat)
#' }
#' 
#' @export all_options
all_options = function(dat, NEGs = NULL, POS.rep = TRUE){
  len_dat = length(dat)
  geno.list = as.list(rep(NA, len_dat))
  
  if(length(NEGs) == 0 | is.null(NEGs)){
    NEGs = "NEG"
  }
  if(length(NEGs) != 1){
    stop("There is more than one NEG...")
  }

  pos_index = NULL
  for(i in 1:len_dat){
    xx = dat[i]
    
    if(is.na(xx)){  ##
      geno.list[[i]] = "POS"
      next
    }
    
    if(xx == "POS"){  
      geno.list[[i]] = "POS"
      pos_index = c(pos_index, i)
      next
    }
    
    if(xx %in% NEGs){  # xx == "NEG"){
      geno.list[[i]] = paste(NEGs, NEGs, sep = "+")
      next
    }
    
    Plus_presence = str_detect(xx, "\\+")
    Bar_presence = str_detect(xx, "\\|")
    Slash_presence = str_detect(xx, "\\/")
    
    
## Genotype has no +, | or /
    if(Plus_presence == FALSE){
      if(Bar_presence == FALSE){
        if(Slash_presence == FALSE){
          geno.list[[i]] = unique(paste(c(xx, NEGs), collapse = "+"))
          next
          
## Genotype has no + or |, but a / 
        } else {
          geno.list[[i]] = unique(multiple_single_options(xx, sign = "/", NEGs = NEGs))
          next
        }
      }
    }
    
## Genotype has no + or /, but a |
    if(Plus_presence == FALSE){
      if(Slash_presence == FALSE){
        geno.list[[i]] = unique(multiple_single_options(xx, sign = "|", NEGs = NEGs))
        next
        
## Genotype has no +, but a | and a /        
      } else {
        splt = unlist(strsplit(xx, "\\|"))
        
        mult_out = NULL
        for(j in 1:length(splt)){
          mult_out = c(mult_out, multiple_single_options(splt[j], sign = "/", NEGs = NEGs))
        }
        
        geno.list[[i]] = unique(mult_out)
        next
      }
    }
    
## Genotype has no / or |, but a +
    if(Plus_presence == TRUE){
      if(Slash_presence == FALSE){
        if(Bar_presence == FALSE){
          geno.list[[i]] = unique(multiple_gene_copies(xx, NEGs = NEGs))
          next
          
## Genotype has no /, but a + and a |            
        } else {
          splt = unlist(strsplit(xx, "\\|"))
          len_splt = length(splt)
          
          mult_out = NULL
          for(j in 1:len_splt){
            mult_out = c(mult_out, multiple_gene_copies(splt[j], NEGs = NEGs))
          }
        }
        
        geno.list[[i]] = unique(mult_out)
        next
        
## Genotype has no |, but a + and a /    
      } else {
        if(Bar_presence == FALSE){
          splt = splitting_slashes(xx)
          len_splt = length(splt)  
        
          mult_out = NULL
          for(j in 1:len_splt){
            mult_out = c(mult_out, multiple_gene_copies(xx = splt[j], NEGs = NEGs))
          }
          
          geno.list[[i]] = unique(mult_out)
          next
          
## Genotype has a +, a / and a |        
        } else {
          splt = unlist(strsplit(xx, "\\|"))
          len_splt = length(splt)
          
          mult_out = NULL
          for(j in 1:len_splt){
            splttd = splitting_slashes(splt[j])
            len_splttd = length(splttd)
            
            for(k in 1:len_splttd){
              mult_out = c(mult_out, multiple_gene_copies(splttd[k], NEGs = NEGs))
            }
          }
          
          geno.list[[i]] = unique(mult_out)
          next
        }
      }
    }
    
    
    stop("For individual ", i, " there is some unforseen problem with the all_options function")
  }
  
  
  if(POS.rep & !is.null(pos_index)){
    all_genos = unique(unlist(geno.list)); all_genos = all_genos[!("POS" == all_genos)]
    
    geno.list[pos_index] = lapply(geno.list[pos_index], function(x){x = all_genos})
  }
  
  
  return(geno.list)
}
