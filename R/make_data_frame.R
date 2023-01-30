#' @title Makes a data frame from a list
#'
#' @description Support function for \code{\link{splitting_slashes}}.
#'
#' @param lst A list. The genotype options splitted (\emph{e.g.} via \code{\link[base]{strsplit}}), on the allelic ambiguities (`/`).
#' 
#' @return A data_frame with the options of each gene copy as a different row. The number of columns are equal to the gene copy with the most options, other rows are supplemented with \code{NA}.
#'
# #' @examples 
# #' dat = "001/002+003/004/005"
# #' dat_splt = unlist(strsplit(dat, "\\+")); len_dat_splt = length(dat_splt)
# #' 
# #' lst = as.list(NULL)
# #' for(j in 1:len_dat_splt){
# #'   lst[[j]] = sort(unlist(strsplit(dat_splt[j], "\\/")))
# #' }
# #' \dontrun{
# #' make_data_frame(lst)
# #' }
#' 
make_data_frame = function(lst){
  
  ncols = max(lengths(lst))
  nrows = length(lst)
  
  dat = matrix(NA, nrow = nrows, ncol = ncols)
  
  for(i in 1:nrows){
    dat[i, ] = t(lst[[i]][1:ncols])
  }
  
  return(as.data.frame(dat))
}
