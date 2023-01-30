#' Function for rounding of numbers
#'
#' @details Values that are a 1000 or higher will be rounded scientifically with the number of digits of ndigits + 1
#'
#' @param x A vector of numerical values that need to be rounded
#' @param digits An intiger. Indicating the number of decimal places
#' @param nsmall An intiger. The minimum number of digits to the right of the decimal point in formatting real/complex numbers in non-scientific formats
#'
# #' @examples
# #' \dontrun{
# #' rounding(c(2.456, 14.346, 1032.43))
# #' }
#' 
rounding = function(x, digits = 2, nsmall = 2){
  rounded = x; len_x = length(x)
  for(i in 1:len_x){
    x_i = x[i]
    
    if(floor(x_i) > 999){
      rounded[i] = format(signif(x_i, digits), scientific = TRUE)
    } else {
      rounded[i] = format(round(x_i, digits), nsmall = nsmall)
    }
  }

  return(rounded)
}
