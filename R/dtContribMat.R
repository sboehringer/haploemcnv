#' Gives the order of all haplotype combinations. Support function for \code{\link{Geno2EM_SB}}, \code{\link{EM_algorithm}} and for \code{\link{information_loss}}.
#' 
#' @param Nhts The number of unique haplotypes for the matrix 
#'
#' @return An matrix with the order of haplotype combinations
#' 
# #' @examples
# #' Nhts <- 3
# #' \dontrun{
# #' dtContribMat(Nhts)
# #' }
#'
dtContribMat = function(Nhts){
  
  ##
    Npairs = function(N)( (N^2 + N)/2 )
  ##
  
  Ndts = Npairs(Nhts);
  templ = rep(0, Nhts);
  
  ##
    pairFromGt = function(Nalleles, gt) {
      mx = Npairs2N(gt);
      return(c(gt - Npairs(mx), mx));
    }
    Npairs2N = function(Npairs)floor( -1/2 + sqrt(1/4 + 2*Npairs) )
  
    vector.assign1p = function(v, idcs, ...)vector.assign(v, idcs + 1, ...)
    vector.assign = function(v, idcs, e, na.rm = 0, N) {
      if (!missing(N)) v = rep(v, N);
      v[idcs] = e;
      if (!is.na(na.rm)) v[is.na(v)] = na.rm;
      v
    }
  ##
  
  dt2ht = apply(sapply(0:(Ndts - 1), pairFromGt, Nalleles = Nhts), 2, vector.assign1p, e = 1, v = templ);
  dt2ht2 = t(dt2ht) * 2 / apply(dt2ht, 2, sum);	# normalized to two alleles
  return(dt2ht2);
}
