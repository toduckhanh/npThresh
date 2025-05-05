#' @importFrom Rcpp evalCpp
#' @useDynLib npThresh, .registration = TRUE

#' @export
cdf_kernel <- function(x, X, kernel_type, bwd){
  if(any(is.na(x)) | any(is.na(X))) return(NA)
  else return(cdf_kernel_C(x, X, kernel_type, bwd))
}
