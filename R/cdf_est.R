#' @importFrom Rcpp evalCpp
#' @useDynLib npThresh, .registration = TRUE

#' @export
cdf_kernel <- function(x, X, kernel_type = c("gauss", "epan"), bwd){
  if(any(is.na(x)) | any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  res <- cdf_kernel_C(x, X, kernel_type_value, bwd)
  return(res)
}
