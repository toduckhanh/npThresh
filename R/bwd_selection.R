# plugin method ----
#' @export
bwd_plugin <- function(X, bwd_method = c("RT", "AZZ", "PB"), hc = NULL){
  bwd.method <- match.arg(bwd_method)
  if(bwd_method %in% c("RT", "AZZ") && is.null(hc)) {
    hc <- switch (bwd_method,
      RT = 0.9,
      AZZ = 1.3
    )
  }
  bwd <- switch(bwd_method,
                RT = hc*min(c(sd(X), IQR(X)/1.34))*length(X)^(-1/5),
                AZZ = hc*sd(X)*length(X)^(-1/3)
  )
  return(bwd)
}

# leave-one-out cross validation method ----
#' @export
bwd_cv <- function(X, kernel_type = c("gauss", "epan"),
                   n_bwd = 101, n_X = 151){
  if(any(is.na(x)) | any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  x_grid <- seq(min(x), max(x), length.out = n_X)
  hx <- diff(x_grid)[1]
  res <- cv_bwd_C(x = x_grid, X = X, Ktype = kernel_type_value, hx = hx,
                  bwd_pts = n_bwd)
  return(res)
}
