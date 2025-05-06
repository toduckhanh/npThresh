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
#' @importFrom randtoolbox sobol
#' @export
bwd_cv <- function(X, kernel_type = c("gauss", "epan"),
                   n_bwd = 101, n_X = 151, method_grid = c("grid", "sobol")){
  if(any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  x_grid <- seq(min(x), max(x), length.out = n_X)
  hx <- diff(x_grid)[1]
  range_X <- max(X) - min(X)
  if(method_grid == "grid"){
    bwd_seq <- seq(range_X/200, range_X/2, length.out = n_bwd)
  } else{
    bwd_seq_temp <- sobol(n = n_bwd - 2, dim = 1)
    bwd_seq <- bwd_seq_temp*(range_X/2 - range_X/200) + range_X/200
    bwd_seq <- c(range_X/200, sort(bwd_seq), range_X/2)
  }
  res <- cv_bwd_C(x = x_grid, X = X, Ktype = kernel_type_value, hx = hx,
                  bwd_seq = bwd_seq)
  return(res)
}

# PCO method ----
#' @export
bwd_pco <- function(X, kernel_type = c("gauss", "epan"), n_bwd = 101,
                    h_min = NULL, method_grid = c("grid", "sobol")){
  if(any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  method_grid <- match.arg(method_grid)
  n_X <- length(X)
  range_X <- max(X) - min(X)
  if(is.null(h_min)) h_min <- 1/(sqrt(2*pi)*n_X)
  if(method_grid == "grid"){
    bwd_seq <- seq(h_min, range_X/2, length.out = n_bwd)
  } else{
    bwd_seq_temp <- sobol(n = n_bwd - 2, dim = 1)
    bwd_seq <- bwd_seq_temp*(range_X/2 - h_min) + h_min
    bwd_seq <- c(h_min, sort(bwd_seq), range_X/2)
  }
  res <- PCO_bwd_C(X = X, Ktype = kernel_type_value, bwd_seq = bwd_seq)
  return(res)
}



