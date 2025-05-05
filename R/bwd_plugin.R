## bandwidth
#' @export
bwd_plugin <- function(x, bwd_method = c("RT", "AZZ", "PB"), hc = NULL){
  bwd.method <- match.arg(bwd_method)
  if(bwd_method %in% c("RT", "AZZ") && is.null(hc)) {
    hc <- switch (bwd_method,
      RT = 0.9,
      AZZ = 1.3
    )
  }
  bwd <- switch(bwd_method,
                RT = hc*min(c(sd(x), IQR(x)/1.34))*length(x)^(-1/5),
                AZZ = hc*sd(x)*length(x)^(-1/3)
  )
  return(bwd)
}
