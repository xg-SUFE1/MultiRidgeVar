h2_MM <- function(y,x,identity=F,Sigma=NULL){
  if(is.null(x)) stop("x must not be NA")
  if(is.null(y)) stop("y must not be NA")
  n = nrow(x)
  p = ncol(x)
  if(is.null(p)) p = 1
  para = c(n,p)

  fit <- .Call("h2_MM_test", y,  t(x), Sigma, as.integer(para), as.integer(identity),
               as.integer(is.null(Sigma)))
  fit
}
