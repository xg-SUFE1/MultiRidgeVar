h2_MLE <- function(y,x,max.iter=100,tol=1e-4){
  if(is.null(x)) stop("x must not be NA")
  if(is.null(y)) stop("y must not be NA")
  n = nrow(x)
  p = ncol(x)
  if(is.null(p)) p = 1
  para = c(n,p)

  fit <- .Call("h2_MLE_test", y, t(x), as.integer(para), as.integer(max.iter), as.double(tol))
  fit
}
