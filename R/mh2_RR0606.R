
mh2_RR0606 <-function(y,x,l,eta=NULL,alpha=0.1){
  if(is.null(x)) stop("x must not be NA")
  if(is.null(y)) stop("y must not be NA")
  n = nrow(x)
  p = ncol(x)
  q = ncol(y)
  
  if(is.null(p)) p = 1
  if(is.null(q)) q = 1
  is_eta = ifelse(is.null(eta),0,1)
  
  is.debias = F
  para = c(n,p,q,is_eta)

  fit <- .Call("mh2_RR_test", t(x), y, l, eta, as.integer(para), as.integer(is.debias), alpha)
  fit$sigma2 <- matrix(fit$sigma2,q)
  fit
}
