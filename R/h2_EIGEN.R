library(cccp)
h2_EIGEN <- function(y,X,invsqrtSig=NULL,zero.ind=c()){
  # Author: Lucas Janson (statweb.stanford.edu/~ljanson)
  # Runs EigenPrism procedure for estimating and generating confidence
  #  intervals for variance components in high-dimensional linear model:
  #       y = X%*%beta + e,   rows of X iid~ N(0,Sig),   e iid~ N(0,sigma^2)
  #  Requires cccp package for solving second order cone optimization.
  #  Note confidence interval endpoints may lie outside parameter domain, so it may be appropriate
  #   to clip them after the fact.
  # 
  # Inputs:
  #  y: response vector of length n (will automatically be centered)
  #  X: n by p design matrix; columns will automatically be centered and scaled to variance 1;
  #      should not contain intercept column, since both y and X will be centered
  #  invsqrtSig: if columns of X not independent, p by p positive definite matrix which is the square-root
  #               of the inverse of Sig, where Sig is the *correlation* matrix of the X (default is identity)
  #  alpha: significance level for confidence interval (default = 0.05)
  #  target: target of estimation/inference
  #		  'beta2' (default) is the squared 2-norm of the coefficient vector: sum(beta^2)
  #           'sigma2' is the noise variance sigma^2
  #           'heritability' is the fraction of variance of y explained by X%*%beta: t(beta)%*%Sig%*%beta/var(y)
  #  zero.ind: vector of which indices of the weight vector w to constrain to zero (default is none)
  #  diagnostics: boolean (default = T) for whether to generate diagnostic plots for the V_i, lambda_i, and w_i
  #  
  # Outputs:
  #  estimate: unbiased estimate of the target (for heritability, only approximately unbiased)
  #  CI: 100*(1-alpha)% confidence interval for target
  
  # Get dimensionality of problem
  n = nrow(X)
  p = ncol(X)
  
  # Transform y and X to proper form
  y = y-mean(y)
  X = scale(X,T,T)*n/(n-1)
  if(!is.null(invsqrtSig)) X = X%*%invsqrtSig
  
  # Take singular value decomposition and rescale singular values
  svd = svd(X)
  lambda = svd$d^2/p
  
  # Defined cone-constrained linear problem to optimize weights; [v; w] is vector of optimization variables
  q = c(1,rep(0,n)) #coefficient vector in objective function
  A = rbind(c(0,rep(1,n)),c(0,lambda)) #matrix for linear constraints
  b = c(1,0) #vector for linear constraints
  # Constrain some weights to be zero if desired
  if(!is.null(zero.ind)){
    A = rbind(A,cbind(rep(0,length(zero.ind)),diag(rep(1,n))[zero.ind,]))
    b = c(b,rep(0,length(zero.ind)))
  }
  # Define second-order cone constraints
  soc1 = socc(diag(c(1/4,rep(1,n))),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
  soc2 = socc(diag(c(1/4,lambda)),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
  prob = dlp(as.vector(q),A,as.vector(b),list(soc1,soc2))
  
  # Solve optimization problem and extract variables
  opt = cps(prob,ctrl(trace=F))
  v = getx(opt)[1]
  w = getx(opt)[-1]
  
  # Compute estimate and y's variance
  est = sum(w*(t(svd$u)%*%y)^2)
  yvar = sum(y^2)/n
  h2 = 1- est/yvar
  # Generate list with results
  result=list()
  result$sigma2 = est
  result$h2 = h2
  
  return(result)
}
