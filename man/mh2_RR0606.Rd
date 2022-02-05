\name{mh2_RR0606}
\alias{mh2_RR0606}
\title{
  Ridge Regression Based Covariance Estimation in High-Dimensional Multivariate Linear Models
}
\description{
  Use the ridge regression based method of Liu et al. (2020) to estimate the residual covariance in multivariate linear model. The proposal is valid under both low- and high-dimensional models, and performs well not only for the non-sparse cases but also for the sparse ones. 
}

\usage{mh2_RR0606(y,x,l,eta=NULL,alpha=0.1)
}

\arguments{
  \item{y}{A \eqn{n\times q} response matrix.}
  
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
 
  \item{l}{Combination coefficients for multiple responses.}
  
  \item{eta}{The tunning parameter using in ridge regression. Default is NULL, and \code{eta} will be selected by the method proposed in Liu et al. (2020) in this case.}

  \item{alpha}{A constant to justify the tunning parameter when eta=NULL. Default is 0.1.}
  
}

\value{
  \item{sigma2}{The estimation of the residual covariance.}
  \item{trA1}{The value of \eqn{n - \mbox{tr}(A_{1n})}, see the details in Liu, X., Zheng, S. and Feng, X. (2020).}
}


\references{Liu, X., Zheng, S. and Feng, X. (2020). Estimation of error variance via ridge regression. Biometrika, 107, 481-488.
}

\examples{ 
n <- 100
p <- 60
q <- 2
eta <- 1e-6
l <- matrix(rnorm(q*1),q,1)
l <- l/norm(l,"F")
B <- matrix(rep(0,p*q),nrow = p,ncol=q)
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
for (j in 1:q){
  B[,j] = beta
}
Sigma_E = diag(q)
E <- matrix(rnorm(n*q),n,q)
X <- matrix(rnorm(n*p),n,p)
Y <- X%*%B + E
fit <- mh2_RR0606(Y,X,l,eta)
}