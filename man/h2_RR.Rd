\name{h2_RR}
\alias{h2_RR}
\title{
  Ridge Regression Based Variance Estimation in High-Dimensional Linear Models
}
\description{
  Use the ridge regression based method of Liu et al. (2020) to estimate the residual variance in (multivariate) linear model. The proposal is valid under both low- and high-dimensional models, and performs well not only for the non-sparse cases but also for the sparse ones. 
}

\usage{h2_RR(y,x,eta=NULL,alpha=0.1)
}

\arguments{
  \item{y}{A \eqn{n\times q} response matrix. \eqn{y} is an \eqn{n}-vector if \eqn{q=1}.}
  
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{eta}{The tunning parameter using in ridge regression. Default is NULL, and \code{eta} will be selected by the method proposed in Liu et al. (2020) in this case.}

  \item{alpha}{A constant to justify the tunning parameter when eta=NULL. Default is 0.1.}
  
}

\value{
  \item{sigma2}{The estimation of the residual variance, which is a \eqn{q} vector.}
  \item{h2}{The estimation of the heritability.}
}


\references{Liu, X., Zheng, S. and Feng, X. (2020). Estimation of error variance via ridge regression. Biometrika, 107, 481-488.
}

\examples{ 
# example 1
n   <- 80
p   <- 100
eta <- 1e-6
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
x <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)
y <- y + x\%*\%beta
fit <- h2_RR(y,x,eta)

# example 2
n   <- 80
p   <- 100
q   <- 5
eta <- 1e-6
rho <- 0.8
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
sig = chol( outer(1:q,1:q,function(i,j) rho^abs(i-j)) )
y <- matrix(rnorm(n*q),n,q) \%*\% sig
x <- matrix(rnorm(n*p),n,p)
y[,1] <- y[,1] + x \%*\% beta
fit <- h2_RR(y,x,eta)
}