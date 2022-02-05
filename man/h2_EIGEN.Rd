\name{h2_EIGEN}
\alias{h2_EIGEN}
\title{
  EigenPrism method in High-Dimensional Linear Models.
}
\description{
  Use EigenPrism method of Janson et al. (2017) to estimate the residual variance in high-dimensional linear model.
}

\usage{h2_EIGEN(y,X,invsqrtSig=NULL,zero.ind=c())
}

\arguments{
  \item{y}{A length \eqn{n} vector of response.}  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  \item{invsqrtSig}{If columns of X not independent, p by p positive definite matrix which is the square-root of the inverse of Sig, where Sig is the *correlation* matrix of the X, default is identity.}
  \item{zero.ind}{Vector of which indices of the weight vector w to constrain to zero, default is none.}
}

\value{
  \item{sigma2}{The estimation of the residual variance.}
  \item{h2}{The estimation of the heritability.}
}

\references{
  JANSON, L., BARBER, R. F. & CANDES, E. (2017). EigenPrism: inference for high dimensional signal-to-noise ratios.
}

\examples{ 
n   <- 80
p   <- 100
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
x <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)
y <- y + x \%*\% beta
fit <-h2_EIGEN(y,x)
}