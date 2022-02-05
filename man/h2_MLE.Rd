\name{h2_MLE}
\alias{h2_MLE}
\title{
  Likelihood Based Variance Estimation in High-Dimensional Linear Models
}
\description{
  Use the maximum likelihood based method of Dicker and Erdogdu (2016) to estimate the residual variance in high-dimensional linear model.
}

\usage{h2_MLE(y,x,max.iter=100,tol=1e-4)
}

\arguments{
  \item{y}{A length \eqn{n} vector of response.}
  
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{max.iter}{Maximum number of iterations. Default is 50.}
  
  \item{tol}{Convergence threshhold. Default is 1e-4.}
}

\value{
  \item{sigma2}{The estimation of the residual variance.}
  \item{h2}{The estimation of the heritability.}
}


\references{Dicker, L. H. and Erdogdu, M. A. (2016). Maximum likelihood for variance estimation in high-dimensional linear models. In Proceedings of the 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016), 159-167. JMLR Workshop \eqn{\&} Conference Proceedings.
}

\examples{ 
n   <- 80
p   <- 100
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
x <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)
y <- y + x \%*\% beta
fit <- h2_MLE(y,x)
}