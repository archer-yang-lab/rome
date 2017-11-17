\name{cv2d.krome}
\alias{cv2d.krome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{2D Cross-validation for krome}
\description{
Does 2D k-fold cross-validation for \code{krome}, and returns the best value for \code{lambda} and \code{sigma} (for Gaussian and Laplace kernels) or \code{scale} for Hyperbolic tanget kernel.
}
\usage{
\method{cv2d}{krome}(x, y, kname = c("rbfdot","laplacedot","tanhdot"), lambda = NULL, sigma = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{y}{response variable.
}
\item{kname}{the name of the kernel to be used. Currently  \code{cv2d.krome} only supports 
three kernels:
\itemize{
\item\code{"rbfdot"} Radial Basis kernel function,
\item\code{"laplacedot"} Laplacian kernel function,
\item\code{"tanhdot"} Hyperbolic tangent kernel function,
}
}
\item{lambda}{a user supplied \code{lambda} sequence. It is better to supply a decreasing sequence of \code{lambda} values, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.
}
\item{sigma}{a user supplied \code{sigma} sequence, as the parameters for the kernels.
}
  \item{\dots}{other arguments that can be passed to \code{krome}.}
}
\details{
The algorithm generates a 2D grid of (\code{lambda}, \code{sigma}) pairs from the user specified \code{lambda} and \code{sigma} sequences. Then for each pair, the function computes the corresponding averaged cross-validation errors. The values and locations of the optimal \code{lambda} and \code{sigma} that give the smallest error are return, as well as the whole 2D grid of cross-validation errors.
}
\value{
an object of class \code{\link{cv2d.krome}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{out_mat}{2D grid of cross validation errors for (\code{lambda}, \code{sigma}) pairs}
		\item{mm.cvm}{the optimal value of cross validation error}
		\item{loc.lambda}{the location of the optimal lambda}
		\item{loc.sigma}{the location of the optimal sigma}
		\item{mm.lambda}{the optimal lambda value}
		\item{mm.sigma}{the optimal sigma value}
}


\references{
Y. Yang, T. Zhang, and H. Zou. (2017) "Flexible Expectile Regression in Reproducing Kernel Hilbert Space." Technometrics. Accepted.
}
\examples{
# generate data
N <- 200
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- 3*runif(N)
SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + X1*X3
sigma <- sqrt(var(Y)/SNR)
Y <- Y + X2*rnorm(N,0,sigma)
X <- cbind(X1,X2,X3)

# define lambda and sigma sequence
lambda <- exp(seq(log(1),log(0.1),len=10))
sigma <- exp(seq(log(0.01),log(0.001),len=10))

# perform 2D CV
cv1 <- cv2d.krome(X, Y, 
		"rbfdot", 					     
		lambda = lambda, sigma = sigma, 
		eps = 1e-6, maxit = 500, 
		delta = 2, gamma = 1e-06)
cat("\n loc.lambda \n ", cv1$loc.lambda)
cat("\n loc.sigma \n ",cv1$loc.sigma)

# use optimal (lambda,sigma) to refit
kern <- rbfdot(sigma=cv1$mm.sigma)
m1 <- krome(X, Y, 
kern, 					     
lambda = cv1$mm.lambda, 
eps = 1e-6, maxit = 500, 
delta = 0.5, gamma = 1e-06)	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
