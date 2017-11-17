\name{krome}
\alias{krome}
\title{
Fits the regularization paths for the kernel huber regression.}
\description{
Fits a regularization path for the kernel huber regression at a sequence of regularization parameters lambda.
}
\usage{
krome(x, y, kern, lambda = NULL, eps = 1e-08, maxit = 1e4,
delta = 2, gamma = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{y}{response variable.
}
\item{kern}{the built-in kernel classes in \pkg{krome}. 
The \code{kern} parameter can be set to any function, of class
kernel, which computes the inner product in feature space between two
vector arguments. \pkg{krome} provides the most popular kernel functions
which can be initialized by using the following
functions:
\itemize{
\item\code{rbfdot} Radial Basis kernel function,
\item\code{polydot} Polynomial kernel function,
\item\code{vanilladot} Linear kernel function,
\item\code{tanhdot} Hyperbolic tangent kernel function,
\item\code{laplacedot} Laplacian kernel function,
\item\code{besseldot} Bessel kernel function,
\item\code{anovadot} ANOVA RBF kernel function,
\item\code{splinedot} the Spline kernel. 
}
Objects can be created by calling the rbfdot, polydot, tanhdot, vanilladot, anovadot, besseldot, laplacedot, splinedot functions etc. (see example.)
}
\item{lambda}{a user supplied \code{lambda} sequence. It is better to supply a decreasing sequence of \code{lambda} values, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.
}
  \item{eps}{
	convergence threshold for majorization minimization algorithm. Each majorization descent loop continues until the relative change in any coefficient \eqn{||alpha(new)-\alpha(old)||_2^2/||\alpha(old)||_2^2} is less than \code{eps}. Defaults value is \code{1e-8}.
	}
  \item{maxit}{
maximum number of loop iterations allowed at fixed lambda value. Default is 1e4. If models do not converge, consider increasing \code{maxit}.
}
  \item{delta}{the parameter \eqn{\delta}{delta} in the huber regression loss. The value must be positive. Default is 2.
}
  \item{gamma}{
 a scalar number. If it is specified, the number will be added to each diagonal element of the kernel matrix as perturbation. The default is \code{1e-06}.
}
}
\details{
Note that the objective function in \code{krome} is
\deqn{Loss(y- \alpha_0 - K * \alpha ) + \lambda * \alpha^T * K * \alpha,}
where the \eqn{\alpha_0} is the intercept, \eqn{\alpha} is the solution vector, and \eqn{K} is the kernel matrix with \eqn{K_{ij}=K(x_i,x_j)}. Users can specify the kernel function to use, options include Radial Basis kernel, Polynomial kernel, Linear kernel, Hyperbolic tangent kernel, Laplacian kernel, Bessel kernel, ANOVA RBF kernel, the Spline kernel. Users can also tweak the penalty by choosing different \eqn{lambda}. 

For computing speed reason, if models are not converging or running slow, consider increasing \code{eps} before increasing \code{maxit}.
}
\value{
An object with S3 class \code{\link{krome}}.
		\item{call}{the call that produced this object.}
		\item{alpha}{a \code{nrow(x)*length(lambda)} matrix of coefficients. Each column is a solution vector corresponding to a lambda value in the lambda sequence.}
		\item{lambda}{the actual sequence of \code{lambda} values used.}
		\item{npass}{total number of loop iterations corresponding to each lambda value.}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}

\author{
Yi Yang, Teng Zhang and Hui Zou\cr
Maintainer: Yi Yang  <yi.yang6@mcgill.ca>
}

\references{
Y. Yang, T. Zhang, and H. Zou. (2017) "Flexible Expectile Regression in Reproducing Kernel Hilbert Space." Technometrics. Accepted.
}

\examples{
# create data
N <- 200
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- 3*runif(N)
SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + X1*X3
sigma <- sqrt(var(Y)/SNR)
Y <- Y + X2*rnorm(N,0,sigma)
X <- cbind(X1,X2,X3)

# set gaussian kernel 
kern <- rbfdot(sigma=0.1)

# define lambda sequence
lambda <- exp(seq(log(0.5),log(0.01),len=10))

# run krome
m1 <- krome(x=X, y=Y, kern=kern, lambda = lambda, delta = 2) 

# plot the solution paths
plot(m1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
