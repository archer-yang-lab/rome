\name{plot.cv.lrome}
\alias{plot.cv.lrome}
\title{plot the cross-validation curve produced by cv.lrome}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv} function from the \code{glmnet} package.}
\usage{
\method{plot}{cv.lrome}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.lrome}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{A plot is produced.}



\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{https://github.com/emeryyi/fastcox.git}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.lrome}}.}
\examples{
# fit an elastic net penalized logistic regression 
# with lambda2 = 1 for the L2 penalty. Use the 
# logistic loss as the cross validation 
# prediction loss. Use five-fold CV to choose 
# the optimal lambda for the L1 penalty.
data(FHT)
set.seed(2011)
cv=cv.lrome(FHT$x, FHT$y, lambda2 = 1, nfolds=5)
plot(cv)
}
\keyword{models}
\keyword{regression}
