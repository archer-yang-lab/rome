\name{cv.lrome}
\alias{cv.lrome}
\alias{cv.hregpath}
\title{Cross-validation for lrome}
\description{Does k-fold cross-validation for lrome, produces a plot,
and returns a value for \code{lambda}. This function is modified based on the \code{cv} function from the \code{glmnet} package.}
\usage{
cv.lrome(x, y, lambda = NULL, nfolds = 5, foldid, delta=2,...)
}
\arguments{
		\item{x}{\code{x} matrix as in \code{\link{lrome}}.}
		\item{y}{response variable or class label \code{y} as in \code{\link{lrome}}.}
		\item{lambda}{optional user-supplied lambda sequence; default is
		\code{NULL}, and \code{\link{lrome}} chooses its own sequence.}
		\item{nfolds}{number of folds - default is 5. Although \code{nfolds}
		can be as large as the sample size (leave-one-out CV), it is not
		recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
		\item{foldid}{an optional vector of values between 1 and \code{nfold}
		identifying what fold each observation is in. If supplied,
		\code{nfold} can be missing.}
		\item{delta}{parameter \eqn{\delta}{delta} used in lasso huber regression for computing huber loss function}
		\item{\dots}{other arguments that can be passed to lrome.}

}
\details{The function runs \code{\link{lrome}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the
folds are computed. 
}
\value{an object of class \code{\link{cv.lrome}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating type of measure (for plotting
		purposes).}
		\item{lrome.fit}{a fitted \code{\link{lrome}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is
		within 1 standard error of the minimum.}
}


\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{https://github.com/emeryyi/fastcox.git}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}


\seealso{\code{\link{lrome}}, \code{\link{plot.cv.lrome}}, \code{\link{predict.cv.lrome}}, and \code{\link{coef.cv.lrome}} methods.}
\examples{
# fit an elastic net penalized HHSVM 
# with lambda2 = 0.1 for the L2 penalty. Use the 
# misclassification rate as the cross validation 
# prediction loss. Use five-fold CV to choose 
# the optimal lambda for the L1 penalty.

data(FHT)
set.seed(2011)
cv=cv.lrome(FHT$x, FHT$y, lambda2=0.1, 
	nfolds=5, delta=1.5)
plot(cv)

# fit an elastic net penalized least squares 
# with lambda2 = 0.1 for the L2 penalty. Use the 
# least square loss as the cross validation 
# prediction loss. Use five-fold CV to choose 
# the optimal lambda for the L1 penalty.

set.seed(2011)
cv1=cv.lrome(FHT$x, FHT$y_reg,  
lambda2=0.1, nfolds=5)
plot(cv1)

# To fit a LASSO penalized logistic regression
# we set lambda2 = 0 to disable the L2 penalty. Use the 
# logistic loss as the cross validation 
# prediction loss. Use five-fold CV to choose 
# the optimal lambda for the L1 penalty.

set.seed(2011)
cv2=cv.lrome(FHT$x, FHT$y, lambda2 = 0, nfolds=5)
plot(cv2)
}
\keyword{models}
\keyword{regression}

