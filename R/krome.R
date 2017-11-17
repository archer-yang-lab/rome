krome <- function(x, y, kern,
                 lambda = NULL, eps = 1e-08, maxit = 1e+04,
                 delta = 2, gamma = 1e-06) {
  #####################################
  #data setup
  this.call <- match.call()
  y <- drop(y)
  y <- as.double(y)
  x <- as.matrix(x)
  Kmat <- kernelMatrix(kern,x)
  diag(Kmat) <- diag(Kmat) + gamma
  np <- dim(x)
  nobs <- as.integer(np[1])
  if (length(y) != nobs)
    stop("x and y have different number of observations")
  #parameter setup
  if (delta <= 0)
    stop("delta must be positive")
  delta <- as.double(delta)
  eigen_result <- eigen(Kmat, symmetric = TRUE)
  Umat <- eigen_result$vectors
  Dvec <- eigen_result$values
  Ksum <- colSums(Kmat)
  maxit <- as.integer(maxit)
  eps <- as.double(eps)
  #lambda setup
  if (is.null(lambda)) {
    stop("user must provide a lambda sequence")
  } else {
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ################################################################################
  #call Fortran core
    fit <- .Fortran(
      "core_light", delta,
      as.double(Kmat), as.double(Umat),
      as.double(Dvec),
      nobs, as.double(y), nlam, ulam, eps, maxit, anlam = integer(1),
      npass = integer(nlam), jerr = integer(1),
      alpmat = double((nobs + 1) * nlam),
      PACKAGE = "krome"
    )

  ################################################################################
  # output
  errmsg <- err(fit$jerr, maxit)
  if (paste(errmsg$n) == '-1')
    print(errmsg$msg, call. = FALSE)
  anlam <- fit$anlam
  vnames <- paste("a", seq(nobs + 1) - 1, sep = "")
  stepnames <- paste("L", seq(anlam), sep = "")  
  alpha <-
    matrix(fit$alpmat[seq((nobs + 1) * anlam)], nobs + 1, anlam, dimnames =
             list(vnames,stepnames))
  outlist <-
    list(
      alpha = alpha, lambda = ulam[seq(anlam)], npass = fit$npass[seq(anlam)], jerr = fit$jerr
    )
  outlist$call <- this.call
  class(outlist) <- "krome"
  outlist
}
