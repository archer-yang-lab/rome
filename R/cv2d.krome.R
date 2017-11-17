cv2d.krome <-
  function(x, y, kname = c("rbfdot","laplacedot","tanhdot"), lambda = NULL, sigma = NULL, ...) {
    if (is.null(sigma)) {
      stop("user must provide a sigma sequence")
    }
    mm.cvm <- Inf
	out_mat <- matrix(NA, length(sigma), length(lambda))
    for (i in seq.int(length(sigma))) {
      if(kname=="tanhdot") kern <- do.call(kname, list(sigma[i], offset = 1))
	  else kern <- do.call(kname, list(sigma[i]))
      cv_out <- cv.krome(x, y, kern, lambda = lambda, ...)
	  out_mat[i,] <- cv_out$cvm
      if (mm.cvm > cv_out$cvm.min) {
        mm.cvm <- cv_out$cvm.min
        mm.lambda <- cv_out$lambda.min
        loc.sigma <- i
      }
      cat("sigma ", i, " completed.\n")
    }
	rownames(out_mat) <- paste("S", seq(length(sigma)), sep = "")
	colnames(out_mat) <- paste("L", seq(length(lambda)), sep = "")
    loc.lambda <- which(mm.lambda == lambda)
    list(out_mat=out_mat,
      mm.cvm = mm.cvm, loc.lambda = loc.lambda,
      loc.sigma = loc.sigma, mm.lambda = mm.lambda,
      mm.sigma = sigma[loc.sigma]
    )
  }
