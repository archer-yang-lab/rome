cv.lrome <- function(x, y, lambda = NULL, nfolds = 5, foldid, delta = 2, ...) {
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    lrome.object <- lrome(x, y, lambda = lambda, delta = delta, 
        ...)
    lambda <- lrome.object$lambda
    # predict -> coef
    nz <- sapply(coef(lrome.object, type = "nonzero"), length)
    if (missing(foldid)) 
        foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <- lrome(x = x[!which, , drop = FALSE], 
            y = y_sub, lambda = lambda, delta = delta, ...)
    }
    ###What to do depends on the model fit
    fun <- paste("cv", class(lrome.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, delta))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
        cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, lrome.fit = lrome.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.lrome"
    obj
} 
