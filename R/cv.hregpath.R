cv.hregpath <- function(outlist, lambda, x, y, foldid, delta) {    
    ###Turn y into c(0,1)
    y <- as.double(y)
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE])
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- huberreg(y-predmat, delta)
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = "Huber Loss")
} 
