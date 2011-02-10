##' @include games.r
NULL

##
## Calculates bootstrap results for a strategic model.
## 
gameBoot <- function(boot, report = TRUE, estimate, y, a = NULL, regr, fn, gr,
                      fixed, ...)
{
    bootMatrix <- matrix(NA, nrow = boot, ncol = length(estimate))
    failedBoot <- logical(boot)
    if (report) {
        cat("\nRunning bootstrap iterations...\n")
        pb <- txtProgressBar(min = 1, max = boot)
    }
    for (i in seq_len(boot)) {
        bootSamp <- sample(seq_len(length(y)), replace = TRUE)
        newy <- y[bootSamp]
        newa <- a[bootSamp]  ## for the ultimatum model
        newregr <- lapply(regr, function(x) x[bootSamp, , drop = FALSE])
        bootResults <- maxBFGS(fn = fn, grad = gr, start = estimate, fixed =
                               fixed, y = newy, acc = newa, regr = newregr, ...)
        if (bootResults$code) {
            warning("bootstrap iteration ", i,
                    "failed to converge and will be removed")
            failedBoot[i] <- TRUE
        }
        bootMatrix[i, ] <- bootResults$estimate

        if (report)
            setTxtProgressBar(pb, i)
    }
    if (report)
        cat("\n")
    bootMatrix <- bootMatrix[!failedBoot, , drop = FALSE]
    colnames(bootMatrix) <- names(estimate)
    return(bootMatrix)
}

##
## Calculates the variance-covariance matrix for a fitted model, including a
## procedure for catching the error (and returning a matrix of NAs) in case the
## Hessian is non-invertible.
##
getGameVcov <- function(hessian, fixed)
{
    hes <- hessian[!fixed, !fixed, drop = FALSE]
    vv <- tryCatch(solve(-hes), error = function(e) e)
    if (inherits(vv, "error")) {
        warning("variance-covariance matrix could not be calculated: ",
                vv$message)
        vv <- matrix(NA, nrow(hes), nrow(hes))
    }
    ans <- hessian
    ans[] <- NA
    ans[!fixed, !fixed] <- vv
    return(ans)
}

##
## Ensures that estimated probabilities aren't numerically equal to 1 or 0, in
## order to ensure no -Infs or 0s in log-likelihoods.
##
finiteProbs <- function(x)
{
    x <- replace(x, x < .Machine$double.eps, .Machine$double.eps)
    x <- replace(x, x > 1 - .Machine$double.neg.eps,
                 1 - .Machine$double.neg.eps)
    return(x)
}

finitize <- function(x)
{
    x <- ifelse(is.finite(x), x, sign(x) * .Machine$double.xmax)
    return(x)
}

##
## Used to ensure that the "formulas" argument of each fitting function contains
## a valid type of object and coerces it to "Formula" class.  Returns an error
## if the function isn't a formula, Formula, or list of formulas.
##
checkFormulas <- function(f, argname = "formulas")
{
    if (inherits(f, "list")) {
        f <- do.call(as.Formula, f)
    } else if (inherits(f, "formula")) {
        f <- as.Formula(f)
    } else {
        stop(argname, " must be a list of formulas or a formula")
    }

    return(f)
}

##
## Takes a list of vectors and finds their intersection
##
intersectAll <- function(...)
{
    x <- list(...)
    ans <- x[[1]]
    for (i in 1:length(x)) ans <- intersect(ans, x[[i]])

    return(ans)
}

##
## Makes the names of the variables for egame12 and egame122 models.
## 
makeVarNames <- function(varNames, prefixes, link, sdterms)
{
    vname <- if (link == "logit") "log(lambda" else "log(sigma"
    if (sdterms == 1L) {
        prefixes <- c(prefixes, paste(vname, ")", sep = ""))
    } else if (sdterms > 1L) {
        prefixes <- c(prefixes, paste(vname, 1:sdterms, ")", sep = ""))
    }

    hasColon <- sapply(varNames, function(x) length(x) > 0 &&
                       !all(x == "(Intercept)"))
    names(hasColon) <- prefixes
    for (i in seq_along(varNames)) {
        if (hasColon[i]) {
            varNames[[i]] <- paste(prefixes[i], varNames[[i]], sep = ":")
        } else {
            varNames[[i]] <- prefixes[i][length(varNames[[i]])]
        }
    }
    varNames <- unlist(varNames)
    ans <- list(varNames = varNames, hasColon = hasColon)
    return(ans)
}

##
## Calculates utilities from the given list of regressors (regr) and vector of
## coefficients (b); returns the regressors and the "remaining" coefficients
## (those for scale parameters, if any) in a list.
##
## nutils specifies which elements of regr are for utilities rather than scale
## parameters; e.g., if nutils = 4, then the first four matrices of regr are
## used to create utility vectors, and the rest are used for scale terms.
##
## unames, which must be of length nutils, specifies names for the utility terms
## in the list returned.
##
makeUtils <- function(b, regr, nutils, unames, finit = TRUE)
{
    utils <- vector("list", nutils)
    if (!missing(unames)) {
        if (length(unames) != nutils)
            stop("length(unames) must equal nutils")
        names(utils) <- unames
    }

    rcols <- sapply(regr, ncol)
    for (i in 1:nutils) {
        if (rcols[i] > 0) {
            utils[[i]] <- as.numeric(regr[[i]] %*% b[1:rcols[i]])
            if (finit)
                utils[[i]] <- finitize(utils[[i]])
            b <- b[-(1:rcols[i])]
        } else {
            utils[[i]] <- rep(0, nrow(regr[[i]]))
        }
    }

    utils$b <- b
    return(utils)
}

##
## Retrieves the "best" starting values from profile.game output
##
svalsFromProfile <- function(x)
{
    x <- do.call(rbind, x)
    bestrow <- which.max(x[, 1])  ## highest log-likelihood
    xn <- names(x)[-1]
    ans <- as.numeric(x[bestrow, ][-1])
    names(ans) <- xn
    return(ans)
}

##
## INPUT:
## method: character string describing optimization method used
##
## RETURN:
## vector of integer values corresponding to convergence codes indicating
## success (which differ by method in maxLik)
## 
convergenceCriterion <- function(method)
{
    switch(tolower(method),
           `newton-raphson` = c(1L, 2L),
           nr = c(1L, 2L),
           bfgs = 0L,
           bfgsr = c(1L, 2L),
           `bfgs-r` = c(1L, 2L),
           bhhh = c(1L, 2L),
           `nelder-mead` = 0L,
           nm = 0L,
           sann = 0L)
}
