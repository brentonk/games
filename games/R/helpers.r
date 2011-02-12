##' @include games.r
NULL

##
## INPUT:
## boot: number of bootstrap iterations
## report: whether to print status bar
## estimate: original fit coefficients
## y: dependent variable
## a: acceptance vector (for ultimatum only)
## regr: list of regressor matrices
## fn: log-likelihood function
## gr: gradient function (if any)
## fixed: logical vector indicating which parameters are held fixed
## method: optimization routine to use
##
## RETURN:
## matrix of bootstrap results, each row an iteration
## 
gameBoot <- function(boot, report = TRUE, estimate, y, a = NULL, regr, fn, gr,
                     fixed, method, ...)
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
        bootResults <- maxLik(logLik = fn, grad = gr, start = estimate, fixed =
                              fixed, method = method, y = newy, acc = newa, regr
                              = newregr, ...)
        cc <- convergenceCriterion(method)
        if (!(bootResults$code %in% cc)) {
            warning("bootstrap iteration ", i,
                    " failed to converge and will be removed")
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

## TODO: get rid of this function
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
## INPUT:
## x: numeric vector of values in [0, 1]
##
## RETURN:
## numeric vector, ensuring all values of x are numerically inside (0, 1)
##
finiteProbs <- function(x)
{
    x <- replace(x, x < .Machine$double.eps, .Machine$double.eps)
    x <- replace(x, x > 1 - .Machine$double.neg.eps,
                 1 - .Machine$double.neg.eps)
    return(x)
}

##
## INPUT:
## x: numeric vector
##
## RETURN:
## numeric vector, replacing Inf with largest representable values
##
finitize <- function(x)
{
    x <- ifelse(is.finite(x), x, sign(x) * .Machine$double.xmax)
    return(x)
}

##
## INPUT:
## f: object inheriting from class "formula", or a list of such objects
## argname: character string specifying the name of the argument being checked
## in the original function (in order to give an informative error message in
## case of failure)
##
## RETURN:
## object of class "Formula", combining supplied formulas (if 'f' is a list)
## into a big one with multiple right-hand sides
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
## INPUT:
## ...: vectors of any standard class
##
## RETURN:
## vector of elements contained in all vectors in '...'
##
intersectAll <- function(...)
{
    x <- list(...)
    ans <- x[[1]]
    for (i in 1:length(x)) ans <- intersect(ans, x[[i]])

    return(ans)
}

##
## INPUT:
## varNames: list of character vectors, each containing the variable names for
## one utility or variance equation
## prefixes: character vector containing names of utility equations (but not
## variance terms)
## link: "logit" or "probit"
## sdterms: number of variance equations
##
## RETURN:
## varNames: character vector of variable names
## hasColon: logical vector indicating which utility/variance equations are not
## fixed to 0 or contain only a constant
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
## INPUT:
## b: parameter vector
## regr: list of regressor matrices
## nutils: number of utility equations
## unames: names of utility equations
## finit: whether to coerce utilities not to contain any Infs
##
## RETURN:
## list of numeric vectors (named according to 'unames') of fitted utilities,
## along with element 'b' containing unused parameters (those pertaining to
## variance terms)
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
## INPUT:
## x: output from profile.game
##
## RETURN:
## vector of parameters giving the highest profiled log-likelihood
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
