################################################################################
### strat.r
################################################################################

library(Formula)
library(maxLik)

finiteProbs <- function(x)
{
    x <- replace(x, x < .Machine$double.eps, .Machine$double.eps)
    x <- replace(x, x > 1 - .Machine$double.neg.eps,
                 1 - .Machine$double.neg.eps)
    return(x)
}

makeUtils22 <- function(b, formulas, mf)
{
    utils <- vector("list", 4)
    names(utils) <- c("u11", "u13", "u14", "u24")

    for (i in 1:4) {
        X <- model.matrix(formulas, data = mf, rhs = i)
        if (ncol(X) > 0) {
            utils[[i]] <- as.numeric(X %*% b[1:ncol(X)])
            b <- b[-(1:ncol(X))]
        } else {
            utils[[i]] <- rep(0, nrow(X))
        }
    }

    utils$b <- b
    return(utils)
}

makeSDs22 <- function(b, formulas, mf)
{
    sds <- vector("list", 4)

    for (i in 1:4) {
        X <- model.matrix(formulas, data = mf, rhs = 4 + i)
        if (ncol(X) > 0) {
            sds[[i]] <- exp(as.numeric(X %*% b[1:ncol(X)]))
            b <- b[-(1:ncol(X))]
        } else {
            sds[[i]] <- rep(1, nrow(X))
        }
    }

    return(sds)
}

makeProbs22 <- function(b, formulas, mf, link, type)
{
    utils <- makeUtils22(b, formulas, mf)

    ## length(utils$b) == 0 means no terms left for the variance components
    if (length(utils$b) == 0) {
        sds <- list(1, 1, 1, 1)
    } else {
        sds <- makeSDs22(b, formulas, mf)
    }

    linkfcn <- switch(link,
                      logit = function(x, ...) { 1 / (1 + exp(-x)) },
                      probit = pnorm)

    sd4 <- if (type == "private") sds[[4]] else sqrt(sds[[3]]^2 + sds[[4]]^2)
    p4 <- finiteProbs(linkfcn(utils$u24, sd = sd4))
    p3 <- 1 - p4

    if (type == "private") {
        sd2 <- sqrt(p3^2 * sds[[2]]^2 + p4^2 * sds[[3]]^2 + sds[[1]]^2)
    } else {
        sd2 <- sqrt(sds[[1]]^2 + sds[[2]]^2)
    }
    p2 <- p3 * utils$u13 + p4 * utils$u14 - utils$u11
    p2 <- finiteProbs(linkfcn(p2, sd = sd2))
    p1 <- 1 - p2

    return(cbind(p1, p2, p3, p4))
}

logLik22 <- function(b, y, formulas, mf, link, type)
{
    probs <- makeProbs22(b, formulas, mf, link, type)
    logProbs <- cbind(log(probs[, 1]), log(probs[, 2]) + log(probs[, 3]),
                      log(probs[, 2]) + log(probs[, 4]))
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

strat22 <- function(formulas, data, subset, na.action,
                    varformulas,
                    link = c("probit", "logit"),
                    type = c("private", "agent"),
                    startvals = c("sbi", "zero"),
                    ...)
{
    cl <- match.call()

    link <- match.arg(link)
    type <- match.arg(type)
    startvals <- match.arg(startvals)

    if (inherits(formulas, "list")) {
        formulas <- do.call(as.Formula, formulas)
    } else if (inherits(formulas, "formula")) {
        formulas <- as.Formula(formulas)
    } else {
        stop("formulas must be a list of formulas or a formula")
    }

    if (!missing(varformulas))
    {
        if (link == "logit") {
            link <- "probit"
            warning("logit link can't be used with regressors in variance; ",
                    "using probit link instead")
        }

        if (inherits(varformulas, "list")) {
            varformulas <- do.call(as.Formula, varformulas)
        } else if (!inherits(varformulas, "formula")) {
            stop("varformulas must be a list of formulas or a formula")
        }
        formulas <- as.Formula(formula(formulas), formula(varformulas))
    }

    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    yf <- as.factor(model.response(mf))
    if (nlevels(yf) != 3) stop("dependent variable must have 3 values")
    y <- as.numeric(yf)

    if (startvals == "zero") {
        numVars <- 0
        for (i in seq_len(length(formulas)[2])) {
            numVars <- numVars +
                ncol(model.matrix(formulas, data = mf, rhs = i))
        }
        sval <- rep(0, numVars)
    }

    varNames <- sapply(1:(length(formulas)[2]), function(x)
                       colnames(model.matrix(formulas, mf, rhs = x)))
    prefixes <- paste(c(rep("u1(", 3), "u2("), c(levels(yf), levels(yf)[3]),
                       ")", sep = "")
    prefixes <- c(prefixes, "v1", "v2", "v3", "v4")
    for (i in seq_len(length(formulas)[2])) {
        if (length(varNames[[i]]))
            varNames[[i]] <- paste(prefixes[i], varNames[[i]], sep = ":")
    }
    names(sval) <- unlist(varNames)

    results <- maxBFGS(fn = logLik22, start = sval, y = y, formulas = formulas,
                       mf = mf, link = link, type = type, ...)

    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- solve(-results$hessian)
    ans$log.likelihood <- results$maximum
    ans$call <- cl
    ans$formulas <- formulas
    ans$model <- mf
    
    return(ans)
}

summaryStrat <- function(object, ...)
{
    cf <- object$coefficients
    se <- sqrt(diag(object$vcov))
    zval <- cf / se
    pval <- 2 * pnorm(-abs(zval))
    ans <- cbind(cf, se, zval, pval)
    colnames(ans) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    printCoefmat(ans)
    invisible(ans)
}
