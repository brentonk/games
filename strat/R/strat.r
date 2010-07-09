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

assignList <- function(x)
{
    for (i in seq_len(length(x)))
        assign(names(x)[i], x[[i]], pos = parent.frame())

    invisible(x)
}

makeUtils22 <- function(b, regr)
{
    utils <- vector("list", 4)
    names(utils) <- c("u11", "u13", "u14", "u24")
    
    rcols <- sapply(regr, ncol)
    for (i in 1:4) {
        if (rcols[i] > 0) {
            utils[[i]] <- as.numeric(regr[[i]] %*% b[1:rcols[i]])
            b <- b[-(1:rcols[i])]
        } else {
            utils[[i]] <- rep(0, nrow(regr[[i]]))
        }
    }

    utils$b <- b
    return(utils)
}

makeSDs22 <- function(b, regr)
{
    sds <- vector("list", 4)

    rcols <- sapply(regr, ncol)
    for (i in 1:4) {
        if (rcols[i+4] > 0) {
            sds[[i]] <- exp(as.numeric(regr[[i+4]] %*% b[1:rcols[i+4]]))
            b <- b[-(1:rcols[i+4])]
        } else {
            sds[[i]] <- rep(1, nrow(regr[[i+4]]))
        }
    }

    return(sds)
}

makeProbs22 <- function(b, regr, link, type)
{
    utils <- makeUtils22(b, regr)

    ## length(utils$b) == 0 means no terms left for the variance components
    if (length(utils$b) == 0) {
        sds <- list(1, 1, 1, 1)
    } else {
        sds <- makeSDs22(b, regr)
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

    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

logLik22 <- function(b, y, regr, link, type)
{
    probs <- makeProbs22(b, regr, link, type)
    probs <- do.call(cbind, probs)
    logProbs <- cbind(log(probs[, 1]), log(probs[, 2]) + log(probs[, 3]),
                      log(probs[, 2]) + log(probs[, 4]))
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

logLikGrad22 <- function(b, y, regr, link, type)
{
    utils <- makeUtils22(b, regr)
    probs <- makeProbs22(b, regr, link, type)
    rcols <- sapply(regr, ncol)
    assignList(utils)
    assignList(probs)

    dp4 <- dnorm(u24)
    dgp4 <- dp4 * regr$Z
    Dp4 <- cbind(matrix(0L, nrow = nrow(regr$X1), ncol = sum(rcols[1:3])),
                   dgp4)
    Dp3 <- -Dp4

    dp1 <- dnorm(u11 - p3 * u13 - p4 * u14, sd = sqrt(1 + p3^2 + p4^2))
    dp1 <- dp1 / sqrt(1 + p3^2 + p4^2)
    dbp1 <- dp1 * cbind(regr$X1, -p3 * regr$X3, -p4 * regr$X4)
    dgp1 <- dp1 * (u13 - u14) * dgp4
    Dp1 <- cbind(dbp1, dgp1)
    Dp2 <- -Dp1

    dL1 <- Dp1 / p1
    dL3 <- Dp2 / p2 + Dp3 / p3
    dL4 <- Dp2 / p2 + Dp4 / p4

    ans <- as.numeric(y == 1) * dL1 + as.numeric(y == 2) * dL3 +
        as.numeric(y == 3) * dL4
    return(ans)
}

strat22 <- function(formulas, data, subset, na.action,
                    varformulas,
                    link = c("probit", "logit"),
                    type = c("private", "agent"),
                    startvals = c("zero", "unif", "sbi"),
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

    if (link == "logit" && type == "private") {
        type <- "agent"
        warning("logit errors cannot be used with private information model; ",
                "changing to agent model with logit link")
    }

    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    yf <- as.factor(model.response(mf))
    if (nlevels(yf) != 3) stop("dependent variable must have 3 values")
    y <- as.numeric(yf)

    regr <- list()
    for (i in seq_len(length(formulas)[2]))
        regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
    names(regr) <- c("X1", "X3", "X4", "Z", "V1", "V2", "V3",
                     "V4")[1:length(regr)]    
    rcols <- sapply(regr, ncol)

    if (startvals == "zero") {
        sval <- rep(0, sum(rcols))
    } else if (startvals == "unif") {
        sval <- runif(sum(rcols), -1, 1)
    }

    varNames <- sapply(regr, colnames)
    prefixes <- paste(c(rep("u1(", 3), "u2("), c(levels(yf), levels(yf)[3]),
                       ")", sep = "")
    prefixes <- c(prefixes, "v1", "v2", "v3", "v4")
    for (i in seq_len(length(formulas)[2])) {
        if (length(varNames[[i]]))
            varNames[[i]] <- paste(prefixes[i], varNames[[i]], sep = ":")
    }
    names(sval) <- unlist(varNames)

    gr <- if (type == "private" && link == "probit") logLikGrad22 else NULL

    results <- maxBFGS(fn = logLik22, grad = gr, start = sval, y = y, regr =
                       regr, link = link, type = type, ...)

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
