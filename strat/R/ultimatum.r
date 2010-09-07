##' @include strat.r
NULL

predict.ultimatum <- function(object, newdata, ...)
{
    if (missing(newdata))
        newdata <- object$model

    mf <- match(c("subset", "na.action"), names(object$call), 0L)
    mf <- object$call[c(1L, mf)]
    mf$formula <- object$formulas
    mf$data <- newdata
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## This is to prevent goofy things from happening with factor variables --
    ## it ensures that the levels of factor variables in "newdata" are the same
    ## as those in the model frame used in the original fitting
    for (i in 1:ncol(mf)) {
        if (is.factor(mf[, i])) {
            iname <- names(mf)[i]
            if (iname %in% names(object$model))
                levels(mf[, i]) <- levels(object$model[, iname])
        }
    }

    X <- model.matrix(object$formulas, data = mf, rhs = 1)
    Z <- model.matrix(object$formulas, data = mf, rhs = 2)

    b <- object$coefficients
    s1 <- exp(b[length(b) - 1])
    s2 <- exp(b[length(b)])
    b <- head(b, length(b) - 2)
    g <- tail(b, ncol(Z))
    b <- head(b, length(b) - ncol(Z))
    fit1 <- as.numeric(X %*% b)
    fit2 <- as.numeric(Z %*% g)
    Q <- object$maxOffer

    Ey <- Q - fit1 - s2 * (1 + LW(exp((Q - fit1 - s2 - fit2) / s2)))
    PrA <- 1 / (1 + exp(-(Ey - fit2) / s2))

    ans <- as.data.frame(cbind(Ey, PrA))
    names(ans) <- c("E(offer)", "Pr(accept)")
    return(ans)
}

offerCDF <- function(y, maxOffer, fit1, fit2, s1, s2)
{
    1 / (1 + exp((maxOffer - y - fit1 - s2 * (1 + exp((y - fit2) / s2))) / s1))
}

offerPDF <- function(y, maxOffer, fit1, fit2, s1, s2)
{
    num <- exp(-(maxOffer - y - fit1 - s2 * (1 + exp((y - fit2) / s2))) / s1)
    denom <- s1 * (1 + num)^2
    ans <- (num / denom) * (1 + exp((y - fit2) / s2))
    return(ans)
}

logLikUlt <- function(b, y, acc, regr, maxOffer, offerOnly, offertol = offertol,
                      ...)
{
    s1 <- exp(b[length(b) - 1])
    s2 <- exp(b[length(b)])
    b <- head(b, length(b) - 2)
    g <- tail(b, ncol(regr$Z))
    b <- head(b, length(b) - ncol(regr$Z))

    fit1 <- as.numeric(regr$X %*% b)
    fit2 <- as.numeric(regr$Z %*% g)

    prAccept <- finiteProbs(1 / (1 + exp(-(y - fit2) / s2)))

    lowball <- finiteProbs(offerCDF(0, maxOffer, fit1, fit2, s1, s2))
    highball <- finiteProbs(1 - offerCDF(maxOffer, maxOffer, fit1, fit2, s1,
                                         s2))
    interior <- offerPDF(y, maxOffer, fit1, fit2, s1, s2)

    isMax <- abs(y - maxOffer) < offertol
    isMin <- abs(y) < offertol

    ans <- ifelse(isMax, highball, ifelse(isMin, lowball, interior))
    ans <- replace(ans, ans < .Machine$double.eps, .Machine$double.eps)
    if (!offerOnly)
        ans <- ans * finiteProbs(ifelse(acc == 1, prAccept, 1 - prAccept))
    ans <- log(ans)
    return(ans)
}

logLikGradUlt <- function(b, y, acc, regr, maxOffer, offerOnly, ...)
{
    s1 <- exp(b[length(b) - 1])
    s2 <- exp(b[length(b)])
    b <- head(b, length(b) - 2)
    g <- tail(b, ncol(regr$Z))
    b <- head(b, length(b) - ncol(regr$Z))

    fit1 <- as.numeric(regr$X %*% b)
    fit2 <- as.numeric(regr$Z %*% g)

    ey <- exp((y - fit2) / s2)
    ey1ey <- ey / (1 + ey)
    Qy <- (-maxOffer + y + fit1 + s2 * (1 + ey)) / s1
    Qy1Qy <- exp(Qy) / (1 + exp(Qy))
    Qy2Qy <- exp(-Qy) / (1 + exp(-Qy))

    dfdb <- ((1 - 2 * Qy1Qy) / s1) * regr$X
    dFdb <- Qy2Qy * (regr$X / s1)
    d1mFdb <- dFdb - regr$X / s1

    dfdg <- ((-ey1ey / s2) + (-1 + 2 * Qy1Qy) * (ey / s1)) * regr$Z
    dFdg <- Qy2Qy * ey * (-regr$Z / s1)
    d1mFdg <- dFdg + (ey / s1) * regr$Z

    dfdlns1 <- 2 * Qy1Qy * Qy - Qy - 1
    dFdlns1 <- Qy2Qy * (-Qy)
    d1mFdlns1 <- dFdlns1 + Qy

    ## change this to lns2!
    dfds2 <- ((1 + ey) / s1) * (1 - (2 * Qy1Qy * exp(-Qy))) + (1 / (1 + ey))
    dfds2 <- dfds2 * ey * ((fit2 - y) / s2^2)
    dFds2 <- Qy2Qy * ((1 + ey) / s1) * ey * ((fit2 - y) / s2^2)
    d1mFds2 <- dFds2 - ((1 + ey) / s1) * ey * ((fit2 - y) / s2^2)

    df <- cbind(dfdb, dfdg, dfdlns1, dfds2)
    dF <- cbind(dFdb, dFdg, dFdlns1, dFds2)
    d1mF <- cbind(d1mFdb, d1mFdg, d1mFdlns1, d1mFds2)

    ans <- df
    ans[y == 0, ] <- dF[y == 0, ]
    ans[y == maxOffer, ] <- d1mF[y == maxOffer, ]

    if (!offerOnly) {
        ey2ey <- exp((fit2 - y) / s2) / (1 + exp((fit2 - y) / s2))

        dPdb <- matrix(0L, nrow = nrow(regr$X), ncol = ncol(regr$X))
        dPdg <- -ey2ey * regr$Z / s2
        dPdlns1 <- 0L
        dPdlns2 <- -ey2ey * (y - fit2) / s2

        d1mPdb <- dPdb
        d1mPdg <- regr$Z / s2 - dPdg
        d1mPdlns1 <- 0L
        d1mPdlns2 <- ((y - fit2) / s2) - dPdlns2

        dAcc <- cbind(dPdb, dPdg, dPdlns1, dPdlns2)
        dRej <- cbind(d1mPdb, d1mPdg, d1mPdlns1, d1mPdlns2)

        ans2 <- dAcc
        ans2[acc == 0, ] <- dRej[acc == 0, ]

        ans <- ans + ans2
    }

    return(ans)
}

##' <description>
##'
##' \preformatted{
##' .       ___ 1 ___
##' .      /         \
##' .     /           \
##' .    / y \in [0,Q] \
##' .   /               \
##' .   -----------------
##' .           /\  2
##' .          /  \
##' .         /    \
##' .        /      \
##' .     Q - y     R1
##' .     y         R2}
##' @title ultimatum game
##' @param formulas e
##' @param data e
##' @param subset e
##' @param na.action e
##' @param maxOffer e
##' @param offertol numeric: offers within \code{offertol} of \code{maxOffer}
##' will be considered to be at the maximum.  If \code{maxOffer} and all
##' observed offers are integer-valued, the value of \code{offertol} should not
##' matter.
##' @param s1 e
##' @param s2 e
##' @param outcome e
##' @param boot e
##' @param bootreport e
##' @param ... other arguments to pass to the fitting function (see
##' \code{\link{maxBFGS}})
##' @param reltol w
##' @return fitted model
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
ultimatum <- function(formulas, data, subset, na.action,
                      maxOffer, offertol = sqrt(.Machine$double.eps),
                      s1 = NULL, s2 = NULL,
                      outcome = c("both", "offer"),
                      boot = 0,
                      bootreport = TRUE,
                      ...,
                      reltol = 1e-12)
{
    cl <- match.call()

    outcome <- match.arg(outcome)
    offerOnly <- outcome == "offer"
    s1null <- is.null(s1)
    s2null <- is.null(s2)

    formulas <- checkFormulas(formulas)

    ## make the model frame
    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ya <- model.part(formulas, mf, lhs = 1, drop = TRUE)
    if (length(dim(ya))) {
        y <- ya[, 1]
        a <- ya[, 2]
        if (!all(unique(a) %in% c(0, 1)))
            stop("acceptance variable must be binary")
    } else {
        if (!offerOnly) {
            stop("a dependent variable for acceptance must be specified if",
                 " `outcome == \"both\"`; see `?stratult`")
        }
        y <- ya
    }
    if (any(y > maxOffer))
        stop("observed offers greater than maxOffer")

    regr <- list()
    regr$X <- model.matrix(formulas, data = mf, rhs = 1)
    regr$Z <- model.matrix(formulas, data = mf, rhs = 2)

    s1 <- if (s1null) log(sd(y)) else log(s1)
    s2 <- if (s2null) log(sd(y)) else log(s2)

    ## suppressing warnings in the logit fitting because fitted probabilities
    ## numerically equal to 0/1 seem to occur often
    aa <- if (exists("a")) a else as.numeric(y >= mean(y))
    m2 <- suppressWarnings(glm.fit(regr$Z, aa, family = binomial(link = "logit"),
                                   intercept = FALSE, offset = -y))
    m1 <- lsfit(regr$X, maxOffer - y, intercept = FALSE)
    sval <- c(m1$coefficients, m2$coefficients, s1, s2)

    firstTry <- logLikUlt(sval, y = y, acc = a, regr = regr, maxOffer =
                          maxOffer, offerOnly = offerOnly, offertol = offertol)
    if (!is.finite(sum(firstTry))) {
        sval <- c(maxOffer - mean(y), rep(0, length(m1$coefficients) - 1),
                  maxOffer - mean(y), rep(0, length(m2$coefficients) - 1),
                  s1, s2)
    }

    names(sval) <- c(paste("R1", colnames(regr$X), sep = ":"),
                     paste("R2", colnames(regr$Z), sep = ":"),
                     "log(s1)", "log(s2)")

    fvec <- logical(length(sval))
    if (!s1null) fvec[length(fvec) - 1] <- TRUE
    if (!s2null) fvec[length(fvec)] <- TRUE
    names(fvec) <- names(sval)

    results <- maxBFGS(fn = logLikUlt, grad = logLikGradUlt, start = sval, fixed
                       = fvec, y = y, acc = a, regr = regr, maxOffer =
                       maxOffer, offerOnly = offerOnly, offertol = offertol,
                       reltol = reltol, ...)
    if (results$code) {
        warning("Model fitting did not converge\nMessage: ",
                results$message)
    }

    if (boot > 0) {
        bootMatrix <- stratBoot(boot, report = bootreport, estimate =
                                results$estimate, y = y, a = a, regr = regr, fn
                                = logLikUlt, gr = logLikGradUlt, fixed = fvec,
                                maxOffer = maxOffer, offerOnly = offerOnly,
                                offertol = offertol, reltol = reltol, ...)
    }

    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getStratVcov(results$hessian, fvec)
    ans$log.likelihood <-
        logLikUlt(results$estimate, y = y, acc = a, regr = regr, maxOffer =
                  maxOffer, offertol = offertol, offerOnly = offerOnly)
    ans$call <- cl
    ans$convergence <- list(code = results$code, message = results$message)
    ans$formulas <- formulas
    ans$link <- "logit"
    ans$type <- "private"  ## double check this!
    ans$model <- mf
    ans$y <- y
    ans$equations <- c("R1", "R2", "log(s1)", "log(s2)")
    attr(ans$equations, "hasColon") <- c(TRUE, TRUE, FALSE, FALSE)
    names(attr(ans$equations, "hasColon")) <- ans$equations
    ans$fixed <- fvec
    if (boot > 0)
        ans$boot.matrix <- bootMatrix
    ans$maxOffer <- maxOffer

    class(ans) <- c("strat", "ultimatum")

    return(ans)
}

## m1 <- ultimatum(OffersP + accepts ~ USmaleS + RmaleS | 1, data = usrussia,
##                 maxOffer = 100, s2 = 1, outcome = "offer")

## m2 <- ultimatum(OffersP + accepts ~ USmaleS + RmaleS + whiteS + SlavS +
##                 ScienceS + BusS | 1, data = usrussia, maxOffer = 100, s2 = 1,
##                 outcome = "offer")

## m3 <- ultimatum(OffersP + accepts ~ US_round2 + US_round3 + US_round4 + US_round5
##                + Russia + russia_round2 + russia_round3 + russia_round4 +
##                russia_round5 + RmaleS + USmaleS | US_round2 + US_round3 +
##                US_round4 + US_round5 + Russia + russia_round2 + russia_round3 +
##                russia_round4 + russia_round5 + RmaleS + USmaleS, data =
##                usrussia, maxOffer = 100, s2 = 1, outcome = "offer")
