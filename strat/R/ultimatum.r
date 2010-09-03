##' @include strat.r
NULL

library(maxLik)
library(Formula)
library(foreign)
usrussia <- read.dta("~/.strat/USRussia_Table4B.dta")

source("~/.strat/strat/R/strat.r")

offerCDF <- function(y, maxOffer, fit1, fit2, s1, s2)
{
    1 / (1 + exp((maxOffer - y - fit1 - s2 * (1 + exp((y - fit2) / s2))) / s1))
}

offerPDF <- function(y, maxOffer, fit1, fit2, s1, s2)
{
    num <- exp(-(maxOffer - y - fit1 - s2 * (1 + exp((y - fit2) / s2))) / s1)
    ## num <- finitize(num)
    denom <- s1 * (1 + num)^2
    ## denom <- finitize(denom)
    ans <- (num / denom) * (1 + exp((y - fit2) / s2))
    return(ans)
}

logLikUlt <- function(b, y, acc, regr, maxOffer, offerOnly, ...)
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
    
    ans <- ifelse(y == maxOffer, highball,
                  ifelse(y == 0, lowball, interior))
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
##' @title 
##' @param formulas 
##' @param data 
##' @param subset 
##' @param na.action 
##' @param maxOffer 
##' @param s1 
##' @param s2 
##' @param outcome 
##' @param ... 
##' @param usegrad 
##' @param reltol 
##' @return 
##' @author Brenton Kenkel
ultimatum <- function(formulas, data, subset, na.action,
                      maxOffer,
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
    } else {
        if (!offerOnly) {
            stop("a dependent variable for acceptance must be specified if",
                 " `outcome == \"both\"`; see `?stratult`")
        }
        y <- ya
    }

    regr <- list()
    regr$X <- model.matrix(formulas, data = mf, rhs = 1)
    regr$Z <- model.matrix(formulas, data = mf, rhs = 2)

    s1 <- if (s1null) log(sd(y)) else log(s1)
    s2 <- if (s2null) log(sd(y)) else log(s2)

    ## suppressing warnings in the logit fitting because fitted probabilities
    ## numerically equal to 0/1 seem to occur often
    m2 <- suppressWarnings(glm.fit(regr$Z, a, family = binomial(link = "logit"),
                                   intercept = FALSE, offset = -y))
    m1 <- lsfit(regr$X, maxOffer - y, intercept = FALSE)
    sval <- c(m1$coefficients, m2$coefficients, s1, s2)

    firstTry <- logLikUlt(sval, y = y, acc = a, regr = regr, maxOffer =
                          maxOffer, offerOnly = offerOnly)
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
                       maxOffer, offerOnly = offerOnly, reltol = reltol, ...)
    if (results$code) {
        warning("Model fitting did not converge\nMessage: ",
                results$message)
    }

    if (boot > 0) {
        bootMatrix <- stratBoot(boot, report = bootreport, estimate =
                                results$estimate, y = y, a = a, regr = regr, fn
                                = logLikUlt, gr = logLikGradUlt, fixed = fvec,
                                maxOffer = maxOffer, offerOnly = offerOnly,
                                reltol = reltol, ...)
    }

    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getStratVcov(results$hessian, fvec)
    ans$log.likelihood <-
        logLikUlt(results$estimate, y = y, acc = a, regr = regr, maxOffer =
                  maxOffer, offerOnly = offerOnly)
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
    
    class(ans) <- c("strat", "stratult")

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
