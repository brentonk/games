##' @include strat.r
NULL

##' Makes predicted probabilities from a discrete strategic model
##'
##' This function uses a fitted strategic model to make predictions for a new
##' set of data.  Useful for cross-validating or for graphical analysis.
##' @title Predicted probabilities for a strat12 model
##' @param object a fitted model of class \code{strat12}, usually the output of a
##' call to the \code{\link{strat12}} fitting function
##' @param newdata data frame of values to make the predicted probabilities for.
##' If this is left empty, the original dataset is used.
##' @param probs whether to provide probabilities for outcomes (L, RL, or RR) or
##' for actions (1 moves L or R, 2 moves L or R given that 1 moves R)
##' @param ... other arguments, currently ignored
##' @return data frame of predicted probabilities
##' @method predict strat12
##' @S3method predict strat12
##' @seealso \code{\link{strat12}}
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
predict.strat12 <- function(object, newdata, probs = c("outcome", "action"), ...)
{
    probs <- match.arg(probs)
    
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

    regr <- list()
    for (i in seq_len(length(object$formulas)[2]))
        regr[[i]] <- model.matrix(object$formulas, data = mf, rhs = i)

    ans <- makeProbs12(object$coefficients, regr = regr, link = object$link, type
                      = object$type)
    ans <- do.call(cbind, ans)

    if (probs == "outcome") {
        ans <- data.frame(cbind(ans[, 1], ans[, 2] * ans[, 3],
                                ans[, 2] * ans[, 4]))
        names(ans) <- paste("Pr(", levels(object$y), ")", sep = "")
    }

    ans <- as.data.frame(ans)
    return(ans)
}

sbi12 <- function(y, regr, link)
{
    ## have to do this because binomial() issues warning if it's not directly
    ## passed a character string to its link argument
    if (link == "probit") {
        fam <- binomial(link = "probit")
    } else {
        fam <- binomial(link = "logit")
    }
    
    Z2 <- regr$Z[y != 1, ]
    y2 <- as.numeric(y == 3)[y != 1]
    m2 <- glm.fit(Z2, y2, family = fam)
    p4 <- as.numeric(regr$Z %*% coef(m2))
    p4 <- if (link == "probit") pnorm(p4) else plogis(p4)

    X1 <- cbind(-regr$X1, (1 - p4) * regr$X3, p4 * regr$X4)
    y1 <- as.numeric(y != 1)
    m1 <- glm.fit(X1, y1, family = fam)

    ans <- sqrt(2) * c(coef(m1), coef(m2))
    return(ans)
}

makeUtils12 <- function(b, regr)
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

makeSDs12 <- function(b, regr, type)
{
    sds <- vector("list", 4)
    rcols <- sapply(regr, ncol)

    if (length(rcols) == 5) {  ## i.e., sdByPlayer == FALSE
        v <- exp(as.numeric(regr[[5]] %*% b))
        for (i in 1:4) sds[[i]] <- v
    } else {
        v1 <- exp(as.numeric(regr[[5]] %*% b[1:rcols[5]]))
        v2 <- exp(as.numeric(regr[[6]] %*% b[(rcols[5]+1):length(b)]))
        if (type == "agent") {
            sds[[1]] <- sds[[2]] <- v1
            sds[[3]] <- sds[[4]] <- v2
        } else if (type == "private") {
            sds[[1]] <- sds[[2]] <- sds[[3]] <- v1
            sds[[4]] <- v2
        }
    }

    return(sds)
}

makeProbs12 <- function(b, regr, link, type)
{
    utils <- makeUtils12(b, regr)

    ## length(utils$b) == 0 means no terms left for the variance components
    if (length(utils$b) == 0) {
        sds <- as.list(rep(1, 4))
    } else {
        sds <- makeSDs12(utils$b, regr, type)
    }

    linkfcn <- switch(link,
                      logit = function(x, sd = 1) plogis(x, scale = sd),
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

actionsToOutcomes12 <- function(probs, log.p = TRUE)
{
    probs <- do.call(cbind, probs)
    ans <- cbind(log(probs[, 1]),
                 log(probs[, 2]) + log(probs[, 3]),
                 log(probs[, 2]) + log(probs[, 4]))

    if (!log.p) ans <- exp(ans)
    return(ans)
}

logLik12 <- function(b, y, regr, link, type, ...)
{
    probs <- makeProbs12(b, regr, link, type)
    logProbs <- actionsToOutcomes12(probs, log.p = TRUE)
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

logLikGrad12 <- function(b, y, regr, link, type, ...)
{
    u <- makeUtils12(b, regr)
    p <- makeProbs12(b, regr, link, type)
    rcols <- sapply(regr, ncol)

    if (link == "probit" && type == "private") {
        dp4 <- dnorm(u$u24)
        dgp4 <- dp4 * regr$Z
        Dp4 <- cbind(matrix(0L, nrow = nrow(regr$X1), ncol = sum(rcols[1:3])),
                     dgp4)
        Dp3 <- -Dp4

        dp1 <- dnorm((u$u11 - p$p3 * u$u13 - p$p4 * u$u14) /
                     sqrt(1 + p$p3^2 + p$p4^2))
        dp1 <- dp1 / sqrt(1 + p$p3^2 + p$p4^2)
        dbp1 <- dp1 * cbind(regr$X1, -p$p3 * regr$X3, -p$p4 * regr$X4)
        dgp1 <- (dp1 * dgp4) / sqrt(1 + p$p3^2 + p$p4^2)
        dgp1 <- dgp1 * ((u$u13 - u$u14) * sqrt(1 + p$p3^2 + p$p4^2) -
                        (u$u11 - p$p3*u$u13 - p$p4*u$u14) *
                        ((p$p4 - p$p3) / sqrt(1 + p$p3^2 + p$p4^2)))
        Dp1 <- cbind(dbp1, dgp1)
        Dp2 <- -Dp1
    } else if (type == "agent") {
        derivCDF <- switch(link,
                           logit = dlogis,
                           probit = dnorm)
        
        dp4 <- derivCDF(u$u24 / sqrt(2)) / sqrt(2)
        Dp4 <- cbind(matrix(0L, nrow = nrow(regr$X1), ncol = sum(rcols[1:3])),
                     dp4 * regr$Z)
        Dp3 <- -Dp4

        dp1 <- derivCDF((u$u11 - p$p3 * u$u13 - p$p4 * u$u14) / sqrt(2)) /
            sqrt(2)
        dbp1 <- dp1 * cbind(regr$X1, -p$p3 * regr$X3, -p$p4 * regr$X4)
        dgp1 <- dp1 * dp4 * (u$u13 - u$u14) * regr$Z
        Dp1 <- cbind(dbp1, dgp1)
        Dp2 <- -Dp1
    }

    dL1 <- Dp1 / p$p1
    dL3 <- Dp2 / p$p2 + Dp3 / p$p3
    dL4 <- Dp2 / p$p2 + Dp4 / p$p4

    ans <- as.numeric(y == 1) * dL1 + as.numeric(y == 2) * dL3 +
        as.numeric(y == 3) * dL4
    return(ans)
}

makeResponse12 <- function(yf)
{
    if (length(dim(yf))) {              # response specified as dummies
        Y <- yf
        if (ncol(Y) > 2)
            warning("only first two columns of response will be used")
        
        Y <- Y[, 1:2]
        if (!identical(sort(unique(unlist(yf))), c(0, 1)))
            stop("dummy responses must be dummy variables")
        
        y <- numeric(nrow(Y))
        y[Y[, 1] == 0] <- 1
        y[Y[, 1] == 1 & Y[, 2] == 0] <- 2
        y[Y[, 1] == 1 & Y[, 2] == 1] <- 3
        yf <- as.factor(y)
        levels(yf) <- c(paste("~", names(Y)[1], sep = ""),
                        paste(names(Y)[1], ",~", names(Y)[2], sep = ""),
                        paste(names(Y)[1], ",", names(Y)[2], sep = ""))
    } else {
        yf <- as.factor(yf)
        if (nlevels(yf) != 3) stop("dependent variable must have 3 values")
    }

    return(yf)
}

##' Fits a strategic model with three terminal nodes, as in the game illustrated
##' below in \dQuote{Details}.
##'
##' The model corresponds to the following extensive-form game, described in
##' Signorino (2003):
##' \preformatted{
##' .     1
##' .     /\
##' .    /  \
##' .   /    \ 2
##' .  u11   /\
##' .       /  \
##' .      /    \
##' .    u13    u14
##' .    0      u24}
##' 
##' If Player 1 chooses L, the game ends and Player 1 receives payoffs of u11.
##' (Player 2's utilities in this case cannot be identified in a statistical
##' model.)  If Player 1 chooses L, then Player 2 can choose L, resulting in
##' payoffs of u13 for Player 1 and 0 for Player 2, or R, with payoffs of u14
##' for 1 and u24 for 2.
##'
##' The four formulas specified in the function's \code{formulas} argument
##' correspond to the regressors to be placed in u11, u13, u14, and u24
##' respectively.  If there is any regressor (including the constant) placed in
##' all of u11, u13, and u14, \code{strat12} will stop and issue an error
##' message, because the model is then unidentified (see Lewis and Schultz
##' 2003).  There are two equivalent ways to express the formulas passed to this
##' argument.  One is to use a list of four formulas, where the first contains
##' the response variable(s) (discussed below) on the left-hand side and the
##' other three are one-sided.  For instance, suppose:
##' \itemize{
##' \item u11 is a function of \code{x1}, \code{x2}, and a constant
##' \item u13 is set to 0
##' \item u14 is a function of \code{x3} and a constant
##' \item u24 is a function of \code{z} and a constant.}
##' The list notation would be \code{formulas = list(y ~ x1 + x2, ~ 0, ~ x3, ~
##' z)}.  The other method is to use the \code{\link{Formula}} syntax, with one
##' left-hand side and four right-hand sides (separated by vertical bars).  This
##' notation would be \code{formulas = y ~ x1 + x2 | 0 | x3 | z}.
##'
##' There are three equivalent ways to specify the outcome in \code{formulas}.
##' One is to use a numeric vector with three unique values, with their values
##' (from lowest to highest) corresponding with the terminal nodes of the game
##' tree illustrated above (from left to right).  The second is to use a factor,
##' with the levels (in order as given by \code{levels(y)}) corresponding to the
##' terminal nodes.  The final way is to use two indicator variables, with the
##' first standing for whether Player 1 moves L (0) or R (1), the second
##' standing for Player 2's choice if Player 1 moves R.  (The values of the
##' second when Player 1 moves L should be set to 0 or 1, \strong{not}
##' \code{NA}, in order to ensure that observations are not dropped from the
##' data when \code{na.action = na.omit}.)  The way to specify \code{formulas}
##' when using indicators variables is, for example, \code{y1 + y2 ~ x1 + x2 | 0
##' | x3 | z}.
##' @title Fit a strategic model with 3 terminal nodes
##' @param formulas a list of four formulas, or a \code{Formula} object with
##' four right-hand sides.  See \dQuote{Details} and \dQuote{Examples}.
##' @param data a data frame
##' @param subset an optional logical vector specifying which observations from
##' \code{data} to use in fitting
##' @param na.action how to deal with \code{NA}s in \code{data}.  Defaults to
##' the \code{na.action} setting of \code{\link{options}}.  See
##' \code{\link{na.omit}}
##' @param varformulas an optional list of four formulas or \code{Formula}
##' object with four right-hand sides.  If omitted, all error terms are assumed
##' to have scale parameter 1.  See \dQuote{Details}.
##' @param fixedUtils numeric vector of values to fix for u11, u13, u14, and u24
##' respectively.  \code{NULL} (the default) indicates that these should be
##' estimated with regressors rather than fixed
##' @param link whether to use a probit (default) or logit link structure
##' @param type whether to use an agent-error (\dQuote{agent}, default) or
##' private-information (\dQuote{private}) stochastic structure
##' @param startvals whether to get starting values for the optimization from
##' statistical backwards induction (\dQuote{sbi}, default), from a uniform
##' distribution (\dQuote{unif}), or to set them all to 0 (\dQuote{zero})
##' @param boot number of bootstrap iterations to perform (if any)
##' @param bootreport logical: whether to print status bar when performing
##' bootstrap iterations
##' @param ... other arguments to pass to the fitting function (see
##' \code{\link{maxBFGS}})
##' @return An object of class \code{c("strat", "strat12")}. A
##' \code{strat} object is a list containing: \describe{
##' \item{\code{coefficients}}{estimated parameters of the model}
##' \item{\code{vcov}}{estimated variance-covariance matrix}
##' \item{\code{log.likelihood}}{vector of individual log likelihoods (left
##' unsummed for use with non-nested model tests)}
##' \item{\code{call}}{the call used to produce the model}
##' \item{\code{formulas}}{the final \code{Formula} object passed to
##' \code{model.frame}}
##' \item{\code{link}}{the link function used}
##' \item{\code{type}}{the stochastic structure used}
##' \item{\code{model}}{the model frame containing all variables used in fitting}
##' \item{\code{y}}{the dependent variable, represented as a factor}
##' }
##' The second class of the returned object, \code{strat12}, is for use with the
##' \code{predict} method.
##' @seealso \code{\link{summary.strat}} and \code{\link{predict.strat12}} for
##' postestimation analysis; \code{\link{Formula}} for formula specification
##' @export
##' @references Jeffrey B. Lewis and Kenneth A Schultz.  2003.
##' \dQuote{Revealing Preferences: Empirical Estimation of a Crisis Bargaining
##' Game with Incomplete Information.}  \emph{Political Analysis} 11:345--367.
##'
##' Curtis S. Signorino.  2003.  \dQuote{Structure and Uncertainty
##' in Discrete Choice Models.}  \emph{Political Analysis} 11:316--344.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
strat12 <- function(formulas, data, subset, na.action,
                    link = c("probit", "logit"),
                    type = c("agent", "private"),
                    startvals = c("sbi", "unif", "zero"),
                    fixedUtils = NULL,
                    sdformula = NULL,
                    sdByPlayer = FALSE,
                    boot = 0,
                    bootreport = TRUE,
                    ...)
{
    cl <- match.call()

    link <- match.arg(link)
    type <- match.arg(type)
    startvals <- match.arg(startvals)

    ## various sanity checks
    formulas <- checkFormulas(formulas)
    if (is.null(fixedUtils) && length(formulas)[2] != 4)
        stop("`formulas` should have four components on the right-hand side")

    if (!is.null(fixedUtils)) {
        if (length(fixedUtils) < 4)
            stop("fixedUtils must have 4 elements (u11, u13, u14, u24)")
        if (length(fixedUtils) > 4) {
            warning("only the first 4 elements of fixedUtils will be used")
            fixedUtils <- fixedUtils[1:4]
        }
        
        formulas <- update(formulas, . ~ 1 | 1 | 1 | 1)

        if (startvals == "sbi")
            startvals <- "zero"

        if (is.null(sdformula))
            sdformula <- if (sdByPlayer) Formula(~ 1 | 1) else Formula(~ 1)
    }

    if (!is.null(sdformula)) {
        sdformula <- checkFormulas(sdformula, argname = "sdformula")
        if (sdByPlayer && length(sdformula)[2] != 2)
            stop("`sdformula` should have two components (one for each player) on the right-hand side when sdByPlayer == TRUE")
        if (!sdByPlayer && length(sdformula)[2] != 1)
            stop("`sdformula` should have exactly one component on the right-hand side")
        formulas <- as.Formula(formula(formulas), formula(sdformula))
    }

    if (sdByPlayer && is.null(sdformula)) {
        warning("to estimate SDs by player, you must specify `sdformula` or `fixedUtils`")
        sdByPlayer <- FALSE
    }

    if (link == "logit" && type == "private") {
        warning("logit link cannot be used with private information model; changing to probit link")
        link <- "probit"
    }

    ## make the model frame
    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
    yf <- makeResponse12(yf)
    y <- as.numeric(yf)

    ## makes a list of the 4 (or 8, if variance formulas specified) matrices of
    ## regressors to be passed to estimation functions
    regr <- list()
    for (i in seq_len(length(formulas)[2]))
        regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
    names(regr) <- character(length(regr))
    names(regr)[1:4] <- c("X1", "X3", "X4", "Z")
    rcols <- sapply(regr, ncol)

    ## makes starting values -- specify "unif" (a numeric vector of length two)
    ## to control the interval from which uniform values are drawn
    if (startvals == "zero") {
        sval <- rep(0, sum(rcols))
    } else if (startvals == "unif") {
        if (!hasArg(unif))
            unif <- c(-1, 1)
        sval <- runif(sum(rcols), unif[1], unif[2])
    } else {
        sval <- sbi12(y, regr, link)
        sval <- c(sval, rep(0, sum(rcols) - length(sval)))
    }

    ## identification check
    varNames <- sapply(regr, colnames)
    idCheck <- (varNames[[1]] %in% varNames[[2]])
    idCheck <- idCheck & (varNames[[1]] %in% varNames[[3]])
    if (is.null(fixedUtils) && any(idCheck)) {
        stop("Identification problem: the following variables appear in all three of player 1's utility equations: ",
             paste(varNames[[1]][idCheck], collapse = ", "))
    }

    ## gives names to starting values
    prefixes <- paste(c(rep("u1(", 3), "u2("), c(levels(yf), levels(yf)[3]),
                       ")", sep = "")
    sdterms <- if (!is.null(sdformula)) { if (sdByPlayer) 2L else 1L } else 0L
    varNames <- makeVarNames(varNames, prefixes, link, sdterms)
    hasColon <- varNames$hasColon
    names(sval) <- varNames$varNames

    ## use the gradient iff no regressors in variance
    gr <- if (is.null(sdformula)) logLikGrad12 else NULL

    ## deals with fixed utilities
    fvec <- rep(FALSE, length(sval))
    names(fvec) <- names(sval)
    if (!is.null(fixedUtils)) {
        sval[1:4] <- fixedUtils
        fvec[1:4] <- TRUE
    }

    results <- maxBFGS(fn = logLik12, grad = gr, start = sval, fixed = fvec, y =
                       y, regr = regr, link = link, type = type, ...)
    if (results$code) {
        warning("Model fitting did not converge\nMessage: ",
                results$message)
    }

    if (boot > 0) {
        bootMatrix <- stratBoot(boot, report = bootreport, estimate =
                                results$estimate, y = y, regr = regr, fn =
                                logLik12, gr = gr, fixed = fvec, link = link,
                                type = type, ...)
    }

    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getStratVcov(results$hessian, fvec)
    ans$log.likelihood <- logLik12(results$estimate, y = y, regr = regr, link =
                                   link, type = type)
    ans$call <- cl
    ans$convergence <- list(code = results$code, message = results$message)
    ans$formulas <- formulas
    ans$link <- link
    ans$type <- type
    ans$model <- mf
    ans$y <- yf
    ans$equations <- prefixes
    attr(ans$equations, "hasColon") <- hasColon
    ans$fixed <- fvec

    if (boot > 0)
        ans$boot.matrix <- bootMatrix
    
    class(ans) <- c("strat", "strat12")
    
    return(ans)
}
