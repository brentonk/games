##' A package for estimating strategic statistical models.
##' 
##' @name strat-package
##' @docType package
##' @examples
##' mean(1)
NULL

##' The default method for printing a \code{strat} object.
##'
##' Prints the call and coefficients of a fitted strategic model.
##' @title Print a strategic model object
##' @method print strat
##' @param x a fitted model of class \code{strat}
##' @param ... other arguments, currently ignored
##' @return \code{x}, invisibly
##' @S3method print strat
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
print.strat <- function(x, ...)
{
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coef)
    cat("\n")
    invisible(x)
}

##' The default method for summarizing a \code{strat} object.
##'
##' Forms the standard regression results table from a fitted strategic model.
##' Normally used interactively, in conjunction with
##' \code{\link{print.summary.strat}}.
##' @title Summarize a strategic model object
##' @method summary strat
##' @param object a fitted model of class \code{strat}
##' @param ... other arguments, currently ignored
##' @return an object of class \code{summary.strat}, containing the coefficient
##' matrix and other information needed for printing
##' @seealso \code{\link{print.summary.strat}}
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
summary.strat <- function(object, ...)
{
    cf <- object$coefficients
    se <- sqrt(diag(object$vcov))
    zval <- cf / se
    pval <- 2 * pnorm(-abs(zval))

    ans <- list()
    ans$coefficients <- cbind(cf, se, zval, pval)
    colnames(ans$coefficients) <- c("Estimate", "Std. Error", "z value",
                                    "Pr(>|z|)")
    ans$call <- object$call
    ans$log.likelihood <- sum(object$log.likelihood)
    ans$nobs <- nrow(object$model)
    class(ans) <- "summary.strat"

    return(ans)
}

##' <description>
##'
##' Prints the standard regression results table from a fitted strategic model,
##' along with the log-likelihood, AIC, and number of observations.
##' @title print summary
##' @param x an object of class \code{summary.strat}, typically produced by
##'running \code{summary} on a fitted model of class \code{strat}
##' @param ... other arguments, currently ignored
##' @return \code{x}, invisibly.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
print.summary.strat <- function(x, ...)
{
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients)
    cat("\nLog-likelihood:", x$log.likelihood)
    cat("\nAIC:", AIC(x))
    cat("\nNo. observations:", x$nobs, "\n\n")
    invisible(x)
}

##' Get coefficients from a strategic model
##'
##' Extracts the coefficient vector from a fitted model of class \code{strat}.
##' @title Coefficients of a strategic model object
##' @method coef strat
##' @param object a fitted model of class \code{strat}
##' @param ... other arguments, currently ignored
##' @return numeric vector containing coefficients
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
coef.strat <- function(object, ...)
{
    object$coefficients
}

##' Get covariance matrix from a strategic model
##'
##' Extracts the covariance matrix from a fitted model of class \code{strat}.
##' @title Covariance matrix of a strategic model object
##' @method vcov strat
##' @param object a fitted model of class \code{strat}
##' @param ... other arguments, currently ignored
##' @return covariance matrix
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
vcov.strat <- function(object, ...)
{
    object$vcov
}

logLik.strat <- function(object, ...)
{
    ans <- sum(object$log.likelihood)
    attr(ans, "df") <- length(object$coefficients)
    attr(ans, "nobs") <- nrow(object$model)
    class(ans) <- "logLik"
    return(ans)
}

logLik.summary.strat <- function(object, ...)
{
    ans <- object$log.likelihood
    attr(ans, "df") <- nrow(object$coefficients)
    attr(ans, "nobs") <- object$nobs
    class(ans) <- "logLik"
    return(ans)
}

predict.strat <- function(object, newdata, probs = c("outcome", "action"), ...)
{
    probs <- match.arg(probs)
    
    if (object$game != "strat3")
        stop("predict not yet implemented for strat models other than strat3")

    if (missing(newdata))
        newdata <- object$model

    mf <- match(c("subset", "na.action"), names(object$call), 0L)
    mf <- object$call[c(1L, mf)]
    mf$formula <- object$formulas
    mf$data <- newdata
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    regr <- list()
    for (i in seq_len(length(object$formulas)[2]))
        regr[[i]] <- model.matrix(object$formulas, data = mf, rhs = i)

    ans <- makeProbs3(object$coefficients, regr = regr, link = object$link, type
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

finiteProbs <- function(x)
{
    x <- replace(x, x < .Machine$double.eps, .Machine$double.eps)
    x <- replace(x, x > 1 - .Machine$double.neg.eps,
                 1 - .Machine$double.neg.eps)
    return(x)
}

sbi3 <- function(y, regr, link)
{
    ## have to do this because binomial() issues warning if it's not directly
    ## passed a character string to its family argument
    if (link == "probit") {
        fam <- binomial(link = "probit")
    } else {
        fam <- binomial(link = "logit")
    }
    
    Z2 <- regr$Z[y != 1, ]
    y2 <- as.numeric(y == 3)[y != 1]
    m2 <- glm.fit(Z2, y2, family = fam)
    class(m2) <- "glm"
    p4 <- as.numeric(regr$Z %*% coef(m2))
    p4 <- if (link == "probit") pnorm(p4) else plogis(p4)

    X1 <- cbind(-regr$X1, (1 - p4) * regr$X3, p4 * regr$X4)
    y1 <- as.numeric(y != 1)
    m1 <- glm.fit(X1, y1, family = fam)

    ans <- sqrt(2) * c(coef(m1), coef(m2))
    return(ans)
}

makeUtils3 <- function(b, regr)
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

makeSDs3 <- function(b, regr)
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

actionsToOutcomes <- function(probs, game, log = TRUE)
{
    probs <- do.call(cbind, probs)
    if (game == "strat3") {
        ans <- cbind(base::log(probs[, 1]),
                     base::log(probs[, 2]) + base::log(probs[, 3]),
                     base::log(probs[, 2]) + base::log(probs[, 4]))
    }

    if (!log) ans <- exp(ans)
    return(ans)
}

makeProbs3 <- function(b, regr, link, type)
{
    utils <- makeUtils3(b, regr)

    ## length(utils$b) == 0 means no terms left for the variance components
    if (length(utils$b) == 0) {
        sds <- list(1, 1, 1, 1)
    } else {
        sds <- makeSDs3(utils$b, regr)
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

logLik3 <- function(b, y, regr, link, type)
{
    probs <- makeProbs3(b, regr, link, type)
    logProbs <- actionsToOutcomes(probs, game = "strat3", log = TRUE)
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

logLikGrad3 <- function(b, y, regr, link, type)
{
    u <- makeUtils3(b, regr)
    p <- makeProbs3(b, regr, link, type)
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

##' Fits a strategic model with three terminal nodes, as in the game illustrated
##' below in \dQuote{Details}.
##'
##' The model corresponds to the following extensive-form game, described in
##' Signorino (2003):
##' \preformatted{
##' .     /\
##' .    /  \
##' .   /    \
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
##' all of u11, u13, and u14, \code{strat3} will stop and issue an error
##' message, because the model is then unidentified (see Lewis and Schultz
##' 2003).  There are two equivalent ways to express the formulas passed to this
##' argument.  One is to use a list of four formulas, where the first contains
##' the response variable(s) (discussed below) on the left-hand side and the
##' other three are one-sided.  For instance, suppose u11 is a function of
##' \code{x1}, \code{x2}, and a constant; u13 is set to 0; u14 is a function of
##' \code{x3} and a constant; and u24 is a function of \code{z} and a constant.
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
##' ?\code{\link{na.omit}}
##' @param varformulas an optional list of four formulas or \code{Formula}
##' object with four right-hand sides.  If omitted, all error terms are assumed
##' to have scale parameter 1.  See \dQuote{Details}.
##' @param link whether to use a probit (default) or logit link structure
##' @param type whether to use an agent-error (\dQuote{agent}, default) or
##' private-information (\dQuote{private}) stochastic structure
##' @param startvals whether to get starting values for the optimization from
##' statistical backwards induction (\dQuote{sbi}, default), from a uniform
##' distribution (\dQuote{unif}), or to set them all to 0 (\dQuote{zero})
##' @param ... other arguments to pass to the fitting function
##' @return An object of class \code{c("strat", "strat3")}. A
##' \code{strat} object is a list containing: \describe{
##' \item{\code{coefficients}}{estimated parameters of the model}
##' \item{\code{vcov}}{estimated variance-covariance matrix}
##' \item{\code{log.likelihood}}{vector of individual log likelihoods (left unsummed
##' for use with non-nested model tests)}
##' \item{\code{call}}{the call used to produce the model}
##' \item{\code{formulas}}{the final \code{Formula} object passed to
##' \code{model.frame}}
##' \item{\code{link}}{the link function used}
##' \item{\code{type}}{the stochastic structure used}
##' \item{\code{model}}{the model frame containing all variables used in fitting}
##' \item{\code{y}}{the dependent variable, represented as a factor}
##' }
##' The second class, \code{strat3}, is for use with the \code{predict} method.
##' @seealso \code{\link{summary.strat}}, \code{\link{Formula}}
##' @export
##' @references Jeffrey B. Lewis and Kenneth A Schultz.  2003.
##' \dQuote{Revealing Preferences: Empirical Estimation of a Crisis Bargaining
##' Game with Incomplete Information.}  \emph{Political Analysis} 11:345--367.
##'
##' Curtis S. Signorino.  2003.  \dQuote{Structure and Uncertainty
##' in Discrete Choice Models.}  \emph{Political Analysis} 11:316--344.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
strat3 <- function(formulas, data, subset, na.action,
                    varformulas,
                    link = c("probit", "logit"),
                    type = c("agent", "private"),
                    startvals = c("sbi", "unif", "zero"),
                    ...)
{
    cl <- match.call()

    link <- match.arg(link)
    type <- match.arg(type)
    startvals <- match.arg(startvals)

    ## various sanity checks
    if (inherits(formulas, "list")) {
        formulas <- do.call(as.Formula, formulas)
    } else if (inherits(formulas, "formula")) {
        formulas <- as.Formula(formulas)
    } else {
        stop("formulas must be a list of formulas or a formula")
    }

    if (!missing(varformulas))
    {
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

    ## make the model frame
    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
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
        y <- as.numeric(yf)
    }

    ## makes a list of the 4 (or 8, if variance formulas specified) matrices of
    ## regressors to be passed to estimation functions
    regr <- list()
    for (i in seq_len(length(formulas)[2]))
        regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
    names(regr) <- c("X1", "X3", "X4", "Z", "V1", "V2", "V3",
                     "V4")[1:length(regr)]    
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
        sval <- sbi3(y, regr, link)
        sval <- c(sval, rep(0, sum(rcols) - length(sval)))
    }

    ## identification check
    varNames <- sapply(regr, colnames)
    if (any(idCheck <- (varNames[[1]] %in% varNames[[2]]) &
            (varNames[[1]] %in% varNames[[3]]))) {
        stop("Identification problem: the following variables appear in all ",
             "three of player 1's utility equations: ",
             paste(varNames[[1]][idCheck], collapse = ", "))
    }

    ## gives names to starting values
    prefixes <- paste(c(rep("u1(", 3), "u2("), c(levels(yf), levels(yf)[3]),
                       ")", sep = "")
    prefixes <- c(prefixes, "v1", "v2", "v3", "v4")
    for (i in seq_len(length(formulas)[2])) {
        if (length(varNames[[i]]))
            varNames[[i]] <- paste(prefixes[i], varNames[[i]], sep = ":")
    }
    names(sval) <- unlist(varNames)

    ## use the gradient iff no regressors in variance
    gr <- if (missing(varformulas)) logLikGrad3 else NULL

    results <- maxBFGS(fn = logLik3, grad = gr, start = sval, y = y, regr =
                       regr, link = link, type = type, ...)

    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- solve(-results$hessian)
    ans$log.likelihood <- logLik3(results$estimate, y = y, regr = regr, link =
                                  link, type = type)
    ans$call <- cl
    ans$formulas <- formulas
    ans$link <- link
    ans$type <- type
    ans$model <- mf
    ans$y <- yf
    class(ans) <- c("strat", "strat3")
    
    return(ans)
}
