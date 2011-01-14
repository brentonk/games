##' @include games.r
##' @include helpers.r
NULL

##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title print
##' @param x wef
##' @param digits qed
##' @param ... qewd
##' @return qdw
##' @author Brenton Kenkel
print.nonnest.test <- function(x, digits = 2, ...)
{
    if (x$test == "vuong") {
        p <- 2 * pnorm(-abs(x$stat))
        pref <- if (x$stat > 0) 1 else 2
    } else if (x$test == "clarke") {
        b <- min(x$stat, x$nobs - x$stat)
        p <- 2 * pbinom(b, x$nobs, 0.5)
        pref <- if (x$stat > x$nobs - x$stat) 1 else 2
    } else {
        stop("unknown test type")
    }

    testname <- switch(x$test, vuong = "Vuong", clarke = "Clarke")
    cat("\n", testname, " test for non-nested models\n", sep = "")

    cat("\nModel 1 log-likelihood:", format(sum(x$loglik1), digits = digits))
    cat("\nModel 2 log-likelihood:", format(sum(x$loglik2), digits = digits))
    cat("\nObservations:", x$nobs)
    cat("\nTest statistic:", format(x$stat, digits = digits))
    if (x$test == "clarke")
        cat(" (", round(100* x$stat / x$nobs), "%)", sep = "")
    cat("\n")

    fp <- format.pval(p, digits = digits)
    if (substr(fp, 1L, 1L) != "<")
        fp <- paste("=", fp)
    if (p < x$level) {
        cat("\nModel ", pref, " is preferred (p ", fp, ")\n\n", sep = "")
    } else {
        cat("\nNeither model is significantly preferred (p ", fp, ")\n\n",
            sep = "")
    }

    invisible(x)
}

##
## Individual log-likelihoods for a model of class "game", "lm", "glm", "vglm"
## (multinomial logit via VGAM or Zelig), "polr" (ordered logit/probit via MASS
## or Zelig)
##
## The "outcome" argument is for comparing a strategic model to a binary logit
## or probit model where the response is one of the possible outcomes (see the
## examples for the "vuong"/"clarke" functions)
## 
indivLogLiks <- function(model, outcome = NULL)
{
    ## add multinomial logit from package mlogit!

    ## allow outcomes to be specified in other models?
    
    if (inherits(model, "game") && is.null(outcome)) {
        ans <- model$log.likelihood
    } else if (inherits(model, "ultimatum")) {
        outcome <- c(offer = 1, accept = 2)[outcome]
        if (outcome == 1) {
            ans <- attr(model$log.likelihood, "offer")
        } else {
            ans <- attr(model$log.likelihood, "accept")
        }
    } else if (inherits(model, "game")) {
        pp <- predict(model)[, outcome]
        ans <- ifelse(as.numeric(model$y) == outcome, pp, 1 - pp)
        ans <- log(ans)
    } else if (inherits(model, "glm") &&
               model$family$family == "binomial") {
        ans <- ifelse(model$y == 1, fitted.values(model),
                      1 - fitted.values(model))
        ans <- log(ans)
    } else if (inherits(model, "lm") && !inherits(model, "glm")) {
        msum <- summary(model)
        ans <- dnorm(msum$residuals, sd = msum$sigma, log = TRUE)
    } else if (inherits(model, "vglm") &&
               "multinomial" %in% model@family@vfamily) {
        if (!all(model@y %in% c(0, 1)))
            stop("multinomial model must be specified in terms of separate observations, not counts")
        y <- apply(model@y, 1, function(x) which(x == 1))
        pp <- fitted.values(model)
        pp <- pp[cbind(1:nrow(pp), y)]  # gets the y[i]'th entry from the i'th
                                        # row of pp
        ans <- log(pp)
    } else if (inherits(model, "polr")) {
        y <- as.numeric(model.response(model$model))
        if (length(y) != model$nobs)
            stop("polr model must be specified in terms of separate observations, not counts")
        pp <- fitted.values(model)
        pp <- pp[cbind(1:nrow(pp), y)]
        ans <- log(pp)
    } else {
        stop("model is not of a supported class")
    }

    return(ans)
}

nparams <- function(model)
{
    ## need to check for the alternative models listed!

    if (inherits(model, "game")) {
        ans <- sum(!model$fixed)
    } else if (inherits(model, "polr")) {
        ans <- model$edf
    } else if (inherits(model, "vglm")) {
        
    } else if (inherits(model, "glm")) {
        ans <- attr(logLik(model), "df")
    } else if (inherits(model, "lm")) {
        ans <- length(coef(model))
    }

    return(ans)
}

nobs <- function(model)
{
    ## need to write methods for the alternative models!

    if (inherits(model, "game")) {
        ans <- nrow(model$model)
    } else if (inherits(model, c("lm", "glm"))) {
        ans <- length(model$residuals)
    } else if (inherits(model, "polr")) {
        ans <- model$n
    }

    return(ans)
}

##
## setup and error checking for clarke and vuong
## 
nonnest <- function(model1, model2, outcome1, outcome2)
{
    if (!inherits(model1, "game") && !is.null(outcome1))
        warning("model1 is not of class \"game\", so outcome1 will be ignored")
    if (!inherits(model2, "game") && !is.null(outcome2))
        warning("model2 is not of class \"game\", so outcome2 will be ignored")

    n <- nobs(model1)
    if (nobs(model2) != n)
        stop("model1 and model2 must have same number of observations")
    
    loglik1 <- indivLogLiks(model1, outcome1)
    loglik2 <- indivLogLiks(model2, outcome2)
    p1 <- nparams(model1)
    p2 <- nparams(model2)

    return(list(n = n, loglik1 = loglik1, loglik2 = loglik2, p1 = p1, p2 = p2))
}

##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title Non-nested model tests
##' @aliases vuong clarke
##' @param model1 A fitted statistical model (see "Details" for the types of
##' models available)
##' @param model2 A fitted statistical model
##' @param outcome1 Optional: if \code{model1} is of class \code{"game"},
##' specify an integer to restrict attention to a particular binary outcome (the
##' corresponding column of \code{predict(model1)}).  See "Details" below for
##' more information on when to specify an outcome.  If \code{model1} is not of
##' class \code{"game"} and \code{outcome1} is non-\code{NULL}, it will be
##' ignored and a warning will be issued.
##' @param outcome2 Optional: same as \code{outcome1}, but corresponding to
##' \code{model2}.
##' @param level Numeric: significance level for the test.
##' @references (cite vuong and kevin!)
##' @return An object of class \code{"nonnest.test"}, which is a list
##' containing: \describe{
##' \item{\code{stat}}{The test statistic}
##' \item{\code{test}}{The type of test (\code{"vuong"} or \code{"clarke"})}
##' \item{\code{level}}{Significance level for the test}
##' \item{\code{loglik1}}{Vector of observationwise log-likelihoods for
##' \code{model1}}
##' \item{\code{loglik2}}{Vector of observationwise log-likelihoods for
##' \code{model2}}
##' \item{\code{nparams}}{Integer vector containing the number of parameters
##' fitted in \code{model1} and \code{model2} respectively}
##' \item{\code{nobs}}{Number of observations of the dependent variable being
##' modeled}}
##' @seealso \code{\link{print.nonnest.test}}
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
vuong <- function(model1, model2, outcome1 = NULL, outcome2 = NULL,
                  level = 0.05)
{
    x <- nonnest(model1, model2, outcome1, outcome2)

    correction <- (x$p1 - x$p2) * (log(x$n) / 2)
    num <- sum(x$loglik1) - sum(x$loglik2) - correction
    denom <- sd(x$loglik1 - x$loglik2) *
        sqrt(x$n / (x$n-1))  # correcting for R's use of n-1 denominator in sd
                             # calculation
    stat <- num / (sqrt(x*n) * denom)

    ans <- list(stat = stat,
                test = "vuong",
                level = level,
                loglik1 = x$loglik1,
                loglik2 = x$loglik2,
                nparams = c(x$p1, x$p2),
                nobs = x$n)
    class(ans) <- "nonnest.test"

    return(ans)
}

clarke <- function(model1, model2, outcome1 = NULL, outcome2 = NULL,
                   level = 0.05)
{
    x <- nonnest(model1, model2, outcome1, outcome2)

    correction <- (x$p1 - x$p2) * (log(x$n) / (2*x$n))
    stat <- sum(x$loglik1 - x$loglik2 > correction)

    ans <- list(stat = stat,
                test = "clarke",
                level = level,
                loglik1 = x$loglik1,
                loglik2 = x$loglik2,
                nparams = c(x$p1, x$p2),
                nobs = x$n)
    class(ans) <- "nonnest.test"

    return(ans)
}
