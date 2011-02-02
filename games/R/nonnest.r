##' @include games.r
##' @include helpers.r
NULL

print.nonnest.test <- function(x, digits = x$digits, ...)
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
## Individual log-likelihoods for a model of class "game", "lm", or "glm"
##
## The "outcome" argument is for comparing a strategic model to a binary logit
## or probit model where the response is one of the possible outcomes (see the
## examples for the "vuong"/"clarke" functions)
## 
indivLogLiks <- function(model, outcome = NULL)
{
    ## get weights (only relevant for lm and glm models)
    if (is.null(wt <- weights(model)))
        wt <- rep(1L, nobs(model))
    
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
        ans <- wt * log(ans)
    } else if (inherits(model, "lm") && !inherits(model, "glm")) {
        msum <- summary(model)
        res <- msum$residuals
        sigma <- msum$sigma
        ans <- log(sqrt(wt) * dnorm(res, sd = sigma))
    } else {
        stop("model is not of a supported class")
    }

    return(ans)
}

nparams <- function(model)
{
    if (inherits(model, "game")) {
        ans <- sum(!model$fixed)
    } else if (inherits(model, "glm")) {
        ans <- attr(logLik(model), "df")
    } else if (inherits(model, "lm")) {
        ans <- length(coef(model))
    }

    return(ans)
}

nobs <- function(model)
{
    if (inherits(model, "game")) {
        ans <- nrow(model$model)
    } else {
        ans <- length(model$residuals)
    }

    return(ans)
}

gety <- function(model, outcome)
{
    if (inherits(model, "ultimatum")) {
        if (!is.null(outcome)) {
            outcome <- c(offer = 1, accept = 2)[outcome]
            ans <- model$y[, outcome]
        } else {
            ans <- model$y
        }
    } else if (inherits(model, "game")) {
        if (!is.null(outcome)) {
            ans <- as.numeric(as.numeric(model$y) == outcome)
        } else {
            ans <- as.numeric(model$y)
        }
    } else if (inherits(model, c("lm", "glm"))) {
        ans <- model$y
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
        stop("model1 and model2 have different numbers of observations")
    
    y1 <- gety(model1, outcome1)
    y2 <- gety(model2, outcome2)

    ## check for equality of dependent variables
    if (length(dim(y1))) {  ## ultimatum model with offer and accept
        if (!length(dim(y2))) {
            stop("models do not have same dependent variable")
        } else {
            if (!isTRUE(all.equal(y1[, 1], y2[, 1])))
                stop("models do not have same dependent variable")
            if (!isTRUE(all.equal(y1[, 2], y2[, 2])))
                stop("models do not have same dependent variable")
        }
    } else {
        if (!isTRUE(all.equal(y1, y2)))
            stop("models do not have same dependent variable")
    }

    if (is.null(w1 <- weights(model1)))
        w1 <- rep(1L, n)
    if (is.null(w2 <- weights(model2)))
        w2 <- rep(1L, n)
    if (!isTRUE(all.equal(w1, w2)))
        stop("model1 and model2 have different weights")
    
    loglik1 <- indivLogLiks(model1, outcome1)
    loglik2 <- indivLogLiks(model2, outcome2)
    p1 <- nparams(model1)
    p2 <- nparams(model2)

    return(list(n = n, loglik1 = loglik1, loglik2 = loglik2, p1 = p1, p2 = p2))
}

##' Perform Vuong's (1989) or Clarke's (2007) test for non-nested model
##' selection.
##'
##' These tests are for comparing two statistical models that have the same
##' dependent variable, where neither model can be expressed as a special case
##' of the other (i.e., they are non-nested).  The null hypothesis is that the
##' estimated models are the same Kullback-Leibler distance from the true
##' model.  To adjust for potential differences in the dimensionality of the
##' models, the test statistic for both \code{vuong} and \code{clarke} is
##' corrected using the Bayesian information criterion (see Clarke 2007 for
##' details).
##'
##' It is crucial that the dependent variable be exactly the same between the
##' two models being tested, including the order the observations are placed
##' in.  Weights, if any are used, must also be the same between models.  The
##' \code{vuong} and \code{clarke} functions check for such discrepancies, and
##' stop with an error if any is found.
##'
##' When comparing a strategic model to a (generalized) linear model, you must
##' take care to ensure that the dependent variable is truly the same between
##' models.  This is where the \code{outcome} arguments come into play.  For
##' example, in an \code{\link{ultimatum}} model where acceptance is observed,
##' the dependent variable for each observation is the vector consisting of the
##' offer size and an indicator for whether it was accepted.  This is not the
##' same as the dependent variable in a least-squares regression of offer size,
##' which is a scalar for each observation.  Therefore, for a proper comparison
##' of \code{model1} of class {"ultimatum"} and \code{model2} of class
##' \code{"lm"}, it is necessary to specify \code{outcome1 = "offer"}.
##' Similarly, consider an \code{\link{egame12}} model on the
##' \code{\link{war1800}} data, where player 1 chooses whether to escalate the
##' crisis and player 2 chooses whether to go to war.  The dependent variable
##' for each observation in this model is the vector of each player's choice.
##' By contrast, in a logistic regression where the dependent variable is
##' whether war occurs, the dependent variable for each observation is a
##' scalar.  To compare these models, it is necessary to specify \code{outcome1
##' = 3}.
##' @title Non-nested model tests
##' @aliases vuong clarke
##' @usage vuong(model1, model2, outcome1=NULL, outcome2=NULL, level=0.05, digits=2)
##' clarke(model1, model2, outcome1=NULL, outcome2=NULL, level=0.05, digits=2)
##' @param model1 A fitted statistical model of class \code{"game"},
##' \code{"lm"}, or \code{"glm"}
##' @param model2 A fitted statistical model of class \code{"game"},
##' \code{"lm"}, or \code{"glm"} whose dependent variable is the same as that of
##' \code{model1}
##' @param outcome1 Optional: if \code{model1} is of class \code{"game"},
##' specify an integer to restrict attention to a particular binary outcome (the
##' corresponding column of \code{predict(model1)}).  For
##' \code{\link{ultimatum}} models, "offer" or "accept" may also be used.  See
##' "Details" below for more information on when to specify an outcome.  If
##' \code{model1} is not of class \code{"game"} and \code{outcome1} is
##' non-\code{NULL}, it will be ignored and a warning will be issued.
##' @param outcome2 Optional: same as \code{outcome1}, but corresponding to
##' \code{model2}.
##' @param level Numeric: significance level for the test.
##' @param digits Integer: number of digits to print
##' @references Quang H. Vuong.  1989.  "Likelihood Ratio Tests for Model
##' Selection and Non-Nested Hypotheses."  \emph{Econometrica} 57(2): 307--333.
##'
##' Kevin Clarke.  2007.  "A Simple Distribution-Free Test for Nonnested
##' Hypotheses."  \emph{Political Analysis} 15(3): 347--363.
##' @return Typical use will be to run the function interactively and examine
##' the printed output.  The functions return an object of class
##' \code{"nonnest.test"}, which is a list containing: \describe{
##' \item{\code{stat}}{The test statistic}
##' \item{\code{test}}{The type of test (\code{"vuong"} or \code{"clarke"})}
##' \item{\code{level}}{Significance level for the test}
##' \item{\code{digits}}{Number of digits to print}
##' \item{\code{loglik1}}{Vector of observationwise log-likelihoods for
##' \code{model1}}
##' \item{\code{loglik2}}{Vector of observationwise log-likelihoods for
##' \code{model2}}
##' \item{\code{nparams}}{Integer vector containing the number of parameters
##' fitted in \code{model1} and \code{model2} respectively}
##' \item{\code{nobs}}{Number of observations of the dependent variable being
##' modeled}}
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @examples
##' data(war1800)
##'
##' ## balance of power model
##' f1 <- esc + war ~ revis1 + s_wt_re1 | 0 | balanc | revis2 + balanc
##' m1 <- egame12(f1, data = war1800, subset = !is.na(regime1) & !is.na(regime2))
##'
##' ## regime type model
##' f2 <- esc + war ~ regime1 + s_wt_re1 | 0 | 1 | regime2
##' m2 <- egame12(f2, data = war1800)
##'
##' 
vuong <- function(model1, model2, outcome1 = NULL, outcome2 = NULL,
                  level = 0.05, digits = 2)
{
    x <- nonnest(model1, model2, outcome1, outcome2)

    correction <- (x$p1 - x$p2) * (log(x$n) / 2)
    num <- sum(x$loglik1) - sum(x$loglik2) - correction
    denom <- sd(x$loglik1 - x$loglik2) *
        sqrt(x$n / (x$n-1))  # correcting for R's use of n-1 denominator in sd
                             # calculation
    stat <- num / (sqrt(x$n) * denom)

    ans <- list(stat = stat,
                test = "vuong",
                level = level,
                digits = digits,
                loglik1 = x$loglik1,
                loglik2 = x$loglik2,
                nparams = c(x$p1, x$p2),
                nobs = x$n)
    class(ans) <- "nonnest.test"

    return(ans)
}

clarke <- function(model1, model2, outcome1 = NULL, outcome2 = NULL,
                   level = 0.05, digits = 2)
{
    x <- nonnest(model1, model2, outcome1, outcome2)

    correction <- (x$p1 - x$p2) * (log(x$n) / (2*x$n))
    stat <- sum(x$loglik1 - x$loglik2 > correction)

    ans <- list(stat = stat,
                test = "clarke",
                level = level,
                digits = digits,
                loglik1 = x$loglik1,
                loglik2 = x$loglik2,
                nparams = c(x$p1, x$p2),
                nobs = x$n)
    class(ans) <- "nonnest.test"

    return(ans)
}
