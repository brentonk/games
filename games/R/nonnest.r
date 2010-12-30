##' @include games.r
##' @include helpers.r
NULL

##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title 
##' @param x 
##' @param digits 
##' @param ... 
##' @return 
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
        y <- apply(model@y, 1, function(x) which(x == 1))
        pp <- fitted.values(model)
        pp <- pp[cbind(1:nrow(pp), y)]  # gets the y[i]'th entry from the i'th
                                        # row of pp
        ans <- log(pp)
    } else if (inherits(model, "polr")) {
        y <- as.numeric(model.response(model$model))
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
    } else {
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
    }

    return(ans)
}

##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title
##' @aliases vuong clarke
##' @param model1 
##' @param model2 
##' @param outcome1 
##' @param outcome2 
##' @param level 
##' @return
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
vuong <- function(model1, model2, outcome1 = NULL, outcome2 = NULL,
                  level = 0.05)
{
    ## be sure to note no checking for equal order, etc etc

    ## clarke and vuong should have combined help page
    
    n <- nobs(model1)
    if (nobs(model2) != n)
        stop("model1 and model2 must have same number of observations")
    
    loglik1 <- indivLogLiks(model1, outcome1)
    loglik2 <- indivLogLiks(model2, outcome2)
    p1 <- nparams(model1)
    p2 <- nparams(model2)

    correction <- (p1 - p2) * (log(n) / 2)
    num <- sum(loglik1) - sum(loglik2) - correction
    denom <- sd(loglik1 - loglik2) * sqrt(n / (n-1))  # correcting for R's use
                                                      # of n-1 denominator in sd
                                                      # calculation
    stat <- num / (sqrt(n) * denom)

    ans <- list(stat = stat,
                test = "vuong",
                level = level,
                loglik1 = loglik1,
                loglik2 = loglik2,
                nparams = c(p1, p2),
                nobs = n)
    class(ans) <- "nonnest.test"

    return(ans)
}

clarke <- function(model1, model2, outcome1 = NULL, outcome2 = NULL,
                   level = 0.05)
{
    n <- nobs(model1)
    if (nobs(model2) != n)
        stop("model1 and model2 must have same number of observations")
    
    loglik1 <- indivLogLiks(model1, outcome1)
    loglik2 <- indivLogLiks(model2, outcome2)
    p1 <- nparams(model1)
    p2 <- nparams(model2)

    correction <- (p1 - p2) * (log(n) / (2*n))
    stat <- sum(loglik1 - loglik2 > correction)

    ans <- list(stat = stat,
                test = "clarke",
                level = level,
                loglik1 = loglik1,
                loglik2 = loglik2,
                nparams = c(p1, p2),
                nobs = n)
    class(ans) <- "nonnest.test"

    return(ans)
}
