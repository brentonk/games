##' @include games.r
##' @include helpers.r
NULL

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

vuong <- function(model1, model2, outcome1 = NULL, outcome2 = NULL)
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
    ans <- num / denom

    return(ans)
}
