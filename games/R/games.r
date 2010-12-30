##' A package for estimating strategic statistical models.
##' 
##' @name games-package
##' @docType package
##' @references
##' Curtis S. Signorino.  2003.  \dQuote{Structure and Uncertainty
##' in Discrete Choice Models.}  \emph{Political Analysis} 11:316--344.
NULL

##' The default method for printing a \code{game} object.
##'
##' Prints the call and coefficients of a fitted strategic model.
##' @title Print a strategic model object
##' @param x a fitted model of class \code{game}
##' @param ... other arguments, currently ignored
##' @return \code{x}, invisibly
##' @method print game
##' @S3method print game
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
print.game <- function(x, ...)
{
    oldx <- x
    cat("\nA fitted strategic model\n\nCALL:\n\n")
    print(x$call)
    cat("\nCOEFFICIENTS:\n")

    for (i in seq_along(x$equations)) {
        eq <- x$equations[i]
        hc <- attr(x$equations, "hasColon")[i]
        cf <- grep(eq, names(x$coefficients), fixed = TRUE)

        cat("\n  ", prefixToString(eq), "\n", sep = "")
        if (length(cf) > 0) {
            isFixed <- all(x$fixed[cf])
            cf <- x$coefficients[cf]

            if (hc) {
                ## this strips out the equation prefix in each term; e.g.,
                ## "u1(war):x1" becomes "x1"
                names(cf) <- sapply(strsplit(names(cf),
                                             paste(eq, ":", sep = ""),
                                             fixed = TRUE), "[", -1)
            } else {  ## i.e., the term is estimated itself, without regressors
                names(cf) <- if (isFixed) "fixed to" else "estimated as"
            }
            
            names(cf) <- paste("     ", names(cf), sep = "")
            cf <- data.frame(as.matrix(cf))
            names(cf) <- " "
            print(cf)
        } else {
            ## this is for cases when there is a utility equation with no terms
            ## estimated
            cat("\n     fixed to 0\n")
        }
    }

    if (x$convergence$code) {
        cat("\nWarning: Model fitting did not converge\nMessage:",
            x$convergence$message)
    }
    
    cat("\n")
    invisible(oldx)
}

##' The default method for summarizing a \code{game} object.
##'
##' Forms the standard regression results table from a fitted strategic model.
##' Normally used interactively, in conjunction with
##' \code{\link{print.summary.game}}.
##' @title Summarize a strategic model object
##' @method summary game
##' @S3method summary game
##' @param object a fitted model of class \code{game}
##' @param useboot logical: use bootstrap estimates (if present) to construct
##' standard error estimates?
##' @param ... other arguments, currently ignored
##' @return an object of class \code{summary.game}, containing the coefficient
##' matrix and other information needed for printing
##' @seealso \code{\link{print.summary.game}}
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
summary.game <- function(object, useboot = TRUE, ...)
{
    useboot <- useboot && !is.null(object$boot.matrix)
    if (useboot)
        object$vcov <- var(object$boot.matrix)
    cf <- object$coefficients[!object$fixed]
    se <- sqrt(diag(object$vcov[!object$fixed, !object$fixed, drop = FALSE]))
    zval <- cf / se
    pval <- 2 * pnorm(-abs(zval))

    ans <- list()
    ans$coefficients <- cbind(cf, se, zval, pval)
    colnames(ans$coefficients) <- c("Estimate", "Std. Error", "z value",
                                    "Pr(>|z|)")
    ans$call <- object$call
    ans$log.likelihood <- sum(object$log.likelihood)
    ans$nobs <- nrow(object$model)
    ans$fixed.terms <- object$coefficients[object$fixed]
    ans$convergence <- object$convergence
    ans$useboot <- useboot
    class(ans) <- "summary.game"

    return(ans)
}

##' Print output from \code{summary.game}
##'
##' Prints the standard regression results table from a fitted strategic model,
##' along with the log-likelihood, AIC, and number of observations.
##' @title Print strategic model summary
##' @method print summary.game
##' @S3method print summary.game
##' @param x an object of class \code{summary.game}, typically produced by
##'running \code{summary} on a fitted model of class \code{game}
##' @param ... other arguments, currently ignored
##' @return \code{x}, invisibly.
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
print.summary.game <- function(x, ...)
{
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients)
    if (x$useboot) {
        cat("\nStandard errors estimated from bootstrap results\n")
    } else {
        cat("\nStandard errors estimated from inverse Hessian\n")
    }
    if (length(x$fixed.terms)) {
        cat("\nFixed terms:\n")
        print(x$fixed)
    }
    cat("\nLog-likelihood:", x$log.likelihood)
    cat("\nAIC:", AIC(x))
    cat("\nNo. observations:", x$nobs, "\n\n")
    if (x$convergence$code) {
        cat("\nWarning: Model fitting did not converge\nMessage:",
            x$convergence$message, "\n")
    }
    invisible(x)
}

coef.game <- function(object, ...)
{
    object$coefficients
}

vcov.game <- function(object, ...)
{
    object$vcov
}

logLik.game <- function(object, ...)
{
    ans <- sum(object$log.likelihood)
    attr(ans, "df") <- length(object$coefficients) - sum(object$fixed)
    attr(ans, "nobs") <- nrow(object$model)
    class(ans) <- "logLik"
    return(ans)
}

logLik.summary.game <- function(object, ...)
{
    ans <- object$log.likelihood
    attr(ans, "df") <- nrow(object$coefficients)
    attr(ans, "nobs") <- object$nobs
    class(ans) <- "logLik"
    return(ans)
}

##' Makes predicted probabilities from a strategic model .
##'
##' This method uses a fitted strategic model to make predictions for a new
##' set of data.  Useful for cross-validating or for graphical analysis.
##'
##' For many uses, such as analyzing the marginal effect of a particular
##' independent variable, the function \code{\link{predProbs}} will be more
##' convenient.
##' @title Predicted probabilities for strategic models
##' @aliases predict.game predict.egame12 predict.egame122 predict.ultimatum
##' @usage
##' \method{predict}{game}(object, ...)
##'
##' \method{predict}{egame12}(object, newdata, probs=c("outcome", "action"), ...)
##' \method{predict}{egame122}(object, newdata, probs=c("outcome", "action"), ...)
##' \method{predict}{ultimatum}(object, newdata, ...)
##' @param object a fitted model of class \code{game}.
##' @param newdata data frame of values to make the predicted probabilities for.
##' If this is left empty, the original dataset is used.
##' @param probs whether to provide probabilities for outcomes (e.g., L, RL, or
##' RR in \code{egame12}) or for actions (e.g., whether 2 moves L or R given
##' that 1 moved R).
##' @param ... other arguments, currently ignored.
##' @return A data frame of predicted probabilities.
##' @method predict game
##' @S3method predict game
##' @export
##' @seealso \code{\link{predProbs}} provides a more full-featured and
##' user-friendly wrapper, including plots and confidence bands.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
predict.game <- function(object, ...)
{
    NextMethod("predict", object, ...)
}

##
## Takes an equation prefix from a strategic model (e.g., "u1(war)") and
## translates it into plain English (e.g., "Player 1's utility for war").
##
prefixToString <- function(x)
{
    first <- substr(x, 1, 1)
    
    if (first == "u") {
        n <- 3
        while (substr(x, n, n) != "(") n <- n + 1

        player <- substr(x, 2, n - 1)
        outcome <- substr(x, n + 1, nchar(x) - 1)

        x <- paste("Player ", player, "'s utility for ", outcome, ":", sep = "")
    } else if (first == "l") {
        player <- substr(x, nchar(x) - 1, nchar(x) - 1)
        if (player == "1" || player == "2") {
            x <- paste("Logged scale parameter for player ", player, ":", sep
                       = "")
        } else {
            x <- "Logged scale parameter:"
        }
    } else if (first == "R") {
        x <- paste("Player ", substr(x, 2, nchar(x)), "'s reservation value:",
                   sep = "")
    }

    return(x)
}
