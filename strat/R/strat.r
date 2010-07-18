##' A package for estimating strategic statistical models.
##' 
##' @name strat-package
##' @docType package
##' @references
##' Curtis S. Signorino.  2003.  \dQuote{Structure and Uncertainty
##' in Discrete Choice Models.}  \emph{Political Analysis} 11:316--344.
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
##' @S3method summary strat
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

##' Print output from \code{summary.strat}
##'
##' Prints the standard regression results table from a fitted strategic model,
##' along with the log-likelihood, AIC, and number of observations.
##' @title Print strategic model summary
##' @method print summary.strat
##' @S3method print summary.strat
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
##' @S3method coef strat
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
##' @S3method vcov strat
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

##' Makes a LaTeX table of strategic model results.
##'
##' \code{latexTable} prints LaTeX code for the presentation of results from a
##' strategic model.  Each row contains one regressor, and each column contains
##' one of the utility (or variance term) equations in the model.  For example,
##' in a model fit with \code{\link{strat3}}, the four columns will be u11, u13,
##' u14, and u24 respectively.  Each cell contains the estimated parameter, atop
##' its standard error in parentheses.  Cells that are empty because a regressor
##' is not included in a particular equation are filled with the string
##' specified in the option \code{blankfill}.  Signorino and Tarar (2003,
##' p. 593) contains a table of this variety.
##' 
##' The \code{digits} option does not yet work seamlessly; you may have to
##' resort to trial and error.
##' @title LaTeX table for strategic models
##' @param x a fitted model of class \code{strat}
##' @param digits number of digits to print
##' @param blankfill text to fill empty cells (those where a certain variable
##' did not enter the given equation)
##' @param math.style.negative whether negative signs should be "math style" or
##' plain hyphens.  Defaults to \code{TRUE}
##' @param file file to save the output in.  Defaults to \code{""}, which prints
##' the table to the R console
##' @param floatplace where to place the table float; e.g., for
##' \code{\\begin\{table\}[htp]}, use \code{floatplace = "htp"}
##' @param rowsep amount of space (in points) to put between rows
##' @return \code{x}, invisibly
##' @export
##' @references Curtis S. Signorino and Ahmer Tarar.  2006.  \dQuote{A Unified
##' Theory and Test of Extended Immediate Deterrence.}  \emph{American Journal
##' of Political Science} 50(3):586--605.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
latexTable <- function(x, digits = max(3, getOption("digits") - 2),
                       blankfill = "", math.style.negative = TRUE,
                       file = "", floatplace = "htbp", rowsep = 2)
{
    if (any(grepl(":", levels(x$y)))) {
        warning("names of outcomes (", paste(levels(x$y), collapse = ", "),
                ") contain colons; undefined behavior may result")
    }
    
    lcat <- function(...) cat(..., file = file, sep = "", append = TRUE)

    n <- names(coef(x))
    cf <- summary(x)$coefficients[, 1:2]
    cf <- format(cf, digits = digits, trim = TRUE)
    if (math.style.negative)
        cf <- sub("-", "$-$", cf)
    
    eqNames <- unique(sapply(strsplit(n, ":"), '[', 1))
    varNames <- sapply(sapply(strsplit(n, ":"), '[', -1), paste, collapse = ":")
    varNames <- unique(varNames)

    lcat("\n%% latex table generated in R ", as.character(getRversion()),
         " by strat package\n")
    lcat("%% ", date(), "\n")
    lcat("%% remember to include \\usepackage{multirow} in your preamble\n\n")
    lcat("\\begin{table}[", floatplace, "]\n")
    lcat("\\begin{center}\n")
    lcat("\\begin{tabular}{",
         paste(c("l", rep("c", length(eqNames))), collapse = ""), "}\n")
    lcat("\\hline\n")
    lcat(paste(c("", eqNames), collapse = " & "), " \\\\\n")
    lcat("\\hline\n")

    for (i in varNames) {
        lcat("\\multirow{2}{*}{", i, "} & ")
        for (J in 1:length(eqNames)) {
            j <- eqNames[J]
            ji <- paste(j, i, sep = ":")
            if (ji %in% n) {
                lcat(cf[ji, 1])
            } else {
                lcat("\\multirow{2}{*}{", blankfill, "}")
            }
            if (J < length(eqNames))
                lcat(" & ")
        }
        lcat(" \\\\\n & ")
        for (J in 1:length(eqNames)) {
            j <- eqNames[J]
            ji <- paste(j, i, sep = ":")
            if (ji %in% n)
                lcat("(", cf[ji, 2], ")")
            if (J < length(eqNames))
                lcat(" & ")
        }
        lcat(" \\\\[", rowsep, "pt]\n")

    }

    lcat("\\hline\n")
    lcat("\\end{tabular}\n")
    lcat("\\end{center}\n")
    lcat("\\end{table}\n")

    invisible(x)
}

finiteProbs <- function(x)
{
    x <- replace(x, x < .Machine$double.eps, .Machine$double.eps)
    x <- replace(x, x > 1 - .Machine$double.neg.eps,
                 1 - .Machine$double.neg.eps)
    return(x)
}
