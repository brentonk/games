##' A package for estimating strategic statistical models.
##' 
##' @name games-package
##' @docType package
##' @references
##' Curtis S. Signorino.  2003.  \dQuote{Structure and Uncertainty
##' in Discrete Choice Models.}  \emph{Political Analysis} 11:316--344.
NULL

##' Dataset of militarized international disputes between 1816 and 1899.
##'
##' The dataset is taken from the Correlates of War project.  The unit of
##' observation is the dyad-year, and the variables are: \describe{
##' \item{\code{ccode1}}{Initiator's COW country code}
##' \item{\code{ccode2}}{Respondent's COW country code}
##' \item{\code{year}}{Year of dispute}
##' \item{\code{cap_1}}{Initiator's military capabilities (as percent of total
##' system capabilities)}
##' \item{\code{cap_2}}{Respondent's military capabilities (as percent of total
##' system capabilities)}
##' \item{\code{balanc}}{Balance of dyadic capabilities possessed by the
##' initiator (i.e., \code{cap_1 / (cap_1 + cap_2)})}
##' \item{\code{s_wt_re1}}{Dyadic S-score (see Signorino and Ritter 1998),
##' weighted by initiator's region}
##' \item{\code{s_wt_re2}}{Dyadic S-score, weighted by respondent's region}
##' \item{\code{dem1}}{Initiator's Polity score}
##' \item{\code{dem2}}{Respondent's Polity score}
##' \item{\code{distance}}{Distance (in miles) between initiator and respondent}
##' \item{\code{peaceyrs}}{Years since last dispute in this dyad}
##' \item{\code{midnum}}{Dispute's number in the MID data set}
##' \item{\code{revis1}}{Whether the initiator had \dQuote{revisionist} aims}
##' \item{\code{revis2}}{Whether the respondent had \dQuote{revisionist} aims}
##' \item{\code{sq}}{Indicator for status quo outcome}
##' \item{\code{capit}}{Indicator for capitulation outcome}
##' \item{\code{war}}{Indicator for war outcome}
##' \item{\code{esc}}{Indicator for escalation (i.e., either capitulation or war
##' occurs)}
##' \item{\code{regime1}}{Initiator's regime type (calculated from \code{dem1})}
##' \item{\code{regime2}}{Respondent's regime type (calculated from \code{dem2})}
##' }
##' @name war1800
##' @usage data(war1800)
##' @title 19th-century international disputes
##' @docType data
##' @references Daniel M. Jones, Stuart A. Bremer and J. David Singer.  1996.
##' \dQuote{Militarized Interstate Disputes, 1816-1992: Rationale, Coding Rules,
##' and Empirical Patterns.} \emph{Conflict Management and Peace Science}
##' 15(2): 163--213.
##' @seealso \code{\link{egame12}}
##' @keywords data
NULL

##' Simulated data for illustrating \code{\link{egame122}}.
##'
##' The variables are: \describe{
##' \item{\code{f1}, \code{f2}}{Factors with levels \dQuote{a}, \dQuote{b},
##' \dQuote{c}}
##' \item{\code{x1}--\code{x5}}{Numeric variables entering Player 1's utilities}
##' \item{\code{z1}--\code{z3}}{Numeric variables entering Player 2's utilities}
##' \item{\code{a1}}{Indicator for Player 1's move (L or R)}
##' \item{\code{a2}}{Indicator for Player 2's move (L or R)}
##' \item{\code{y}}{Factor containing outcome}
##' }
##' @name data_122
##' @usage data(data_122)
##' @title Simulated egame122 data
##' @docType data
##' @seealso \code{\link{egame122}}
##' @keywords data
NULL

##' Simulated data for illustrating \code{\link{ultimatum}}.
##'
##' The variables are: \describe{
##' \item{\code{offer}}{The offer made by Player 1}
##' \item{\code{accept}}{Whether Player 2 accepted the offer (0 for rejection, 1
##' for acceptance)}
##' \item{\code{w1}, \code{w2}}{Variables entering both players' reservation values}
##' \item{\code{x1}--\code{x4}}{Variables entering Player 1's reservation value}
##' \item{\code{z1}--\code{z4}}{Variables entering Player 2's reservation value}
##' }
##' The maximum offer size is 15.
##' @name data_ult
##' @usage data(data_ult)
##' @title Simulated ultimatum data
##' @docType data
##' @seealso \code{\link{ultimatum}}
##' @keywords data
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
## Converts a character vector for use in LaTeX by inserting escape sequences
## where appropriate.  Not comprehensive, but should catch most common
## problems.
## 
latexEsc <- function(x)
{
    x <- gsub("{", "\\{", x, fixed = TRUE)
    x <- gsub("}", "\\}", x, fixed = TRUE)
    x <- gsub("_", "\\_", x, fixed = TRUE)
    x <- gsub("#", "\\#", x, fixed = TRUE)
    x <- gsub("$", "\\$", x, fixed = TRUE)
    x <- gsub("%", "\\%", x, fixed = TRUE)
    x <- gsub("^", "\\^", x, fixed = TRUE)
    x <- gsub("&", "\\&", x, fixed = TRUE)
    x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
    return(x)
}

##' Makes a LaTeX table of strategic model results.
##'
##' \code{latexTable} prints LaTeX code for the presentation of results from a
##' strategic model.  Each row contains one regressor, and each column contains
##' one of the utility (or variance term) equations in the model.  For example,
##' in a model fit with \code{\link{egame12}}, the four columns will be u11,
##' u13, u14, and u24 respectively.  Each cell contains the estimated parameter,
##' atop its standard error in parentheses.  Cells that are empty because a
##' regressor is not included in a particular equation are filled with the
##' string specified in the option \code{blankfill}.  Signorino and Tarar (2006,
##' p. 593) contains a table of this variety.
##'
##' The table generated depends on the \pkg{multirow} package for LaTeX, so
##' make sure to include \code{\\usepackage{multirow}} in the preamble of your
##' document.
##' 
##' The \code{digits} option does not yet work seamlessly; you may have to
##' resort to trial and error.
##' @title LaTeX table for strategic models
##' @param x a fitted model of class \code{game}.
##' @param digits number of digits to print.
##' @param scientific logical or integer value to control use of scientific
##' notation.  See \code{\link{format}}.
##' @param blankfill text to fill empty cells (those where a certain variable
##' did not enter the given equation).
##' @param math.style.negative whether negative signs should be "math style" or
##' plain hyphens.  Defaults to \code{TRUE}.
##' @param file file to save the output in.  Defaults to \code{""}, which prints
##' the table to the R console.
##' @param floatplace where to place the table float; e.g., for
##' \code{\\begin\{table\}[htp]}, use \code{floatplace = "htp"}.
##' @param rowsep amount of space (in points) to put between rows.
##' @param useboot logical: use bootstrap estimates (if available) to calculate
##' standard errors?
##' @return \code{x}, invisibly.
##' @export
##' @references Curtis S. Signorino and Ahmer Tarar.  2006.  \dQuote{A Unified
##' Theory and Test of Extended Immediate Deterrence.}  \emph{American Journal
##' of Political Science} 50(3):586--605.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @examples
##' data(war1800)
##' f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | balanc + revis1 | balanc
##' m1 <- egame12(f1, data = war1800)
##'
##' latexTable(m1)
##' latexTable(m1, digits = 8)
##' latexTable(m1, blankfill = "--")  ## dashes in blank cells
##' 
##' \dontrun{
##' latexTable(m1, file = "my_table.tex")  ## write to file}
latexTable <- function(x, digits = max(3, getOption("digits") - 2), scientific =
                       NA, blankfill = "", math.style.negative = TRUE, file =
                       "", floatplace = "htbp", rowsep = 2, useboot = TRUE)
{
    lcat <- function(...) cat(..., file = file, sep = "", append = TRUE)

    n <- names(coef(x))
    cf <- summary(x, useboot = useboot)$coefficients[, 1:2]
    cf <- rbind(cf, c(sum(x$log.likelihood), 0))
    cf <- format(cf, digits = digits, trim = TRUE, scientific = scientific)
    if (math.style.negative)
        cf <- gsub("-", "$-$", cf)

    eqNames <- x$equations[attr(x$equations, "hasColon")]
    varNames <- sapply(sapply(strsplit(n, ":"), '[', -1), paste, collapse = ":")
    varNames <- unique(varNames[nchar(varNames) > 0])
    otherNames <- x$equations[!attr(x$equations, "hasColon") &
                              x$equations %in% n[!x$fixed]]

    lcat("\n%% latex table generated in R ", as.character(getRversion()),
         " by games package\n")
    lcat("%% ", date(), "\n")
    lcat("%% remember to include \\usepackage{multirow} in your preamble\n\n")
    lcat("\\begin{table}[", floatplace, "]\n")
    lcat("\\begin{center}\n")
    lcat("\\begin{tabular}{",
         paste(c("l", rep("c", length(eqNames))), collapse = ""), "}\n")
    lcat("\\hline\n")
    lcat(paste(c("", latexEsc(eqNames)), collapse = " & "), " \\\\\n")
    lcat("\\hline\n")

    for (i in varNames) {
        lcat("\\multirow{2}{*}{", latexEsc(i), "} & ")
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

    if (length(otherNames) > 0) {
        lcat("\\hline\n")
        for (i in otherNames) {
            lcat("\\multirow{2}{*}{", i, "} & ", cf[i, 1], " \\\\\n & (",
                 cf[i, 2], ") \\\\[", rowsep, "pt]\n")
        }
    }

    lcat("\\hline \\hline\n")
    lcat("Log-likelihood & ", cf[nrow(cf), 1], " \\\\\n $N$ & ",
         nrow(x$model), "\\\\\n")
    lcat("\\hline\n")
    
    lcat("\\end{tabular}\n")
    lcat("\\end{center}\n")
    lcat("\\end{table}\n")

    invisible(x)
}

##
## Ensures that estimated probabilities aren't numerically equal to 1 or 0, in
## order to ensure no -Infs or 0s in log-likelihoods.
##
finiteProbs <- function(x)
{
    x <- replace(x, x < .Machine$double.eps, .Machine$double.eps)
    x <- replace(x, x > 1 - .Machine$double.neg.eps,
                 1 - .Machine$double.neg.eps)
    return(x)
}

finitize <- function(x)
{
    x <- ifelse(is.finite(x), x, sign(x) * .Machine$double.xmax)
    return(x)
}

##
## Used to ensure that the "formulas" argument of each fitting function contains
## a valid type of object and coerces it to "Formula" class.  Returns an error
## if the function isn't a formula, Formula, or list of formulas.
##
checkFormulas <- function(f, argname = "formulas")
{
    if (inherits(f, "list")) {
        f <- do.call(as.Formula, f)
    } else if (inherits(f, "formula")) {
        f <- as.Formula(f)
    } else {
        stop(argname, " must be a list of formulas or a formula")
    }

    return(f)
}

##
## Takes a list of vectors and finds their intersection
##
intersectAll <- function(...)
{
    x <- list(...)
    ans <- x[[1]]
    for (i in 1:length(x)) ans <- intersect(ans, x[[i]])

    return(ans)
}

##
## Makes the names of the variables for egame12 and egame122 models.
## 
makeVarNames <- function(varNames, prefixes, link, sdterms)
{
    vname <- if (link == "logit") "log(lambda" else "log(sigma"
    if (sdterms == 1) {
        prefixes <- c(prefixes, paste(vname, ")", sep = ""))
    } else if (sdterms == 2) {
        prefixes <- c(prefixes, paste(vname, 1:2, ")", sep = ""))
    }

    hasColon <- sapply(varNames, function(x) length(x) > 0 &&
                       !all(x == "(Intercept)"))
    names(hasColon) <- prefixes
    for (i in seq_along(varNames)) {
        if (hasColon[i]) {
            varNames[[i]] <- paste(prefixes[i], varNames[[i]], sep = ":")
        } else {
            varNames[[i]] <- prefixes[i][length(varNames[[i]])]
        }
    }
    varNames <- unlist(varNames)
    ans <- list(varNames = varNames, hasColon = hasColon)
    return(ans)
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

##
## Calculates the variance-covariance matrix for a fitted model, including a
## procedure for catching the error (and returning a matrix of NAs) in case the
## Hessian is non-invertible.
##
getGameVcov <- function(hessian, fixed)
{
    hes <- hessian[!fixed, !fixed, drop = FALSE]
    vv <- tryCatch(solve(-hes), error = function(e) e)
    if (inherits(vv, "error")) {
        warning("variance-covariance matrix could not be calculated: ",
                vv$message)
        vv <- matrix(NA, nrow(hes), nrow(hes))
    }
    ans <- hessian
    ans[] <- NA
    ans[!fixed, !fixed] <- vv
    return(ans)
}

##
## Calculates bootstrap results for a strategic model.
## 
gameBoot <- function(boot, report = TRUE, estimate, y, a = NULL, regr, fn, gr,
                      fixed, ...)
{
    bootMatrix <- matrix(NA, nrow = boot, ncol = length(estimate))
    failedBoot <- logical(boot)
    if (report) {
        cat("\nRunning bootstrap iterations...\n")
        pb <- txtProgressBar(min = 1, max = boot)
    }
    for (i in seq_len(boot)) {
        bootSamp <- sample(seq_len(length(y)), replace = TRUE)
        newy <- y[bootSamp]
        newa <- a[bootSamp]  ## for the ultimatum model
        newregr <- lapply(regr, function(x) x[bootSamp, , drop = FALSE])
        bootResults <- maxBFGS(fn = fn, grad = gr, start = estimate, fixed =
                               fixed, y = newy, acc = newa, regr = newregr, ...)
        if (bootResults$code) {
            warning("bootstrap iteration ", i,
                    "failed to converge and will be removed")
            failedBoot[i] <- TRUE
        }
        bootMatrix[i, ] <- bootResults$estimate

        if (report)
            setTxtProgressBar(pb, i)
    }
    if (report)
        cat("\n")
    bootMatrix <- bootMatrix[!failedBoot, , drop = FALSE]
    colnames(bootMatrix) <- names(estimate)
    return(bootMatrix)
}

##' Finds the modal value of a vector of any class.
##'
##' Based on code from the R wiki:
##' \url{http://rwiki.sciviews.org/doku.php?id=tips:stats-basic:modalvalue}.
##' @title Mode of a vector
##' @param x a vector (lists and arrays will be flattened).
##' @param na.rm logical: strip \code{NA} values?
##' @return The value of \code{x} that occurs most often.  If there is a tie,
##' the one that appears first (among those tied) is chosen.
##' @export
##' @author R Wiki contributors
##' @examples
##' x <- c(1, 3, 3, 4)
##' Mode(x)  ## 3
##' x.char <- letters[x]
##' Mode(x.char)  ## "c"
##' x.factor <- as.factor(x.char)
##' Mode(x.factor)  ## "c", with levels "a", "c", "d"
##' x.logical <- x > 3
##' Mode(x.logical)  ## FALSE
##'
##' ## behavior with ties
##' y <- c(3, 3, 1, 1, 2)
##' Mode(y)  ## 3
##' z <- c(2, 1, 1, 3, 3)
##' Mode(z)  ## 1
Mode <- function(x, na.rm = FALSE)
{
    x <- unlist(x);
    if (na.rm)
        x <- x[!is.na(x)]
    u <- unique(x);
    n <- length(u);
    frequencies <- rep(0, n);
    for (i in seq_len(n)) {
        if (is.na(u[i])) {
            frequencies[i] <- sum(is.na(x))
        } else {
            frequencies[i] <- sum(x == u[i], na.rm = TRUE)
        }
    }
    return(u[which.max(frequencies)])
}

##
## Takes a data frame and returns a one-row data frame with the same variables,
## containing a "typical observation" profile from x.  The default is to take
## the means of numeric (non-binary) variables, the medians of ordered
## variables, and the modes of all other types of variables.  These can be
## overridden for individual variables via arguments to "..."; e.g., to set
## variable "z" to its median in the profile, use makeProfile(x, z = median(z)).
##
## This function isn't meant to be called directly by users -- it should just be
## used within the predProbs function.
## 
makeProfile <- function(x, ...)
{
    cl <- match.call(expand.dots = FALSE)

    ans <- x[1, ]
    for (i in seq_len(ncol(x))) {
        xvar <- x[, i]
        isDummy <- all(unique(xvar) %in% c(0, 1))
        if (is.numeric(xvar) && !isDummy) {
            ans[, i] <- mean(xvar)
        } else if (is.ordered(xvar)) {
            ans[, i] <- median(xvar)
        } else if (isDummy) {
            if (sum(xvar == 0) == sum(xvar == 1)) {
                warning(names(x)[i],
                        " has equal number of 0s and 1s; defaulting to 1")
                ans[, i] <- 1
            }
        } else {
            ans[, i] <- Mode(xvar)
        }
    }

    if ("..." %in% names(cl)) {
        
        ## Evaluates the expressions fed to "...", evaluates them within the
        ## supplied data frame, and returns them as a vector.  So if the call is
        ## makeProfile(x, foo = median(foo), bar = quantile(bar, .25)), this
        ## will return a vector with named elements "foo" and "bar".  The
        ## eval(substitute()) business is to ensure that these aren't evaluated
        ## in the global environment rather than "x", which would (most likely)
        ## throw an error when "foo" and "bar" weren't found in the workspace,
        ## or (possibly) obtain the wrong values for these.
        dots <- eval(substitute(list(...)), x)
        dots <- unlist(lapply(dots, unname))

        toReplace <- match(names(dots), names(x), 0L)
        ans[toReplace] <- dots[toReplace != 0L]
    }

    ## Ensures that choices in "..." expressed as characters (e.g., foo = "a"
    ## for a factor foo with levels a, b, c) don't wind up turning factor
    ## variables into characters
    fvars <- sapply(x, inherits, what = "factor")
    for (i in seq_len(ncol(x))) {
        if (fvars[i]) {
            ans[, i] <- factor(ans[, i], levels = levels(x[, i]))
        }
    }
    
    return(ans)
}

##
## Makes bootstrap confidence intervals for predicted probabilities from a
## fitted strategic model.
## 
CIfromBoot <- function(x, newdata, ci = .95, report = TRUE, ...)
{
    n <- nrow(x$boot.matrix)
    forDims <- predict(x, newdata = newdata)
    ans <- vector("list", ncol(forDims))
    for (i in seq_along(ans))
        ans[[i]] <- matrix(nrow = n, ncol = nrow(forDims))

    ## Calculates the predicted values for each observation using each
    ## bootstrapped coefficient vector
    if (report) {
        cat("\nCalculating confidence intervals...\n")
        pb <- txtProgressBar(min = 1, max = n)
    }
    for (i in seq_len(n)) {
        x$coefficients <- x$boot.matrix[i, ]
        xpred <- predict(x, newdata = newdata)
        for (j in seq_along(ans))
            ans[[j]][i, ] <- xpred[, j]
        if (report)
            setTxtProgressBar(pb, i)
    }
    if (report)
        cat("\n")

    q <- .5 - (ci / 2)
    lows <- lapply(ans, function(x) apply(x, 2, quantile, probs = q))
    lows <- do.call(cbind, lows)
    highs <- lapply(ans, function(x) apply(x, 2, quantile, probs = 1 - q))
    highs <- do.call(cbind, highs)

    colnames(lows) <- paste(colnames(forDims), "low", sep = ":")
    colnames(highs) <- paste(colnames(forDims), "high", sep = ":")

    return(list(lows = lows, highs = highs))
}

##' Easy generation and plotting of predicted probabilities from a fitted
##' strategic model.
##'
##' \code{predProbs} provides an easy way to analyze the estimated marginal
##' effect of an independent variable on the probability of particular outcomes,
##' using the estimates returned by a strategic model.  The procedure is
##' designed so that, for a preliminary analysis, the user can simply specify
##' the fitted model and the independent variable of interest, and quickly
##' obtain plots of predicted probabilities.  However, it is flexible enough to
##' allow for finely tuned analysis as well.
##' 
##' The procedure works by varying \code{x}, the variable of interest, across
##' its observed range (or one specified by the user in \code{xlim}) while
##' holding all other independent variables in the model fixed.  The profile
##' created by default is as follows (the same defaults as in the \code{sim}
##' function in the \code{Zelig} package):
##' \itemize{
##' \item numeric, non-binary variables are fixed at their means
##' \item \code{\link{ordered}} variables are fixed at their medians
##' \item all others are fixed at their modes (see \code{\link{Mode}})}
##' However, it is possible to override these defaults for any or all
##' variables.  For example, to set a variable named \code{polity} to its lower
##' quartile, call \code{predProbs} with the argument \code{polity =
##' quantile(polity, 0.25)}.  To set a factor variable to a particular level,
##' provide the name of the level as a character string (in quotes). (Also see
##' the examples below.)
##'
##' Confidence intervals for each predicted point are generated by bootstrap.
##' If \code{model} has a non-null \code{boot.matrix} element (i.e., a bootstrap
##' was performed with the model fitting), then these results are used to
##' make the confidence intervals.  Otherwise, a parametric bootstrap sample is
##' generated by sampling from a multivariate normal distribution around the
##' parameter estimates.  In this case, a warning is issued.
##'
##' For information on plotting the predicted probabilities, see
##' \code{\link{plot.predProbs}}.  The plots are made with base graphics.  If you
##' prefer to use an alternative graphics package, all the information necessary
##' to make the plots is included in the data frame returned.
##' @title User-friendly predicted probability analysis
##' @param model a fitted model of class \code{game}.
##' @param x character string giving the name of the variable to place
##' \dQuote{on the x-axis} while all others are held constant.  Partial matches
##' are accepted.
##' @param xlim numeric, length 2: the range that \code{x} should be varied over
##' (if \code{x} is continuous).  Defaults to the observed range of \code{x}.
##' @param n integer: the number of observations to generate (if \code{x} is
##' continuous).
##' @param ci numeric: width of the confidence interval to estimate around each
##' predicted probability.  Set to \code{0} to estimate no confidence intervals.
##' @param makePlots logical: whether to automatically make the default plot
##' for the returned object.  See \code{\link{plot.predProbs}}.
##' @param report logical: whether to print a status bar while obtaining the
##' confidence intervals for the predicted probabilities.
##' @param ... used to set values for variables other than \code{x} in the
##' profile of observations.  See \dQuote{Details} and \dQuote{Examples}.
##' @return An object of class \code{predProbs}.  This is a data frame containing
##' each hypothetical observation's predicted probability, the upper and lower
##' bounds of the confidence interval, and the value of each regressor.
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com}).  Code for
##' escaping special regex characters was taken from the \code{Hmisc} package's
##' function \code{escapeRegex}, written by Charles Dupont.
##' @seealso \code{\link{predict.game}} for somewhat more flexible (but fussier)
##' generation of predicted probabilities.
##' @examples
##' data(war1800)
##' f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | regime1 | balanc + regime2
##' m1 <- egame12(f1, data = war1800, boot = 10)
##'
##' pp1 <- predProbs(m1, x = "s_wt_re1", n = 5)
##' print(pp1)  ## the hypothetical observations and their predicted probs
##' plot(pp1, which = 2)  ## see ?plot.predProbs for more plot examples
##'
##' ## changing the profile used
##' pp2 <- predProbs(m1, x = "s_wt_re1", n = 5, revis1 = 1, balanc = 0.7)
##' pp3 <- predProbs(m1, x = "s_wt_re1", n = 5, regime1 = "dem")
##' pp4 <- predProbs(m1, x = "s_wt_re1", n = 5, balanc = median(balanc))
##'
##' ## variable names (other than x) must match exactly!
##' \dontrun{
##' pp5 <- predProbs(m1, x = "s_wt_re1", bal = 0.7)  ## error will result}
##'
##' ## x can be a factor too
##' pp6 <- predProbs(m1, x = "regime1")
##' 
##' ## predProbs (despite the name) also provides predictions for the optimal
##' ## offer in ultimatum models
##' data(data_ult)
##' f2 <- offer + accept ~ x1 + x2 + x3 + x4 + w1 + w2 | z1 + z2 + z3 + z4 + w1 + w2
##' m2 <- ultimatum(f2, data = data_ult, maxOffer = 15, boot = 10)
##' pp7 <- predProbs(m2, x = "w1", n = 5)
##' print(pp7)
##'
##' op <- par(mfrow = c(2, 1))
##' plot(pp7)
##' par(op)
predProbs <- function(model, x, xlim = c(min(x), max(x)), n = 100, ci = .95,
                      makePlots = FALSE, report = TRUE, ...)
{
    xc <- charmatch(x, names(model$model))
    if (is.na(xc)) {
        stop(dQuote(x),
             " does not match any variable names in the supplied model: ",
             paste(names(model$model), collapse = ", "))
    } else if (xc == 0) {
        ## The regular expression in the next line is taken from the escapeRegex
        ## function in the Hmisc package (v3.8-2, 2010-06-22), written by
        ## Charles Dupont, licensed under GPL
        xc <- paste("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", x), sep =
                    "")
        xc <- grep(xc, names(model$model), value = TRUE)
        stop(dQuote(x),
             " matches multiple variable names in the supplied model: ",
             paste(xc, collapse = ", "))
    } else {
        x <- model$model[, xc]
    }

    if (all(unique(x) %in% c(0, 1))) {
        ## If x is binary, just use 0 and 1
        xs <- c(0, 1)
    } else if (is.numeric(x)) {
        ## If x is continuous (or treated as such in the model fitting; e.g.,
        ## count variables), make a sequence along its values controlled by the
        ## arguments xlim and n
        xs <- seq(xlim[1], xlim[2], length.out = n)
    } else if (is.factor(x)) {
        ## If x is a factor, use each of its levels
        xs <- rep(x[1], nlevels(x))
        xs[] <- levels(x)
    } else if (is.logical(x)) {
        ## If x is logical, use the logical values
        xs <- c(FALSE, TRUE)
    }

    prof <- makeProfile(model$model, ...)
    profData <- prof[rep(1, length(xs)), , drop = FALSE]
    profData[, xc] <- xs
    rownames(profData) <- seq_along(xs)

    ans <- predict(model, newdata = profData)
    if (is.list(ans))
        ans <- do.call(cbind, ans)

    if (ci > 0 && is.null(model$boot.matrix)) {
        warning("Bootstrap values unavailable; using normal resampling to generate confidence intervals")
        vcv <- vcov(model)[!model$fixed, !model$fixed, drop = FALSE]
        bm <- mvrnorm(1000, mu = coef(model)[!model$fixed], Sigma = vcv)
        bbm <- matrix(NA, nrow = 1000, ncol = length(coef(model)))
        bbm[, !model$fixed] <- bm
        for (i in seq_along(model$fixed)) {
            if (model$fixed[i])
                bbm[, i] <- coef(model)[i]
        }
        model$boot.matrix <- bbm
    }

    if (ci > 0)
        CIvals <- CIfromBoot(model, newdata = profData, ci = ci)

    ## Saving the columns that the probabilities, confidence bounds, and
    ## variable of interest are in, and making them attributes of the output, in
    ## order for the plotting function to use them
    probcols <- 1:ncol(ans)
    if (ci > 0) {
        ans <- cbind(ans, do.call(cbind, unname(CIvals)))
        lowcols <- max(probcols) + probcols
        highcols <- max(lowcols) + probcols
    } else {
        lowcols <- highcols <- numeric(0)
    }
    xcol <- max(probcols, lowcols, highcols) + xc
    ans <- cbind(ans, profData)
    attr(ans, "probcols") <- probcols
    attr(ans, "lowcols") <- lowcols
    attr(ans, "highcols") <- highcols
    attr(ans, "xcol") <- xcol
    class(ans) <- c("predProbs", "data.frame")

    if (makePlots)
        plot(ans)
    invisible(ans)
}

##
## The next two functions are designed to mimic the behavior of the "plot"
## method for "gam" objects.  See the source of "plot.gam" and
## "plot.preplot.gam" for more.  The "gam" package is licensed under the GPL,
## and some of the code below closely matches it.
## 

##' Plots predicted probabilities and associated confidence bands, using the
##' data returned from a call to \code{\link{predProbs}}.
##'
##' Most \code{predProbs} objects will be associated with multiple plots: one for
##' each outcome in the estimated model.  These are the three or four terminal
##' nodes for a \code{\link{egame12}} or \code{\link{egame122}} model
##' respectively; for an \code{\link{ultimatum}} model, these are the expected
##' offer and the probability of acceptance.  By default, \code{plot.predProbs}
##' produces plots for all of them, so only the last will be visible unless the
##' graphics device is set to have multiple figures (e.g., by setting
##' \code{par(mfrow = ...)}).  The argument \code{ask} displays a menu to select
##' among the possible plots for a given object, and \code{which} allows for
##' this to be done non-interactively.
##' @title Plot predicted probabilities
##' @param x an object of class \code{predProbs} (i.e., a data frame returned by
##' \code{\link{predProbs}}).
##' @param which optional integer specifying which plot (as numbered in the menu
##' displayed when \code{ask == TRUE}) to make.  If none is given, all available
##' plots are printed in succession.
##' @param ask logical: display interactive menu with options for which plot to
##' make?
##' @param ... further arguments to pass to the plotting function.  See
##' \code{\link{plot.default}} (when the variable on the x-axis is continuous)
##' or \code{\link{bxp}} (when it is discrete).
##' @return an object of class \code{preplot.predProbs}, invisibly.  This contains
##' the raw information used by lower-level plotting functions.
##' @method plot predProbs
##' @S3method plot predProbs
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @examples
##' data(war1800)
##' f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | regime1 | balanc + regime2
##' m1 <- egame12(f1, data = war1800, boot = 10)
##' pp1 <- predProbs(m1, x = "balanc", n = 5)  ## continuous x
##' pp2 <- predProbs(m1, x = "regime1")  ## discrete x
##'
##' ## if "ask" is FALSE and "which" isn't specified, all plots are printed
##' op <- par(mfrow = c(2, 2))
##' plot(pp1)
##' par(op)
##'
##' \dontrun{
##' plot(pp1, ask = TRUE)
##'   # displays menu:
##'   # Make a plot selection (or 0 to exit):
##'   #   1: plot: Pr(~esc)
##'   #   2: plot: Pr(esc,~war)
##'   #   3: plot: Pr(esc,war)
##'   #   4: plot all terms}
##'
##' ## to change line type for confidence bounds, use argument "lty.ci"
##' plot(pp1, which = 3, lty.ci = 3)
##'
##' ## all the standard plotting options work too
##' plot(pp1, which = 3, xlab = "Capabilities", ylab = "Probability", main = "Title")
##'
##' ## discrete x variables are plotted via R's boxplot functionality
##' plot(pp2, which = 3)
plot.predProbs <- function(x, which = NULL, ask = FALSE, ...)
{
    probs <- x[, attr(x, "probcols"), drop = FALSE]
    if (length(attr(x, "lowcols"))) {
        lows <- x[, attr(x, "lowcols"), drop = FALSE]
        highs <- x[, attr(x, "highcols"), drop = FALSE]
    } else {
        lows <- highs <- NULL
    }
    xvar <- x[, attr(x, "xcol")]
    probnames <- names(x)[attr(x, "probcols")]

    preplotObj <- vector("list", ncol(probs))
    class(preplotObj) <- "preplot.predProbs"
    for (i in 1:ncol(probs)) {
        preplotObj[[i]] <- list()
        preplotObj[[i]]$x <- xvar
        preplotObj[[i]]$y <- probs[, i]
        preplotObj[[i]]$low <- lows[, i]
        preplotObj[[i]]$high <- highs[, i]
        preplotObj[[i]]$xlab <- names(x)[attr(x, "xcol")]
        preplotObj[[i]]$ylab <- probnames[i]
        class(preplotObj[[i]]) <- "preplot.predProbs"
    }

    if (is.null(which) && ask) {
        tmenu <- c(paste("plot:", probnames), "plot all terms")
        pick <- 1
        while (pick > 0 && pick <= length(tmenu)) {
            pick <- menu(tmenu, title =
                         "Make a plot selection (or 0 to exit):\n")
            if (pick > 0 && pick < length(tmenu)) {
                plot.preplot.predProbs(preplotObj[[pick]], ...)
            } else if (pick == length(tmenu)) {
                plot.preplot.predProbs(preplotObj, ...)
            }
        }
    } else if (!is.null(which)) {
        plot.preplot.predProbs(preplotObj[[which]], ...)
    } else {
        plot.preplot.predProbs(preplotObj, ...)
    }

    invisible(preplotObj)
}

plot.preplot.predProbs <- function(x, xlab = x$xlab, ylab = x$ylab,
                                 ylim = c(min(x$y, x$low), max(x$y, x$high)),
                                 type = "l", lty.ci = 2, ...)
{
    listof <- inherits(x[[1]], "preplot.predProbs")
    cl <- match.call()
    
    if (listof) {
        for (i in seq_along(x)) {
            icall <- cl
            icall$x <- x[[i]]
            eval(icall, parent.frame())
        }
    } else if (is.factor(x$x)) {
        ## If x is a factor, then we need to make boxplots "manually"
        ## via the bxp function (see ?boxplot and ?bxp)
        boxStats <- list()
        boxStats$stats <- matrix(x$y, nrow = 5, ncol = length(x$y), byrow =
                                 TRUE)
        if (!is.null(x$low)) {
            boxStats$stats[1, ] <- x$low
            boxStats$stats[5, ] <- x$high
        }
        boxStats$n <- rep(1, length(x$x))
        boxStats$conf <- boxStats$stats[c(1, 5), ]
        boxStats$out <- numeric(0)
        boxStats$group <- numeric(0)
        boxStats$names <- as.character(x$x)

        bxp(z = boxStats, xlab = xlab, ylab = ylab, ...)
    } else {
        plot(x$x, x$y, type = type, xlab = xlab, ylab = ylab, ylim = ylim, ...)
        if (!is.null(x$low)) {
            lines(x$x, x$low, lty = lty.ci)
            lines(x$x, x$high, lty = lty.ci)
        }
    }
}

##' Solves for W in the equation \eqn{W e^W = x}{W * exp(W) = x}.
##'
##' The function is based on the code given in Barry et al. (1995).  It is used
##' to calculate fitted values for the \code{\link{ultimatum}} model.
##'
##' If negative values of \code{x} are supplied, \code{NaN}s will likely be
##' returned.
##' @title Lambert's W
##' @param x vector of values to solve for.
##' @return Solutions to Lambert's W for each value in \code{x}.
##' @export
##' @references D.A. Barry, P.J. Culligan-Hensley, and S.J. Barry.  1995.
##' \dQuote{Real Values of the W-Function.}  \emph{ACM Transactions on
##' Mathematical Software} 21(2):161--171.
##' @author Curt Signorino (\email{curt.signorino@@rochester.edu})
##' @examples
##' x <- rexp(10)
##' w <- LW(x)
##' all.equal(x, w * exp(w))
LW <- function(x)
{
    ## Note: there is a Lambert's W package, but it depends on the gsl (GNU
    ## Scientific Library) package, which is hellish for Windows users to
    ## install, hence our hand-rolled Lambert's W.
    
    eW <- function(x, W)
    {
        zn <- log(x/W)-W
        first <- zn/(1+W)
        common <- 2*(1+W)*(1+W+2*zn/3)
        first*(common-zn)/(common-2*zn)
    }
    
    lx <- log(x)
    wlt <- (x*(1+4/3*x))/(1+x*(7/3+5/6*x))
    wgt <- lx-(24*(lx*(lx+2)-3))/(lx*(7*lx+58)+127)
    W1 <- ifelse(x < 0.7385, wlt, wgt)
    
    xge20 <- x >= 20
    i <- seq_along(x)
    ixge20 <- i[xge20]
    a1 <- 1.124491989777808
    b1 <- .4225028202459761
    xg <- x[ixge20]
    h <- exp(-a1/(b1+log(xg)))
    Wp2 <- log(xg/log(xg/(log(xg))^h))  # for x>20
    W1[ixge20] <- Wp2;                   

    ## iteration for improved accuracy
    W2 <- W1*(1+eW(x,W1))
    W3 <- W2*(1+eW(x,W2))
    W4 <- W3*(1+eW(x,W3))               # enough

    return(W4)
}

##
## Individual log-likelihoods for a model of class "game", "glm", "vglm"
## (multinomial logit via VGAM or Zelig), "polr" (ordered logit/probit via MASS
## or Zelig)
##
## The "outcome" argument is for comparing a strategic model to a binary logit
## or probit model where the response is one of the possible outcomes (see the
## examples for the "vuong"/"clarke" functions)
## 
indivLogLiks <- function(model, outcome)
{
    ## need to deal with weights?

    ## change how it works for the ultimatum model?  model$y == outcome will
    ## give the wrong results.  OTOH if only outcome of interest is the offer,
    ## that should specified -- perhaps ask curt?
    
    if (inherits(model, "game") && is.null(outcome)) {
        ans <- model$log.likelihood
    } else if (inherits(model, "game")) {
        pp <- predict(model)[, outcome]
        ans <- ifelse(as.numeric(model$y) == outcome, pp, 1 - pp)
        ans <- log(ans)
    } else if (inherits(model, "glm")) {
        ans <- ifelse(model$y == 1, fitted.values(model),
                      1 - fitted.values(model))
        ans <- log(ans)
    } else if (inherits(model, "vglm") && "multinomial" %in%
               model@family@vfamily) {
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

vuong <- function(model1, model2, outcome1 = NULL, outcome2 = NULL)
{
}
