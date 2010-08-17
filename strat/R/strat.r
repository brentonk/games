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
    oldx <- x
    cat("\nA fitted strategic model\n\nCALL:\n\n")
    print(x$call)
    cat("\nCOEFFICIENTS:\n")

    names(x$coefficients)[x$fixed] <-
        paste(names(x$coefficients)[x$fixed], "fixed to", sep = ":")

    for (eq in x$equations) {
        cf <- grep(eq, names(x$coefficients), fixed = TRUE)
        cat("\n  ", prefixToString(eq), "\n", sep = "")
        if (length(cf) > 0) {
            cf <- x$coefficients[cf]

            ## this strips out the equation prefix in each term; e.g.,
            ## "u1(war):x1" becomes "x1"
            names(cf) <- sapply(strsplit(names(cf), paste(eq, ":", sep = ""),
                                         fixed = TRUE), "[", -1)

            ## this is a hack for the estimated variance terms in `stratult`;
            ## the strsplit code just above returns character(0) in these cases
            ## since the names of these ("log(s1)", "log(s2)") don't contain a
            ## colon
            names(cf)[names(cf) == "character(0)"] <- "estimated as"
            
            names(cf) <- paste("     ", names(cf), sep = "")
            cf <- data.frame(as.matrix(cf))
            names(cf) <- " "
            print(cf)
        } else {
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

##' The default method for summarizing a \code{strat} object.
##'
##' Forms the standard regression results table from a fitted strategic model.
##' Normally used interactively, in conjunction with
##' \code{\link{print.summary.strat}}.
##' @title Summarize a strategic model object
##' @method summary strat
##' @S3method summary strat
##' @param object a fitted model of class \code{strat}
##' @param useboot logical: use bootstrap estimates (if present) to construct
##' standard error estimates?
##' @param ... other arguments, currently ignored
##' @return an object of class \code{summary.strat}, containing the coefficient
##' matrix and other information needed for printing
##' @seealso \code{\link{print.summary.strat}}
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
summary.strat <- function(object, useboot = TRUE, ...)
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
    attr(ans, "df") <- length(object$coefficients) - sum(object$fixed)
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
    ## TODO: use `equations` attribute of strat objects instead
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

## Ensures that estimated probabilities aren't numerically equal to 1 or 0, in
## order to ensure no -Infs or 0s in log-likelihoods
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
## Takes an equation prefix from a strategic model (e.g., "u1(war)") and
## translates it into plain English (e.g., "Player 1's utility for war")
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
    } else if (first == "v" || first == "s") {
        x <- paste("Standard dev. term ", substr(x, 2, nchar(x)), " (logged):",
                   sep = "")
    } else if (first == "R") {
        x <- paste("Player ", substr(x, 2, nchar(x)), "'s reservation value:",
                   sep = "")
    }

    return(x)
}

##' <description>
##'
##' <details>
##' @title Mode of a vector
##' @param x A vector of class \code{numeric}, \code{logical}, \code{character},
##' \code{factor}, or \code{ordered}
##' @return The mode of \code{x}, in a vector of the same class
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
Mode <- function(x)
{
    makechar <- is.character(x)
    if (makechar)
        x <- as.factor(x)

    makenum <- is.numeric(x)
    if (makenum)
        x <- as.factor(x)

    makelog <- is.logical(x)

    xtab <- table(x)
    ans <- which(xtab == max(xtab))
    if (length(ans) > 1) {
        warning("multiple modes: ", paste(names(ans), collapse = ", "),
                "; ", names(ans)[1], " will be used")
        ans <- ans[1]
    }
    ans <- x[x == names(ans)][1]
    
    if (makenum)
        ans <- as.numeric(as.character(ans))

    if (makelog)
        ans <- as.logical(ans)

    if (makechar)
        ans <- as.character(ans)

    return(ans)
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
                warning(names(x)[i], " has equal number of 0s and 1s; ",
                        "defaulting to 1")
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

CIfromBoot <- function(x, newdata, ci = .95, report = TRUE, ...)
{
    n <- nrow(x$boot.matrix)
    forDims <- predict(x, newdata = newdata)
    ans <- vector("list", ncol(forDims))
    for (i in seq_along(ans))
        ans[[i]] <- matrix(nrow = n, ncol = nrow(forDims))

    ## Gets the predicted values for each observation using each bootstrapped
    ## coefficient vector
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

    q <- (1 - ci) / 2
    lows <- lapply(ans, function(x) apply(x, 2, quantile, probs = q))
    lows <- do.call(cbind, lows)
    highs <- lapply(ans, function(x) apply(x, 2, quantile, probs = 1 - q))
    highs <- do.call(cbind, highs)

    colnames(lows) <- paste(colnames(forDims), "low", sep = ":")
    colnames(highs) <- paste(colnames(forDims), "high", sep = ":")

    return(list(lows = lows, highs = highs))
}

predProbs <- function(x, model, xlim = c(min(x), max(x)), n = 100, ci = .95,
                      makePlots = TRUE, report = TRUE, ...)
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
        warning("Bootstrap values unavailable; using normal resampling to ",
                "generate confidence intervals")
        vcv <- vcov(model)
        vcv[model$fixed, ] <- vcv[, model$fixed] <- 0L
        model$boot.matrix <- mvrnorm(1000, mu = coef(model), Sigma = vcv)
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
    class(ans) <- c("stratpp", "data.frame")

    invisible(ans)
}

##
## The next two functions are designed to mimic the behavior of the "plot"
## method for "gam" objects.  See the source of "plot.gam" and
## "plot.preplot.gam" for more.  The "gam" package is licensed under the GPL,
## and some of the code below closely matches it.
## 

plot.stratpp <- function(x, which = NULL, ask = FALSE, ...)
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
    for (i in 1:ncol(probs)) {
        preplotObj[[i]] <- list()
        preplotObj[[i]]$x <- xvar
        preplotObj[[i]]$y <- probs[, i]
        preplotObj[[i]]$low <- lows[, i]
        preplotObj[[i]]$high <- highs[, i]
        preplotObj[[i]]$xlab <- names(x)[attr(x, "xcol")]
        preplotObj[[i]]$ylab <- probnames[i]
        class(preplotObj[[i]]) <- "preplot.stratpp"
    }

    if (is.null(which) && ask) {
        tmenu <- c(paste("plot:", probnames), "plot all terms")
        pick <- 1
        while (pick > 0 && pick <= length(tmenu)) {
            pick <- menu(tmenu, title =
                         "Make a plot selection (or 0 to exit):\n")
            if (pick > 0 && pick < length(tmenu)) {
                plot.preplot.stratpp(preplotObj[[pick]], ...)
            } else if (pick == length(tmenu)) {
                plot.preplot.stratpp(preplotObj, ...)
            }
        }
    } else if (!is.null(which)) {
        plot.preplot.stratpp(preplotObj[[which]], ...)
    } else {
        plot.preplot.stratpp(preplotObj, ...)
    }

    invisible(x)
}

plot.preplot.stratpp <- function(x, xlab = x$xlab, ylab = x$ylab,
                                 ylim = c(min(x$y, x$low), max(x$y, x$high)),
                                 type = "l", lty.ci = 2, ...)
{
    listof <- inherits(x[[1]], "preplot.stratpp")
    
    if (listof) {
        for (i in seq_along(x))
            plot.preplot.stratpp(x[[i]], ...)
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
