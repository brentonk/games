##' @include games.r
##' @include helpers.r
NULL

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
    ans <- structure(cbind(ans, profData), probcols = probcols, lowcols =
                     lowcols, highcols = highcols, xcol = xcol, class =
                     c("predProbs", "data.frame"))

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
