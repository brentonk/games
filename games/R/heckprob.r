##' @include games.r
##' @include helpers.r
NULL

require(pbivnorm)  # remove this, and change DESCRIPTION, when building final
                   # version of the package (since there may need to be quick
                   # bug fixes to 1.0 before heckprob can be introduced)

##
## INPUT:
##   x: vector of values on the real line
##
## RETURN:
##   corresponding values between -1 and 1
## 
invrhobit <- function(x)
{
    num <- finitize(exp(x) - 1)
    denom <- finitize(exp(x) + 1)
    ans <- num / denom
    return(ans)
}

##
## finiteRho:
##   prevent numerical issues in likelihood evaluation for the heckprob model
##   when the value of rho is numerically equal to -1 or 1
##
## INPUT:
##   x: vector of values in [-1, 1]
##
## RETURN:
##   corresponding values numerically within (-1, 1)
##
finiteRho <- function(x)
{
    if (any(x < -1) || any(x > 1))
        stop("all values of x must be within [-1, 1]")

    lowerBound <- -(1 - .Machine$double.neg.eps)
    upperBound <- 1 - .Machine$double.neg.eps

    x <- replace(x, x < lowerBound, lowerBound)
    x <- replace(x, x > upperBound, upperBound)
    return(x)
}

## TODO: split this up into the usual sequence of functions so that predProbs
## and its friends will still work?  or can that be accomplished just be writing
## the predict function correctly?
logLikHeckprob <- function(b, y, regr, ...)
{
    names(regr) <- c("Xs", "Xo")
    bs <- b[1:ncol(regr$Xs)]
    bo <- b[(ncol(regr$Xs)+1):(length(b)-1)]
    rho <- finiteRho(invrhobit(b[length(b)]))

    fits <- as.numeric(tcrossprod(Xs, t(bs)))
    fito <- as.numeric(tcrossprod(Xo, t(bo)))

    ## for certain values of rho very close to -1 or 1, genz's algorithm returns
    ## negative probabilities (around -1e-18 or so), so the "pmax" calls ensure
    ## that these are converted to 0, and then finiteProbs ensures that these
    ## are numerically nonzero
    p1 <- finiteProbs(1 - pnorm(fits))
    p2 <- finiteProbs(pmax(pbivnorm(fits, -fito, -rho), 0))
    p3 <- finiteProbs(pmax(pbivnorm(fits, fito, rho), 0))

    ## construct the log-likelihood
    ans <- ifelse(y == 1, log(p1), ifelse(y == 2, log(p2), log(p3)))
    return(ans)
}
