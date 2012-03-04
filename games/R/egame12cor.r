##' @include games.r
##' @include helpers.r
##' @include egame12.r
NULL

##
## Transformations to map between the reals and [-1, 1]
##
rhobit <- function(x)
    log((1 + x) / (1 - x))
invrhobit <- function(x)
    (finitize(exp(x)) - 1) / (finitize(exp(x)) + 1)

##
## Conditional distribution of e2 given e1 when (e1, e2) are bivariate standard
## normal correlated at rho
## 
pnormCond <- function(x, e1, rho)
{
    pnorm(x, mean = rho * e1, sd = sqrt(1 - rho^2))
}

##
## Find approximate roots via linear interpolation. 'x' and 'y' must be
## two-column matrices with x[, 2] > x[, 1] and y[, 2] > 0 > y[, 1]
##
approxRoot <- function(x, y)
{
    b <- (y[, 2] - y[, 1]) / (x[, 2] - x[, 1])
    a <- y[, 1] - b * x[, 1]
    -a / b
}

##
## Outcome probabilities for given parameters
## 
makeProbs12cor <- function(b, regr, gridsize)
{
    ## Compute utilities
    utils <- makeUtils(b, regr, nutils = 4,
                       unames = c("u11", "u13", "u14", "u24"))
    rho <- finiteProbs(invrhobit(b[length(b)]))

    ## Make grid of evaluation points to find 1's indifference point(s)
    n <- nrow(regr[[1]])
    x <- seq(-1, 1, length.out = gridsize)  # standardized grid
    x <- matrix(x, nrow = n, ncol = gridsize, byrow = TRUE)
    x <- x * with(utils, pmax(abs(u11 - u13), abs(u11 - u14)) + .1)

    ## Evaluate the function whose roots are 1's cutpoints over grid
    p3 <- with(utils, pnormCond(-u24, x, rho))
    y <- with(utils, p3 * u13 + (1 - p3) * u14 + x - u11)

    ## Compute number and grid location of cutpoints for each observation
    sydiff <- sign(y)[, 2:gridsize] - sign(y)[, 1:(gridsize-1)]
    cutgrid <- apply(sydiff, 1, function(m) which(m != 0))
    ncut <- sapply(cutgrid, length)
    only1 <- ncut == 1

    ## Compute approximate cutpoints via linear interpolation for observations
    ## with 1 cutpoint
    if (any(only1)) {
        ind1 <- unlist(cutgrid[only1])
        ind1a <- cbind(which(only1), ind1)
        ind1b <- cbind(which(only1), ind1 + 1)
        x1 <- cbind(x[ind1a], x[ind1b])
        y1 <- cbind(y[ind1a], y[ind1b])
        cut1 <- approxRoot(x1, y1)
    }

    ## Do the same for those with 3 cutpoints
    if (any(!only1)) {
        ind3 <- do.call(rbind, cutgrid[!only1])
        ind31a <- cbind(which(!only1), ind3[, 1])
        ind31b <- cbind(which(!only1), ind3[, 1] + 1)
        ind32a <- cbind(which(!only1), ind3[, 2])
        ind32b <- cbind(which(!only1), ind3[, 2] + 1)
        ind33a <- cbind(which(!only1), ind3[, 3])
        ind33b <- cbind(which(!only1), ind3[, 3] + 1)
        x31 <- cbind(x[ind31a], x[ind31b])
        y31 <- cbind(y[ind31a], y[ind31b])
        x32 <- cbind(x[ind32a], x[ind32b])
        y32 <- cbind(y[ind32a], y[ind32b])
        x33 <- cbind(x[ind33a], x[ind33b])
        y33 <- cbind(y[ind33a], y[ind33b])
        cut3 <- cbind(approxRoot(x31, y31), approxRoot(x32, y32),
                      approxRoot(x33, y33))
    }

    ## Compute outcome probabilities for each observation with 1 cutpoint
    ans <- matrix(NA, nrow = n, ncol = 3)
    if (any(only1)) {
        ans[only1, 1] <- pnorm(cut1)
        ans[only1, 2] <- pbivnorm(-cut1, -utils$u24[only1], -rho)
    }

    ## Compute outcome probabilities for each observation with 3 cutpoints
    if (any(!only1)) {
        ans[!only1, 1] <- pnorm(cut3[, 1]) + pnorm(cut3[, 3]) - pnorm(cut3[, 2])
        ans[!only1, 2] <- {
            pbivnorm(cut3[, 2], -utils$u24[!only1], rho)
            - pbivnorm(cut3[, 1], -utils$u24[!only1], rho)
            + pbivnorm(-cut3[, 3], -utils$u24[!only1], -rho)
        }
    }

    ans[, 1] <- finiteProbs(ans[, 1])
    ans[, 2] <- finiteProbs(ans[, 2])
    ans[, 3] <- finiteProbs(1 - ans[, 1] - ans[, 2])
    return(ans)
}

##
## Log-likelihood function
##
logLik12cor <- function(b, y, regr, gridsize, ...)
{
    logProbs <- log(makeProbs12cor(b, regr, gridsize))
    ans <- logProbs[cbind(seq_len(nrow(logProbs)), y)]
    return(ans)
}

egame12cor <- function(formulas, data, subset, na.action,
                       gridsize = 100,
                       startvals = c("sbi", "unif", "zero"),
                       boot = 0,
                       bootreport = TRUE,
                       profile,
                       method = "BFGS",
                       ...)
{
    cl <- match.call()
    startvals <- match.arg(startvals)

    ## Convert 'formulas' to appropriate class and check validity
    formulas <- checkFormulas(formulas)
    if (length(formulas)[2] != 4)
        stop("'formulas' should have four components on the right-hand side")

    ## Make the model frame (same as in 'egame12')
    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## Get the response and store it as factor 'yf' and numeric 'y' (same as in
    ## 'egame12')
    yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
    yf <- makeResponse12(yf)
    y <- as.numeric(yf)

    ## Create list of regressor matrices
    regr <- list()
    for (i in 1:4)
        regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
    rcols <- sapply(regr, ncol)

    ## Create starting values
    if (missing(profile) || is.null(profile)) {
        if (startvals == "zero") {
            sval <- rep(0, sum(rcols) + 1)
        } else if (startvals == "unif") {
            if (!hasArg(unif))
                unif <- c(-1, 1)
            sval <- runif(sum(rcols) + 1, unif[1], unif[2])
        } else {
            ## For SBI, must divide by sqrt(2) since we're now only including
            ## one error term per player
            sval <- sbi12(y, regr, link = "probit") / sqrt(2)
            sval <- c(sval, 0)
        }
    } else {
        sval <- svalsFromProfile(profile)
    }

    ## Identification check (same as in 'egame12')
    varNames <- lapply(regr, colnames)
    idCheck <- (varNames[[1]] %in% varNames[[2]])
    idCheck <- idCheck & (varNames[[1]] %in% varNames[[3]])
    if (any(idCheck)) {
        stop("Identification problem: the following variables appear in all three of player 1's utility equations: ",
             paste(varNames[[1]][idCheck], collapse = ", "))
    }

    ## Variable naming
    prefixes <- paste(c(rep("u1(", 3), "u2("), c(levels(yf), levels(yf)[3]),
                       ")", sep = "")
    varNames <- makeVarNames(varNames, prefixes, utils = 1:4, link = "probit",
                             sdterms = 0)
    hasColon <- c(varNames$hasColon, "rho (transformed)" = FALSE)
    names(sval) <- c(varNames$varNames, "rho (transformed)")
    fvec <- rep(FALSE, length(sval))  # No fixed terms

    ## Compute MLE
    results <- maxLik(logLik = logLik12cor, start = sval, fixed = fvec, method =
                      method, y = y, regr = regr, gridsize = gridsize,
                      tol.uniroot = tol.uniroot, mmm = mmm, ...)

    ## Check convergence and local identification
    cc <- convergenceCriterion(method)
    if (!(results$code %in% cc)) {
        warning("Model fitting did not converge\nCode:", results$code,
                "\nMessage: ", results$message)
    }
    lid <- checkLocalID(results$hessian, fvec)
    if (!lid)
        warning("Hessian is not negative definite; coefficients may not be locally identified")

    ## Bootstrap
    if (boot > 0) {
        bootMatrix <-
            gameBoot(boot, report = bootreport, estimate = results$estimate, y =
                     y, regr = regr, fn = logLik12cor, gr = NULL, fixed = fvec,
                     method = method, gridsize = gridsize, ...)
    }
    
    ## Store results
    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getGameVcov(results$hessian, fvec)
    ans$log.likelihood <- logLik12cor(results$estimate, y = y, regr = regr,
                                      gridsize = gridsize)
    ans$call <- cl
    ans$convergence <- list(method = method, iter = nIter(results), code =
                            results$code, message = results$message, gradient =
                            FALSE)
    ans$formulas <- formulas
    ans$link <- "probit"
    ans$type <- "agent"
    ans$model <- mf
    ans$xlevels <- .getXlevels(attr(mf, "terms"), mf)
    ans$y <- yf
    ans$equations <- names(hasColon)
    attr(ans$equations, "hasColon") <- hasColon
    ans$fixed <- fvec
    if (boot > 0)
        ans$boot.matrix <- bootMatrix
    ans$localID <- lid
    ans$gridsize <- gridsize
    class(ans) <- c("game", "egame12cor")

    return(ans)
}


## simple simulation
library(MASS)
library(pbivnorm)
library(maxLik)
library(Formula)

n <- 100
x24 <- runif(n)
u11 <- 0.5
u13 <- 0
u14 <- 1
u24 <- x24
e <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, .2, .2, 1), 2, 2))
p3 <- pnormCond(-u24, e[, 1], .2)
y1star <- p3 * u13 + (1 - p3) * u14 + e[, 1] - u11
y2star <- u24 + e[, 2]
y <- ifelse(y1star < 0, 1, ifelse(y2star < 0, 2, 3))
dat <- data.frame(y, x11, x14, x24)
m2 <- egame12cor(y ~ 1 | 0 | 1 | x24, data = dat, print.level = 3,
                 gridsize = 100)

## system.time({
## m2 <- egame12cor(y ~ x11 | 0 | x14 | x24, data = dat, print.level = 3,
##                  mmm = "a")
## })

## cf <- m1$estimate
## cfdiff <- seq(-1, 1, by = 0.1)
## ans <- numeric(length(cfdiff))
## for (i in seq_along(cfdiff)) {
##     newcf <- cf
##     newcf[length(newcf)] <- newcf[length(newcf)] + cfdiff[i]
##     ans[i] <- sum(logLik12cor(newcf, y, regr, 100, .Machine$double.eps^0.5))
## }

## k <- 0
## maxit <- 50
## while (!m5$localID && k < maxit)
## {
##     k <- k + 1
##     p5 <- profile(m5, dist = 0.1, steps = 25, use.se = FALSE)
##     p5x <- do.call(rbind, p5)
##     if (max(p5x[, 1], na.rm = TRUE) <= logLik(m5))
##         break
##     m5 <- update(m5, profile = p5)
## }
