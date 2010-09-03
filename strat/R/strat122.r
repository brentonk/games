##' @include strat.r
NULL

##' @S3method predict strat122
predict.strat122 <- function(object, newdata, probs = c("outcome", "action"), ...)
{
    probs <- match.arg(probs)

    if (missing(newdata))
        newdata <- object$model

    mf <- match(c("subset", "na.action"), names(object$call), 0L)
    mf <- object$call[c(1L, mf)]
    mf$formula <- object$formulas
    mf$data <- newdata
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## This is to prevent goofy things from happening with factor variables --
    ## it ensures that the levels of factor variables in "newdata" are the same
    ## as those in the model frame used in the original fitting
    for (i in 1:ncol(mf)) {
        if (is.factor(mf[, i])) {
            iname <- names(mf)[i]
            if (iname %in% names(object$model))
                levels(mf[, i]) <- levels(object$model[, iname])
        }
    }

    regr <- list()
    for (i in seq_len(length(object$formulas)[2]))
        regr[[i]] <- model.matrix(object$formulas, data = mf, rhs = i)

    ans <- makeProbs122(object$coefficients, regr = regr, link = object$link, type
                        = object$type)

    if (probs == "outcome") {
        ans <- data.frame(actionsToOutcomes122(ans, log.p = FALSE))
        names(ans) <- paste("Pr(", levels(object$y), ")", sep = "")
    }

    ans <- do.call(cbind, ans)
    ans <- as.data.frame(ans)
    return(ans)
}

sbi122 <- function(y, regr, link)
{
    if (link == "probit") {
        fam <- binomial(link = "probit")
    } else {
        fam <- binomial(link = "logit")
    }

    ZL <- regr$Z2[y == 1 | y == 2, , drop = FALSE]
    yL <- as.numeric(y == 2)[y == 1 | y == 2]
    mL <- suppressWarnings(glm.fit(ZL, yL, family = fam))
    p2 <- as.numeric(regr$Z2 %*% coef(mL))
    p2 <- if (link == "probit") pnorm(p2) else plogis(p2)

    ZR <- regr$Z4[y == 3 | y == 4, , drop = FALSE]
    yR <- as.numeric(y == 4)[y == 3 | y == 4]
    mR <- suppressWarnings(glm.fit(ZR, yR, family = fam))
    p4 <- as.numeric(regr$Z4 %*% coef(mR))
    p4 <- if (link == "probit") pnorm(p4) else plogis(p4)

    X1 <- cbind(-(1-p2) * regr$X1, -p2 * regr$X2, (1 - p4) * regr$X3, p4 *
                regr$X4)
    y1 <- as.numeric(y == 3 | y == 4)
    m1 <- glm.fit(X1, y1, family = fam)

    ans <- sqrt(2) * c(coef(m1), coef(mL), coef(mR))
    return(ans)
}

makeUtils122 <- function(b, regr)
{
    utils <- vector("list", 6)
    names(utils) <- c("u11", "u12", "u13", "u14", "u22", "u24")

    rcols <- sapply(regr, ncol)
    for (i in 1:6) {
        if (rcols[i] > 0) {
            utils[[i]] <- as.numeric(regr[[i]] %*% b[1:rcols[i]])
            b <- b[-(1:rcols[i])]
        } else {
            utils[[i]] <- rep(0, nrow(regr[[i]]))
        }
    }

    utils$b <- b
    return(utils)
}

makeSDs122 <- function(b, regr)
{
    sds <- vector("list", 6)

    rcols <- sapply(regr, ncol)
    for (i in 1:6) {
        if (rcols[i+6] > 0) {
            sds[[i]] <- exp(as.numeric(regr[[i+6]] %*% b[1:rcols[i+6]]))
            b <- b[-(1:rcols[i+6])]
        } else {
            sds[[i]] <- rep(1, nrow(regr[[i+6]]))
        }
    }

    return(sds)
}

makeProbs122 <- function(b, regr, link, type)
{
    private <- type == "private"

    utils <- makeUtils122(b, regr)

    if (length(utils$b) == 0) {
        sds <- as.list(rep(1, 6))
    } else {
        sds <- makeSDs122(utils$b, regr)
    }

    linkfcn <- switch(link,
                      logit = function(x, sd = 1) plogis(x, scale = sd),
                      probit = pnorm)

    sd6 <- if (private) sds[[6]] else sqrt(sds[[5]]^2 + sds[[6]]^2)
    p6 <- finiteProbs(linkfcn(utils$u24, sd = sd6))
    p5 <- 1 - p6

    sd4 <- if (private) sds[[5]] else sqrt(sds[[3]]^2 + sds[[4]]^2)
    p4 <- finiteProbs(linkfcn(utils$u22, sd = sd4))
    p3 <- 1 - p4

    sd2 <- if (private) {
        sqrt(p3^2 * sds[[1]]^2 + p4^2 * sds[[2]]^2 + p5^2 * sds[[3]]^2 +
             p6^2 * sds[[4]]^2)
    } else sqrt(sds[[1]]^2 + sds[[2]]^2)
    p2 <- p5 * utils$u13 + p6 * utils$u14 - p3 * utils$u11 - p4 * utils$u12
    p2 <- finiteProbs(linkfcn(p2, sd = sd2))
    p1 <- 1 - p2

    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
}

actionsToOutcomes122 <- function(probs, log.p = TRUE)
{
    ans <- cbind(log(probs$p1) + log(probs$p3),
                 log(probs$p1) + log(probs$p4),
                 log(probs$p2) + log(probs$p5),
                 log(probs$p2) + log(probs$p6))
    if (!log.p) ans <- exp(ans)
    return(ans)
}

logLik122 <- function(b, y, regr, link, type, ...)
{
    probs <- makeProbs122(b, regr, link, type)
    logProbs <- actionsToOutcomes122(probs, log.p = TRUE)
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

logLikGrad122 <- function(b, y, regr, link, type, ...)
{
    u <- makeUtils122(b, regr)
    p <- makeProbs122(b, regr, link, type)
    rcols <- sapply(regr, ncol)
    n <- nrow(regr$X1)

    if (link == "probit" && type == "private") {
        dp4db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp4dg2 <- dnorm(u$u22) * regr$Z2
        dp4dg4 <- matrix(0L, nrow = n, ncol = rcols[6])
        dp4 <- cbind(dp4db, dp4dg2, dp4dg4)
        dp3 <- -dp4

        dp6db <- dp4db
        dp6dg2 <- matrix(0L, nrow = n, ncol = rcols[5])
        dp6dg4 <- dnorm(u$u24) * regr$Z4
        dp6 <- cbind(dp6db, dp6dg2, dp6dg4)
        dp5 <- -dp6

        num2 <- p$p5 * u$u13 + p$p6 * u$u14 - p$p3 * u$u11 - p$p4 * u$u12
        denom2 <- sqrt(p$p3^2 + p$p4^2 + p$p5^2 + p$p6^2)
        dn2 <- dnorm(num2 / denom2)
        dp2db1 <- dn2 * (-p$p3) * regr$X1 / denom2
        dp2db2 <- dn2 * (-p$p4) * regr$X2 / denom2
        dp2db3 <- dn2 * p$p5 * regr$X3 / denom2
        dp2db4 <- dn2 * p$p6 * regr$X4 / denom2
        dp2dg2 <- dn2 * ((u$u11 - u$u12) * denom2 - num2 * (p$p4 - p$p3) / denom2)
        dp2dg2 <- (dp2dg2 * dp4dg2) / denom2^2
        dp2dg4 <- dn2 * ((u$u14 - u$u13) * denom2 - num2 * (p$p6 - p$p5) / denom2)
        dp2dg4 <- (dp2dg4 * dp6dg4) / denom2^2
        dp2 <- cbind(dp2db1, dp2db2, dp2db3, dp2db4, dp2dg2, dp2dg4)
        dp1 <- -dp2
    } else if (type == "agent") {
        dlink <- switch(link,
                        logit = dlogis,
                        probit = dnorm)

        dp4db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp4dg2 <- dlink(u$u22 / sqrt(2)) * regr$Z2 / sqrt(2)
        dp4dg4 <- matrix(0L, nrow = n, ncol = rcols[6])
        dp4 <- cbind(dp4db, dp4dg2, dp4dg4)
        dp3 <- -dp4

        dp6db <- dp4db
        dp6dg2 <- matrix(0L, nrow = n, ncol = rcols[5])
        dp6dg4 <- dlink(u$u24 / sqrt(2)) * regr$Z4 / sqrt(2)
        dp6 <- cbind(dp6db, dp6dg2, dp6dg4)
        dp5 <- -dp6

        dn2 <- dlink((p$p5 * u$u13 + p$p6 * u$u14 - p$p3 * u$u11 - p$p4 * u$u12)
                     / sqrt(2))
        dp2db1 <- dn2 * (-p$p3 / sqrt(2)) * regr$X1
        dp2db2 <- dn2 * (-p$p4 / sqrt(2)) * regr$X2
        dp2db3 <- dn2 * (p$p5 / sqrt(2)) * regr$X3
        dp2db4 <- dn2 * (p$p6 / sqrt(2)) * regr$X4
        dp2dg2 <- dn2 * (u$u11 - u$u12) * dlink(u$u22 / sqrt(2)) * regr$Z2 / 2
        dp2dg4 <- dn2 * (u$u14 - u$u13) * dlink(u$u24 / sqrt(2)) * regr$Z4 / 2
        dp2 <- cbind(dp2db1, dp2db2, dp2db3, dp2db4, dp2dg2, dp2dg4)
        dp1 <- -dp2
    }

    dL1 <- (1 / p$p1) * dp1 + (1 / p$p3) * dp3
    dL2 <- (1 / p$p1) * dp1 + (1 / p$p4) * dp4
    dL3 <- (1 / p$p2) * dp2 + (1 / p$p5) * dp5
    dL4 <- (1 / p$p2) * dp2 + (1 / p$p6) * dp6

    ans <- matrix(NA, nrow = n, ncol = sum(rcols[1:6]))
    ans[y == 1, ] <- dL1[y == 1, ]
    ans[y == 2, ] <- dL2[y == 2, ]
    ans[y == 3, ] <- dL3[y == 3, ]
    ans[y == 4, ] <- dL4[y == 4, ]

    return(ans)
}

makeResponse122 <- function(yf)
{
    if (length(dim(yf))) {
        Y <- yf

        if (ncol(Y) > 3) {
            warning("only first three columns of response will be used")
            Y <- Y[, 1:3]
        }

        if (!identical(sort(unique(unlist(Y))), c(0, 1)))
            stop("dummy responses must be dummy variables")

        if (ncol(Y) == 3) {
            Y[, 2] <- ifelse(Y[, 1] == 1, Y[, 3], Y[, 2])
            ylevs <- c(paste("~", names(Y)[2], sep = ""),
                       names(Y)[2],
                       paste("~", names(Y)[3], sep = ""),
                       names(Y)[3])
        } else {
            ylevs <- c(paste("~", names(Y)[1], ",~", names(Y)[2], sep = ""),
                       paste("~", names(Y)[1], ",", names(Y)[2], sep = ""),
                       paste(names(Y)[1], ",~", names(Y)[2], sep = ""),
                       paste(names(Y)[1], ",", names(Y)[2], sep = ""))
        }

        y <- numeric(nrow(Y))
        y[Y[, 1] == 0 & Y[, 2] == 0] <- 1
        y[Y[, 1] == 0 & Y[, 2] == 1] <- 2
        y[Y[, 1] == 1 & Y[, 2] == 0] <- 3
        y[Y[, 1] == 1 & Y[, 2] == 1] <- 4

        yf <- as.factor(y)
        levels(yf) <- ylevs
    } else {
        yf <- as.factor(yf)
        if (nlevels(yf) != 4) stop("dependent variable must have four values")
    }

    return(yf)
}

##' <description>
##'
##' \preformatted{
##' .        ___ 1 ___
##' .       /         \
##' .      /           \
##' .   2 /             \ 2
##' .    / \           / \
##' .   /   \         /   \
##' .  /     \       /     \
##' . u11    u12    u13    u14
##' . 0      u22    0      u24}
##' @title
##' @param formulas
##' @param data
##' @param subset
##' @param na.action
##' @param varformulas
##' @param fixedUtils
##' @param link
##' @param type
##' @param startvals
##' @param ...
##' @return
##' @author Brenton Kenkel
strat122 <- function(formulas, data, subset, na.action,
                     link = c("probit", "logit"),
                     type = c("agent", "private"),
                     startvals = c("sbi", "unif", "zero"),
                     fixedUtils = NULL,
                     sdformula = NULL,
                     sdByPlayer = FALSE,
                     boot = 0,
                     bootreport = TRUE,
                     ...)
{
    cl <- match.call()

    link <- match.arg(link)
    type <- match.arg(type)
    startvals <- match.arg(startvals)

    formulas <- checkFormulas(formulas)

    if (!is.null(fixedUtils)) {
        if (length(fixedUtils) < 6)
            stop("fixedUtils must have 6 elements (u11, u12, u13, u14, u22, u24)")
        if (length(fixedUtils) > 6) {
            warning("only the first 6 elements of fixedUtils will be used")
            fixedUtils <- fixedUtils[1:6]
        }

        if (startvals == "sbi")
            startvals <- "zero"

        if (is.null(sdformula))
            sdformula <- if (sdByPlayer) Formula(~ 1 | 1) else Formula(~ 1)
    }

    if (!is.null(sdformula)) {
        sdformula <- checkFormulas(sdformula, argname = "sdformula")
        if (sdByPlayer && length(sdformula)[2] != 2)
            stop("`sdformula` should have two components (one for each player) on the right-hand side when sdByPlayer == TRUE")
        if (!sdByPlayer && length(sdformula)[2] != 1)
            stop("`sdformula` should have exactly one component on the right-hand side")
        formulas <- as.Formula(formula(formulas), formula(sdformula))
    }

    if (sdByPlayer && is.null(sdformula)) {
        warning("to estimate SDs by player, you must specify `sdformula` or `fixedUtils`")
        sdByPlayer <- FALSE
    }

    if (link == "logit" && type == "private") {
        warning("logit link cannot be used with private information model; changing to probit link")
        link <- "probit"
    }

    ## make the model frame
    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
    yf <- makeResponse122(yf)
    y <- as.numeric(yf)

    regr <- list()
    for (i in seq_len(length(formulas)[2]))
        regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
    names(regr) <- character(length(regr))
    names(regr)[1:6] <- c("X1", "X2", "X3", "X4", "Z2", "Z4")
    rcols <- sapply(regr, ncol)

    ## starting values
    if (startvals == "zero") {
        sval <- rep(0, sum(rcols))
    } else if (startvals == "unif") {
        if (!hasArg(unif))
            unif <- c(-1, 1)
        sval <- runif(sum(rcols), unif[1], unif[2])
    } else {
        sval <- sbi122(y, regr, link)
        sval <- c(sval, rep(0, sum(rcols) - length(sval)))
    }

    ## identification check
    varNames <- sapply(regr, colnames)
    idCheck <- do.call(intersectAll, varNames[1:4])
    if (is.null(fixedUtils) && length(idCheck) > 0) {
        stop("Identification problem: the following variables appear in all four of player 1's utility equations: ",
             paste(idCheck, collapse = ", "))
    }

    prefixes <- paste(c(rep("u1(", 4), rep("u2(", 2)),
                      c(levels(yf), levels(yf)[2], levels(yf)[4]), ")",
                      sep = "")
    sdterms <- if (!is.null(sdformula)) { if (sdByPlayer) 2L else 1L } else 0L
    varNames <- makeVarNames(varNames, prefixes, link, sdterms)
    hasColon <- varNames$hasColon
    names(sval) <- varNames$varNames

    gr <- if (is.null(sdformula)) logLikGrad122 else NULL

    fvec <- rep(FALSE, length(sval))
    names(fvec) <- names(sval)
    if (!is.null(fixedUtils)) {
        sval[1:6] <- fixedUtils
        fvec[1:6] <- TRUE
    }

    results <- maxBFGS(fn = logLik122, grad = gr, start = sval, fixed = fvec,
                       y = y, regr = regr, link = link, type = type, ...)
    if (results$code) {
        warning("Model fitting did not converge\nMessage: ",
                results$message)
    }

    if (boot > 0) {
        bootMatrix <- stratBoot(boot, report = bootreport, estimate =
                                results$estimate, y = y, regr = regr, fn =
                                logLik122, gr = gr, fixed = fvec, link = link,
                                type = type, ...)
    }

    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getStratVcov(results$hessian, fvec)
    ans$log.likelihood <-
        logLik122(results$estimate, y = y, regr = regr, link = link, type = type)
    ans$call <- cl
    ans$convergence <- list(code = results$code, message = results$message)
    ans$formulas <- formulas
    ans$link <- link
    ans$type <- type
    ans$model <- mf
    ans$y <- yf
    ans$equations <- prefixes
    attr(ans$equations, "hasColon") <- hasColon
    names(attr(ans$equations, "hasColon")) <- ans$equations
    ans$fixed <- fvec
    if (boot > 0)
        ans$boot.matrix <- bootMatrix

    class(ans) <- c("strat", "strat122")

    return(ans)
}
