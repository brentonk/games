##' @include games.r
##' @include helpers.r
NULL

sbi123 <- function(y, regr, link)
{
    names(regr) <- character(length(regr))
    names(regr)[1:8] <- c("X1", "X3", "X5", "X6", "Z3", "Z5", "Z6", "W6")

    if (link == "probit") {
        fam <- binomial(link = "probit")
        linkfcn <- pnorm
    } else {
        fam <- binomial(link = "logit")
        linkfcn <- plogis
    }

    reg3 <- regr$W6[y == 3 | y == 4, , drop = FALSE]
    y3 <- as.numeric(y == 4)[y == 3 | y == 4]
    m3 <- suppressWarnings(glm.fit(reg3, y3, family = fam))
    p6 <- as.numeric(regr$W6 %*% coef(m3))
    p6 <- linkfcn(p6)

    reg2 <- cbind(-regr$Z3, (1-p6) * regr$Z5, p6 * regr$Z6)
    reg22 <- reg2[y != 1, , drop = FALSE]
    y22 <- as.numeric(y != 2)[y != 1]
    m22 <- suppressWarnings(glm.fit(reg22, y22, family = fam))
    p4 <- as.numeric(reg2 %*% coef(m22))
    p4 <- linkfcn(p4)

    reg1 <- cbind(-regr$X1, (1-p4) * regr$X3, p4 * (1-p6) * regr$X5,
                  p4 * p6 * regr$X6)
    y1 <- as.numeric(y != 1)
    m1 <- glm.fit(reg1, y1, family = fam)

    ans <- sqrt(2) * c(coef(m1), coef(m22), coef(m3))
    return(ans)
}

makeSDs123 <- function(b, regr, type)
{
    sds <- vector("list", if (type == "private") 8L else 6L)
    regr <- regr[-(1:8)]
    rcols <- sapply(regr, ncol)

    if (length(rcols) == 1L) {  ## sdByPlayer == FALSE
        v <- exp(as.numeric(regr[[1]] %*% b))
        for (i in 1:length(sds)) sds[[i]] <- v
    } else {
        b1 <- b[1:rcols[1]]
        b2 <- b[(rcols[1]+1):(rcols[1]+rcols[2])]
        b3 <- b[(rcols[1]+rcols[2]+1):length(b)]
        v1 <- exp(as.numeric(regr[[1]] %*% b1))
        v2 <- exp(as.numeric(regr[[2]] %*% b2))
        v3 <- exp(as.numeric(regr[[3]] %*% b3))

        if (type == "private") {
            sds[[1]] <- sds[[2]] <- sds[[3]] <- sds[[4]] <- v1
            sds[[5]] <- sds[[6]] <- sds[[7]] <- v2
            sds[[8]] <- v3
        } else {
            sds[[1]] <- sds[[2]] <- v1
            sds[[3]] <- sds[[4]] <- v2
            sds[[5]] <- sds[[6]] <- v3
        }
    }

    return(sds)
}

makeProbs123 <- function(b, regr, link, type)
{
    utils <- makeUtils(b, regr, nutils = 8,
                       unames = c("u11", "u13", "u15", "u16", "u23", "u25",
                       "u26", "u36"))

    if (length(utils$b) == 0) {  ## variance unparameterized
        sds <- as.list(rep(1, 8))
    } else {
        sds <- makeSDs123(utils$b, regr, type)
    }

    linkfcn <- switch(link,
                      logit = function(x, sd = 1) plogis(x, scale = sd),
                      probit = pnorm)

    sd6 <- if (type == "private") sds[[8]] else sqrt(sds[[5]]^2 + sds[[6]]^2)
    p6 <- finiteProbs(linkfcn(utils$u36, sd = sd6))
    p5 <- 1 - p6

    if (type == "private") {
        sd4 <- sqrt(p5^2 * sds[[6]]^2 + p6^2 * sds[[7]]^2 + sds[[5]]^2)
    } else {
        sd4 <- sqrt(sds[[3]]^2 + sds[[4]]^2)
    }
    p4 <- p5 * utils$u25 + p6 * utils$u26 - utils$u23
    p4 <- finiteProbs(linkfcn(p4, sd = sd4))
    p3 <- 1 - p4

    if (type == "private") {
        sd2 <- sqrt(p3^2 * sds[[2]]^2 + p4^2 * p5^2 * sds[[3]]^2 +
                    p4^2 * p6^2 * sds[[4]]^2 + sds[[1]]^2)
    } else {
        sd2 <- sqrt(sds[[1]]^2 + sds[[2]]^2)
    }
    p2 <- p3 * utils$u13 + p4 * (p5 * utils$u15 + p6 * utils$u16) - utils$u11
    p2 <- finiteProbs(linkfcn(p2, sd = sd2))
    p1 <- 1 - p2

    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
}

actionsToOutcomes123 <- function(probs, log.p = TRUE)
{
    probs <- log(do.call(cbind, probs))
    ans <- cbind(probs[, 1],
                 probs[, 2] + probs[, 3],
                 probs[, 2] + probs[, 4] + probs[, 5],
                 probs[, 2] + probs[, 4] + probs[, 6])

    if (!log.p) ans <- exp(ans)
    return(ans)
}

logLik123 <- function(b, y, regr, link, type, ...)
{
    probs <- makeProbs123(b, regr, link, type)
    logProbs <- actionsToOutcomes123(probs, log.p = TRUE)
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

makeResponse123 <- function(yf)
{
    if (length(dim(yf))) {  ## yf is a matrix of dummies
        if (ncol(yf) == 2) {
            stop("response must be specified as a single vector or three dummy variables")
        } else if (ncol(yf) > 3) {
            warning("only first three columns of response will be used")
            yf <- yf[, 1:3]
        }

        if (!(all(unlist(yf) %in% c(0L, 1L))))
            stop("dummy responses must be dummy variables")

        ylevs <- c(paste("~", names(yf)[1], sep = ""),
                   paste(names(yf)[1], ",~", names(yf)[2], sep = ""),
                   paste(names(yf)[1], ",", names(yf)[2], ",~", names(yf)[3],
                         sep = ""),
                   paste(names(yf)[1], names(yf)[2], names(yf)[3], sep = ","))

        y <- integer(nrow(yf))
        y[yf[, 1] == 0] <- 1L
        y[yf[, 1] == 1 & yf[, 2] == 0] <- 2L
        y[yf[, 1] == 1 & yf[, 2] == 1 & yf[, 3] == 0] <- 3L
        y[yf[, 1] == 1 & yf[, 2] == 1 & yf[, 3] == 1] <- 4L
        yf <- as.factor(y)
        levels(yf) <- ylevs
    } else {                ## yf is a vector
        yf <- as.factor(yf)
        if (nlevels(yf) != 4) stop("dependent variable must have four values")
    }

    return(yf)
}

##' \preformatted{
##' .     1
##' .     /\
##' .    /  \
##' .   /    \ 2
##' .  u11   /\
##' .       /  \
##' .      /    \
##' .    u13     \ 3
##' .    u23     /\
##' .           /  \
##' .          /    \
##' .         u15   u16
##' .         u25   u26
##' .         0     u36}
