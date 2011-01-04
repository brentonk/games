##' @include games.r
##' @include helpers.r
NULL

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
