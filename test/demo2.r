################################################################################
### File:   demo2.r
### Date:   2010-08-11
### Author: Brenton Kenkel
### 
### New demos: namely of strat4 function.
################################################################################

source("~/.strat/strat/R/strat.r")
source("~/.strat/strat/R/strat122.r")
library(Formula)
library(maxLik)
library(MASS)

n <- 1000
nmc <- 2
ans1 <- ans2 <- matrix(NA, nrow = nmc, ncol = 15)
pb <- txtProgressBar(min = 1, max = nmc)
for (i in 1:nmc) {
    for (j in 1:8) {
        assign(paste("x", j, sep = ""), rnorm(n))
        assign(paste("z", j, sep = ""), rnorm(n))
    }
    mcDataNames <- ls(pattern = "^[xz][0-9].*$")
    mcData <- lapply(mcDataNames, get)
    mcData <- do.call(data.frame, mcData)
    names(mcData) <- mcDataNames

    u11 <- 0.5 - 0.5 * x1 + 0.5 * x2
    u12 <- -0.5 + 0.5 * x3 - 0.5 * x4
    u13 <- 0
    u14 <- 2 + 2 * x5 - 3 * x6
    u22 <- 1 - 1 * z1 + 1 * z2
    u24 <- -1 + 1 * z4 - 1 * z5

    ## agent/logit
    p6 <- plogis(u24 / sqrt(2))
    p5 <- 1 - p6
    p4 <- plogis(u22 / sqrt(2))
    p3 <- 1 - p4
    y1 <- as.numeric(p3 * u11 + p4 * u12 + rlogis(n) < p5 * u13 + p6 * u14 +
                     rlogis(n))
    y2L <- as.numeric(u22 + rlogis(n) > rlogis(n))
    y2R <- as.numeric(u24 + rlogis(n) > rlogis(n))
    y <- rep(1L, n)
    y[y1 == 0 & y2L == 1] <- 2L
    y[y1 == 1 & y2R == 0] <- 3L
    y[y1 == 1 & y2R == 1] <- 4L

    mcDataAL <- mcData
    mcDataAL$y <- y

    ## probit/private
    p6 <- pnorm(u24)
    p5 <- 1 - p6
    p4 <- pnorm(u22)
    p3 <- 1 - p4
    y1 <- as.numeric(p3 * (u11 + rnorm(n)) + p4 * (u12 + rnorm(n)) <
                     p5 * (u13 + rnorm(n)) + p6 * (u14 + rnorm(n)))
    y2L <- as.numeric(u22 + rnorm(n) > 0)
    y2R <- as.numeric(u24 + rnorm(n) > 0)
    y <- rep(1L, n)
    y[y1 == 0 & y2L == 1] <- 2L
    y[y1 == 1 & y2R == 0] <- 3L
    y[y1 == 1 & y2R == 1] <- 4L
    mcDataPP <- mcData
    mcDataPP$y <- y

    m1 <- strat122(y ~ x1 + x2 | x3 + x4 | 0 | x5 + x6 | z1 + z2 | z4 + z5,
                 data = mcDataAL, type = "agent", link = "logit")

    m2 <- strat122(y ~ x1 + x2 | x3 + x4 | 0 | x5 + x6 | z1 + z2 | z4 + z5,
                 data = mcDataPP, type = "private", link = "probit")

    ans1[i, ] <- coef(m1)
    ans2[i, ] <- coef(m2)

    setTxtProgressBar(pb, i)
}
colnames(ans1) <- colnames(ans2) <- names(coef(m1))
