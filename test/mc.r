source("../strat/R/strat.r")

n <- 1000
nmc <- 100

ans1 <- ans2 <- ans3 <- matrix(nrow = nmc, ncol = 9)

for (i in seq_len(nmc)) {
    mcData <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 =
                         rnorm(n), z1 = rnorm(n), z2 = rnorm(n))

    attach(mcData)
    u24 <- 3 - 2 * z1 + z2
    u14 <- 3 * x4
    u13 <- 2 - 3 * x3
    u11 <- 1 - x1 + 2 * x2
    detach(mcData)

    ## probit, agent
    p4 <- pnorm(u24 / sqrt(2))
    y1 <- as.numeric(u11 + rnorm(n) < (1 - p4) * u13 + p4 * u14 + rnorm(n))
    y2 <- as.numeric(u24 + rnorm(n) > rnorm(n))
    y <- rep(1L, n)
    y[y1 == 1] <- 3
    y[y1 == 1 & y2 == 1] <- 4
    mcData$y <- y
    m1 <- strat22(y ~ x1 + x2 | x3 | x4 - 1 | z1 + z2,
                  data = mcData, link = "probit", type = "agent")
    ans1[i, ] <- coef(m1)

    ## logit, agent
    p4 <- plogis(u24 / sqrt(2))
    y1 <- as.numeric(u11 + rlogis(n) < (1 - p4) * u13 + p4 * u14 + rlogis(n))
    y2 <- as.numeric(u24 + rlogis(n) > rlogis(n))
    y <- rep(1L, n)
    y[y1 == 1] <- 3
    y[y1 == 1 & y2 == 1] <- 4
    mcData$y <- y
    m2 <- strat22(y ~ x1 + x2 | x3 | x4 - 1 | z1 + z2,
                  data = mcData, link = "logit", type = "agent")
    ans2[i, ] <- coef(m2)

    ## probit, private info
    p4 <- pnorm(u24)
    y1 <- as.numeric(u11 + rnorm(n) <
                     (1 - p4) * (u13 + rnorm(n)) + p4 * (u14 + rnorm(n)))
    y2 <- as.numeric(u24 + rnorm(n) > 0)
    y <- rep(1L, n)
    y[y1 == 1] <- 3
    y[y1 == 1 & y2 == 1] <- 4
    mcData$y <- y
    m3 <- strat22(y ~ x1 + x2 | x3 | x4 - 1 | z1 + z2,
                  data = mcData, link = "probit", type = "private")
    ans3[i, ] <- coef(m3)

    cat(i, "\n")
}

colMeans(ans1)
colMeans(ans2)
colMeans(ans3)
