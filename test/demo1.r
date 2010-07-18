################################################################################
### File:   demo1.r
### Date:   2010-07-11
### Author: Brenton Kenkel
###
### First demo of "strat" R package.
################################################################################

source("~/.strat/strat/R/strat.r")
library(Formula)
library(maxLik)

## make dataset with logit link and agent error
n <- 1000
x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
z1 <- rnorm(n); z2 <- rnorm(n)
w1 <- rnorm(n); w2 <- rnorm(n)

u24 <- 3 - 2 * z1 + z2
u14 <- 2 * x4
u13 <- -2 + 3 * x3
u11 <- 1 - x1 + 2 * x2

p4 <- plogis(u24 / sqrt(2))
y1 <- as.numeric(u11 + rlogis(n) < (1 - p4) * u13 + p4 * u14 + rlogis(n))
y2 <- as.numeric(u24 + rlogis(n) > rlogis(n))
y <- rep(1L, n)
y[y1 == 1] <- 3
y[y1 == 1 & y2 == 1] <- 4

mcData <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, z1 = z1, z2 = z2, w1 =
                     w1, w2 = w2)
mcData$y1 <- y1
mcData$y2 <- y2
mcData$y <- y
mcData$yf <- factor(y, levels = c(1, 3, 4), labels = c("sq", "cap", "war"))


## the main estimation function is called "strat3" (since there are 3 terminal
## nodes).  I coded up the gradient, so the estimation is relatively fast and
## stable.

## the regressors are specified via four formulas in R:
##   outcome ~ u1(1) vars | u1(3) vars | u1(4) vars | u2(4) vars

## it's also possible to use a list of formulas if you don't like the vertical
## bar notation:
##   list(outcome ~ u1(1) vars,
##        ~ u1(3) vars,
##        ~ u1(4) vars,
##        ~ u2(4) vars)

## the outcome can be a numeric vector with three values, a factor with three
## levels, or a set of two dummy variables (first is L or R for player 1, second
## is L or R for player 2)

## I have the machinery in place to parameterize the variance terms with
## regressors, but I haven't yet tested it and highly doubt that it works
## properly.  so ignore the "varformulas" option for now.


## formulas all together, outcome as numeric
m1 <- strat3(y ~ x1 + x2 | x3 | x4 - 1 | z1 + z2, data = mcData, link =
             "logit", type = "agent", varformulas = ~ w1 | w2 | 0 | 0)
print(m1)
summary(m1)


## formulas all together, outcome as factor (notice that the levels of the
## factor are used as names in the output)
m2 <- strat3(yf ~ x1 + x2 | x3 | x4 - 1 | z1 + z2, data = mcData, link =
             "logit", type = "agent")
print(m2)
summary(m2)


## formulas all together, outcome as dummy variables
m3 <- strat3(y1 + y2 ~ x1 + x2 | x3 | x4 - 1 | z1 + z2, data = mcData, link =
             "logit", type = "agent")
print(m3)
summary(m3)


## formulas as list, outcome as numeric
m4 <- strat3(list(y ~ x1 + x2, ~ x3, ~ x4 - 1, ~ z1 + z2), data = mcData, link
             = "logit", type = "agent")
print(m4)
summary(m4)


## starting values are taken from SBI by default, but you can also have them
## chosen from U(-1, 1) or to be all zero
m5 <- strat3(yf ~ x1 + x2 | x3 | x4 - 1 | z1 + z2, data = mcData, link =
             "logit", type = "agent", startvals = "zero")
m6 <- strat3(yf ~ x1 + x2 | x3 | x4 - 1 | z1 + z2, data = mcData, link =
             "logit", type = "agent", startvals = "unif")
print(m5)
print(m6)


## for simplicity, I'm using an S3 class to store the output from the estimation
## function.  the methods written so far are just the basics (notably, I have
## yet to finish up "predict"):
print(m1)
summary(m1)
coef(m1)
vcov(m1)
logLik(m1)
AIC(m1)
library(stats4); BIC(m1)
