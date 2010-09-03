### demo3.r -- testing the predicted probability functions

setwd("~/.strat/test")
library(Formula)
library(maxLik)
library(MASS)
source("~/.strat/strat/R/strat.r")
source("~/.strat/strat/R/strat12.r")

library(foreign)
war1800 <- read.dta("~/.505/hw/set9/war1800.dta")
war1800$esc <- 1 - war1800$sq

regime1 <- ifelse(war1800$dem1 <= -7, "autoc",
                  ifelse(war1800$dem1 >= 7, "dem", "mixed"))
regime2 <- ifelse(war1800$dem2 <= -7, "autoc",
                  ifelse(war1800$dem2 >= 7, "dem", "mixed"))
war1800$regime1 <- as.factor(regime1)
war1800$regime2 <- as.factor(regime2)

m1 <- strat12(esc + war ~ s_wt_re1 + revis1 | 0 | balanc + regime1 | balanc +
              regime2, sdformula = ~ s_wt_re1 | s_wt_re2, data = war1800, link
              = "logit", type = "agent", sdByPlayer = TRUE)

m2 <- strat12(esc + war ~ 0, fixedUtils = rnorm(4), data = war1800, link
              = "logit", type = "agent", sdByPlayer = TRUE)

m1pp1 <- predProbs(m1, "balanc", makePlot = FALSE)
plot(m1pp1, ask = TRUE)

m1pp2 <- predProbs(m1, "balanc", s_wt_re1 = quantile(s_wt_re1, 0.1), makePlot = FALSE)
plot(m1pp2, ask = TRUE)

m1pp3 <- predProbs(m1, "regime1", makePlot = FALSE)
par(mfrow = c(2, 2))
plot(m1pp3)
