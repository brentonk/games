data("data_122")

## Model formula:
fr1 <- y ~ x1 + x2 | x3 + f1 | 0 | x4 + x5 | z1 + z2 | z3 + f2
##     ^   ^^^^^^^   ^^^^^^^   ^   ^^^^^^^   ^^^^^^^   ^^^^^^^
##     y     u11       u12    u13    u14       u22       u24

m1 <- egame122(fr1, data = data_122)
summary(m1)

## Dummy specification of the dependent variable
fr2 <- update(Formula(fr1), a1 + a2 ~ .)
m2 <- egame122(fr2, data = data_122)
summary(m2)

## Estimation of scale parameters
fr3 <- y ~ x1 | x2 | 0 | x3 | z1 | z2
m3 <- egame122(fr3, data = data_122, sdformula = ~ x4 + z3 - 1)
summary(m3)

## Fixed utilities
utils <- c(0.25, -0.25, 0, 0.25, 0.5, -0.5)
m4 <- egame122(y ~ 1, data = data_122, fixedUtils = utils)
summary(m4)
