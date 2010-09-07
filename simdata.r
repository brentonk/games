## simdata.r -- for simulating datasets included in the strat package

n <- 1000
for (i in 1:5)
    assign(paste("x", i, sep = ""), rnorm(n))
for (i in 1:3)
    assign(paste("z", i, sep = ""), rnorm(n))
for (i in 1:2) {
    assign(paste("f", i, sep = ""),
           as.factor(sample(letters[1:3], n, replace = TRUE)))
    f <- get(paste("f", i, sep = ""))
    for (j in 1:3)
        assign(paste("d", i, letters[j], sep = ""), as.numeric(f == letters[j]))
}


u11 <- 0.5 - 0.5 * x1 + 0.5 * x2
u12 <- -0.5 + 0.5 * x3 - 0.5 * d1b + 0.5 * d1c
u13 <- 0
u14 <- 0.5 - 0.5 * x4 + 0.5 * x5
u22 <- 1 - 1 * z1 + 1 * z2
u24 <- -1 + 1 * z3 - 1 * d2b + 1 * d2c

p6 <- plogis(u24 / sqrt(2))
p5 <- 1 - p6
p4 <- plogis(u22 / sqrt(2))
p3 <- 1 - p4
a1 <- as.numeric(p3 * u11 + p4 * u12 + rlogis(n) < p5 * u13 + p6 * u14 +
                 rlogis(n))
a2L <- as.numeric(u22 + rlogis(n) > rlogis(n))
a2R <- as.numeric(u24 + rlogis(n) > rlogis(n))
a2 <- ifelse(a1 == 1, a2R, a2L)
y <- ifelse(a1 == 1, ifelse(a2 == 1, "RR", "RL"), ifelse(a2 == 1, "LR", "LL"))
y <- as.factor(y)

sim122names <- ls(pattern = "^[fxz][0-9]")
sim122 <- lapply(sim122names, get)
sim122 <- do.call(data.frame, sim122)
names(sim122) <- sim122names
sim122$a1 <- a1
sim122$a2 <- a2
sim122$y <- y

save(sim122, file = "~/.strat/strat/data/sim122.rda")
