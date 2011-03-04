## simdatault.r -- for simulating the dataset that goes with the ultimatum model

n <- 1500
for (i in 1:4) {
    assign(paste("x", i, sep = ""), rnorm(n))
    assign(paste("z", i, sep = ""), rnorm(n))
}
for (i in 1:2)
    assign(paste("w", i, sep = ""), rnorm(n))

s1 <- 5
s2 <- 1
R1s <- 5 - 2 * x1 + 2 * x2 - 2 * x3 + 2 * x4 + w1 - w2
R2s <- 10 + 2 * z1 - 2 * z2 + 2 * z3 - 2 * z4 - w1 + w2
R1 <- R1s + rlogis(n, scale = s1)
R2 <- R2s + rlogis(n, scale = s2)

Q <- 15
y.star <- Q - R1 - s2 * (1 + LW(exp((Q - R1 - s2 - R2s) / s2)))
y <- ifelse(y.star < 0, 0L, ifelse(y.star > Q, Q, y.star))
accept <- as.numeric(y >= R2)

simultnames <- c("y", "accept", ls(pattern = "^[xwz][0-9]"))
simult <- lapply(simultnames, get)
simult <- do.call(data.frame, simult)
names(simult) <- simultnames
names(simult)[1] <- "offer"
data_ult <- simult

save(data_ult, file = "~/.games/games/data/data_ult.rda")
