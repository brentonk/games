## Replicate analysis in Leblang (2003)
data("leblang2003")
m1 <- egame12(outcome ~
              capcont + lreserves + overval + creditgrow + USinterest + service
              + contagion + prioratt - 1 |
              1 |
              1 |
              unifgov + lexports + preelec + postelec + rightgov + realinterest
              + capcont + lreserves,
              data = leblang2003,
              link = "probit",
              type = "private")

summary(m1)
