#!/bin/bash

rm -rf ./strat.roxygen

/usr/bin/R --no-save <<EOF
library(roxygen)
roxygenize("strat", use.Rd2 = TRUE)
q()
EOF

rmdir strat.roxygen/inst/doc
rmdir strat.roxygen/inst

R CMD check strat.roxygen
R CMD build strat.roxygen
