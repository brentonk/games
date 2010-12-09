#!/bin/bash

rm -rf ./games.roxygen

/usr/bin/R --no-save <<EOF
library(roxygen)
roxygenize("games", use.Rd2 = TRUE)
q()
EOF

rmdir games.roxygen/inst/doc
rmdir games.roxygen/inst
rm -f games.roxygen/.Rhistory
rm -f games.roxygen/.Rhistorynew

# for including S3 method directives for undocumented methods
cat ./games/NAMESPACE >> ./games.roxygen/NAMESPACE

R CMD check games.roxygen
R CMD build games.roxygen
