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
rm -f games.roxygen/R/.Rhistory

# for including S3 method directives for undocumented methods
cat ./games/NAMESPACE >> ./games.roxygen/NAMESPACE

# add acknowledgments
echo "\section{Acknowledgments}{We thank the Wallis Institute of Political Economy for financial support.}" >> ./games.roxygen/man/games-package.Rd

R CMD check games.roxygen
R CMD build games.roxygen
