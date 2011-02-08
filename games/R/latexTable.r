##' @include games.r
##' @include helpers.r
NULL

##
## Converts a character vector for use in LaTeX by inserting escape sequences
## where appropriate.  Not comprehensive, but should catch most common
## problems.
## 
latexEsc <- function(x)
{
    x <- gsub("{", "\\{", x, fixed = TRUE)
    x <- gsub("}", "\\}", x, fixed = TRUE)
    x <- gsub("_", "\\_", x, fixed = TRUE)
    x <- gsub("#", "\\#", x, fixed = TRUE)
    x <- gsub("$", "\\$", x, fixed = TRUE)
    x <- gsub("%", "\\%", x, fixed = TRUE)
    x <- gsub("^", "\\^", x, fixed = TRUE)
    x <- gsub("&", "\\&", x, fixed = TRUE)
    x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
    return(x)
}

##' Makes a LaTeX table of strategic model results.
##'
##' \code{latexTable} prints LaTeX code for the presentation of results from a
##' strategic model.  Each row contains one regressor, and each column contains
##' one of the utility (or variance term) equations in the model.  For example,
##' in a model fit with \code{\link{egame12}}, the four columns will be u11,
##' u13, u14, and u24 respectively.  Each cell contains the estimated parameter,
##' atop its standard error in parentheses.  Cells that are empty because a
##' regressor is not included in a particular equation are filled with the
##' string specified in the option \code{blankfill}.  Signorino and Tarar (2006,
##' p. 593) contains a table of this variety.
##'
##' The table generated depends on the \pkg{multirow} package for LaTeX, so
##' make sure to include \code{\\usepackage{multirow}} in the preamble of your
##' document.
##' 
##' The \code{digits} option does not yet work seamlessly; you may have to
##' resort to trial and error.
##' @title LaTeX table for strategic models
##' @param x a fitted model of class \code{game}.
##' @param digits number of digits to print.
##' @param scientific logical or integer value to control use of scientific
##' notation.  See \code{\link{format}}.
##' @param blankfill text to fill empty cells (those where a certain variable
##' did not enter the given equation).
##' @param math.style.negative whether negative signs should be "math style" or
##' plain hyphens.  Defaults to \code{TRUE}.
##' @param file file to save the output in.  Defaults to \code{""}, which prints
##' the table to the R console.
##' @param floatplace where to place the table float; e.g., for
##' \code{\\begin\{table\}[htp]}, use \code{floatplace = "htp"}.
##' @param rowsep amount of space (in points) to put between rows.
##' @param useboot logical: use bootstrap estimates (if available) to calculate
##' standard errors?
##' @return \code{x}, invisibly.
##' @export
##' @references Curtis S. Signorino and Ahmer Tarar.  2006.  "A Unified Theory
##' and Test of Extended Immediate Deterrence."  \emph{American Journal of
##' Political Science} 50(3):586--605.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @examples
##' data(war1800)
##' f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | balanc + revis1 | balanc
##' m1 <- egame12(f1, data = war1800)
##'
##' latexTable(m1)
##' latexTable(m1, digits = 8)
##' latexTable(m1, blankfill = "--")  ## dashes in blank cells
##' 
##' \dontrun{
##' latexTable(m1, file = "my_table.tex")  ## write to file}
latexTable <- function(x, digits = max(3, getOption("digits") - 2), scientific =
                       NA, blankfill = "", math.style.negative = TRUE, file =
                       "", floatplace = "htbp", rowsep = 2, useboot = TRUE)
{
    lcat <- function(...) cat(..., file = file, sep = "", append = TRUE)

    n <- names(coef(x))
    cf <- summary(x, useboot = useboot)$coefficients[, 1:2]
    cf <- rbind(cf, c(sum(x$log.likelihood), 0))
    cf <- format(cf, digits = digits, trim = TRUE, scientific = scientific)
    if (math.style.negative)
        cf <- gsub("-", "$-$", cf)

    eqNames <- x$equations[attr(x$equations, "hasColon")]
    varNames <- sapply(sapply(strsplit(n, ":"), '[', -1), paste, collapse = ":")
    varNames <- unique(varNames[nchar(varNames) > 0])
    otherNames <- x$equations[!attr(x$equations, "hasColon") &
                              x$equations %in% n[!x$fixed]]

    lcat("\n%% latex table generated in R ", as.character(getRversion()),
         " by games package\n")
    lcat("%% ", date(), "\n")
    lcat("%% remember to include \\usepackage{multirow} in your preamble\n\n")
    lcat("\\begin{table}[", floatplace, "]\n")
    lcat("\\begin{center}\n")
    lcat("\\begin{tabular}{",
         paste(c("l", rep("c", length(eqNames))), collapse = ""), "}\n")
    lcat("\\hline\n")
    lcat(paste(c("", latexEsc(eqNames)), collapse = " & "), " \\\\\n")
    lcat("\\hline\n")

    for (i in varNames) {
        lcat("\\multirow{2}{*}{", latexEsc(i), "} & ")
        for (J in 1:length(eqNames)) {
            j <- eqNames[J]
            ji <- paste(j, i, sep = ":")
            if (ji %in% n) {
                lcat(cf[ji, 1])
            } else {
                lcat("\\multirow{2}{*}{", blankfill, "}")
            }
            if (J < length(eqNames))
                lcat(" & ")
        }
        lcat(" \\\\\n & ")
        for (J in 1:length(eqNames)) {
            j <- eqNames[J]
            ji <- paste(j, i, sep = ":")
            if (ji %in% n)
                lcat("(", cf[ji, 2], ")")
            if (J < length(eqNames))
                lcat(" & ")
        }
        lcat(" \\\\[", rowsep, "pt]\n")

    }

    if (length(otherNames) > 0) {
        lcat("\\hline\n")
        for (i in otherNames) {
            lcat("\\multirow{2}{*}{", i, "} & ", cf[i, 1], " \\\\\n & (",
                 cf[i, 2], ") \\\\[", rowsep, "pt]\n")
        }
    }

    lcat("\\hline \\hline\n")
    lcat("Log-likelihood & ", cf[nrow(cf), 1], " \\\\\n $N$ & ",
         nrow(x$model), "\\\\\n")
    lcat("\\hline\n")
    
    lcat("\\end{tabular}\n")
    lcat("\\end{center}\n")
    lcat("\\end{table}\n")

    invisible(x)
}
