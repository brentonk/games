##' @include games.r
##' @include helpers.r
NULL

makeFormulas <- function(model, yvals, yname = "y")
{
    ## convert 'model' to character if supplied without quotation marks
    ischar <- tryCatch(is.character(model) && length(model) == TRUE, error =
                       identity)
    if (inherits(ischar, "error"))
        ischar <- FALSE
    if (!ischar)
        model <- deparse(substitute(model))

    ## extract outcome names from 'yvals' if not directly supplied
    if (!is.character(yvals))
        yvals <- levels(as.factor(yvals))

    ## which players can have utilities in which equations
    eqByPlayer <- switch(model,
                         egame12 = list(rep(TRUE, 3), c(FALSE, FALSE, TRUE)),
                         egame122 = list(rep(TRUE, 4), c(FALSE, TRUE, FALSE,
                         TRUE)),
                         egame123 = list(rep(TRUE, 4), c(FALSE, rep(TRUE, 3)),
                         c(rep(FALSE, 3), TRUE)),
                         ultimatum = list(c(TRUE, FALSE), c(FALSE, TRUE)))
    if (is.null(eqByPlayer))
        stop("'", model, "' is not a model in the 'games' package; see 'help(package=\"games\")'")

    ## set of equations relevant to identification for a player's utilities
    ## (i.e. no equations for actions before their move)
    idByPlayer <- switch(model,
                         egame12 = list(1:3, 2:3),
                         egame122 = list(1:4, 1:4),
                         egame123 = list(1:4, 2:4, 3:4),
                         ultimatum = list(1:2, 1:2))  # not actually relevant
                                                      # for ultimatum

    nplayers <- length(eqByPlayer)
    neqs <- length(eqByPlayer[[1]])
    ans <- vector("list", nplayers)

    choices <- c("fix to 0", "intercept only", "regressors, no intercept",
                 "regressors with intercept")

    for (i in seq_len(nplayers)) {
        ## storing selections so options with intercepts can be removed when
        ## necessary for identification.  the last is set to 2 to avoid making
        ## 'allowCond' (below) even more convoluted
        sel <- c(rep(1, neqs - 1), 2)
        noconstant <- rep(FALSE, neqs)
        ans[[i]] <- vector("list", neqs)
        
        for (j in seq_len(neqs)) {
            if (!eqByPlayer[[i]][j])
                next

            if (model == "ultimatum") {
                title <- paste("\nPlayer ", i, "'s reservation value:", sep = "")
            } else {
                title <- paste("\nEquation for player ", i, "'s utility from ",
                               yvals[j], ":", sep = "")
            }
            allowCond <- j < neqs || !all(sel[idByPlayer[[i]]] %in% c(2, 4))
            allowed <- if (allowCond) 1:4 else c(1,3)
            eqtype <- menu(choices[allowed], title = title)

            if (eqtype == 0L) {
                stop("stopped by user")
            } else if (eqtype == 1L) {
                ans[[i]][[j]] <- "0"
            } else if (eqtype == 2L && allowCond) {
                ans[[i]][[j]] <- "1"
            } else {
                if (eqtype == 3L || !allowCond)
                    noconstant[j] <- TRUE
                
                repeat {
                    varnames <-
                        readline("\nEnter variable names (separated by spaces):\n")
                    varnames <- str_split(varnames, " ")[[1]]
                    varnames <- varnames[varnames != ""]
                    varnames <- varnames[!duplicated(varnames)]
                    ans[[i]][[j]] <- varnames
                    bad <- do.call(intersectAll, ans[[i]][idByPlayer[[i]]])
                    if (j == neqs && length(bad)) {
                        cat("The following variables cannot be used due to identification problems: ",
                            paste(bad, collapse = ", "), ". Try again.\n", sep =
                            "")
                    } else {
                        break
                    }
                }
            }

            sel[j] <- eqtype
        }

        ## need to close the old loop and start a new one before combining
        ## regressors into equation form, else the identification check won't
        ## work
        for (j in seq_len(neqs)) {
            if (length(ans[[i]][[j]]) >= 1) {
                ans[[i]][[j]] <- paste(ans[[i]][[j]], collapse = " + ")
                if (noconstant[j])
                    ans[[i]][[j]] <- paste(ans[[i]][[j]], "- 1")
            }
        }
    }

    ans <- unlist(ans)
    ans <- paste(ans, collapse = " | ")
    ans <- paste(yname, "~", ans)
    ans <- as.Formula(ans)

    return(ans)
}
