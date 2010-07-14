\name{strat3}
\alias{strat3}
\title{Fit a strategic model}
\usage{strat3(formulas, data, subset, na.action, varformulas, link=c("probit", "logit"),
    type=c("agent", "private"), startvals=c("sbi", "unif", "zero"), ...)
}
\description{<description>}
\details{<details>}
\value{A \code{strat} object}
\references{Signorino 2003}
\author{Brenton Kenkel (\email{brenton.kenkel@gmail.com})}
\arguments{\item{formulas}{set of formulas describing utilities}
\item{data}{a data frame}
\item{subset}{logical vector, optional, specifying subset of \code{data} to
use}
\item{na.action}{how to deal with \code{NA}s in \code{data}}
\item{varformulas}{set of formulas for regressors in variance, optional}
\item{link}{whether to use probit or logit errors (default probit)}
\item{type}{whether to use an agent-error or private-information stochastic
structure (default agent)}
\item{startvals}{whether to get starting values for the optimization from
statistical backwards induction ("sbi," the default), from a uniform
distribution ("unif"), or to set them all to 0 ("zero")}
\item{...}{other arguments, currently ignored}
}

