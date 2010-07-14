\name{strat3}
\alias{strat3}
\title{Fit a strategic model with 3 terminal nodes}
\usage{strat3(formulas, data, subset, na.action, varformulas, link=c("probit",
    "logit"), type=c("agent", "private"), startvals=c("sbi", "unif",
    "zero"), ...)
}
\description{Fits a strategic model with three terminal nodes, as in the game illustrated
below in \dQuote{Details}.}
\details{The model corresponds to the following extensive-form game, described in
Signorino (2003):
\preformatted{
.     /\
.    /  \
.   /    \
.  u11   /\
.       /  \
.      /    \
.    u13    u14
.    0      u24}

If Player 1 chooses L, the game ends and Player 1 receives payoffs of u11.
(Player 2's utilities in this case cannot be identified in a statistical
model.)  If Player 1 chooses L, then Player 2 can choose L, resulting in
payoffs of u13 for Player 1 and 0 for Player 2, or R, with payoffs of u14
for 1 and u24 for 2.

The four formulas specified in the function's \code{formulas} argument
correspond to the regressors to be placed in u11, u13, u14, and u24
respectively.  If there is any regressor (including the constant) placed in
all of u11, u13, and u14, \code{strat3} will stop and issue an error
message, because the model is then unidentified (see Lewis and Schultz
2003).  There are two equivalent ways to express the formulas passed to this
argument.  One is to use a list of four formulas, where the first contains
the response variable(s) (discussed below) on the left-hand side and the
other three are one-sided.  For instance, suppose u11 is a function of
\code{x1}, \code{x2}, and a constant; u13 is set to 0; u14 is a function of
\code{x3} and a constant; and u24 is a function of \code{z} and a constant.
The list notation would be \code{formulas = list(y ~ x1 + x2, ~ 0, ~ x3, ~
z)}.  The other method is to use the \code{\link{Formula}} syntax, with one
left-hand side and four right-hand sides (separated by vertical bars).  This
notation would be \code{formulas = y ~ x1 + x2 | 0 | x3 | z}.

There are three equivalent ways to specify the outcome in \code{formulas}.
One is to use a numeric vector with three unique values, with their values
(from lowest to highest) corresponding with the terminal nodes of the game
tree illustrated above (from left to right).  The second is to use a factor,
with the levels (in order as given by \code{levels(y)}) corresponding to the
terminal nodes.  The final way is to use two indicator variables, with the
first standing for whether Player 1 moves L (0) or R (1), the second
standing for Player 2's choice if Player 1 moves R.  (The values of the
second when Player 1 moves L should be set to 0 or 1, \strong{not}
\code{NA}, in order to ensure that observations are not dropped from the
data when \code{na.action = na.omit}.)  The way to specify \code{formulas}
when using indicators variables is, for example, \code{y1 + y2 ~ x1 + x2 | 0
| x3 | z}.}
\value{An object of class \code{c("strat", "strat3")}. A
\code{strat} object is a list containing: \describe{
\item{\code{coefficients}}{estimated parameters of the model}
\item{\code{vcov}}{estimated variance-covariance matrix}
\item{\code{log.likelihood}}{vector of individual log likelihoods (left unsummed
for use with non-nested model tests)}
\item{\code{call}}{the call used to produce the model}
\item{\code{formulas}}{the final \code{Formula} object passed to
\code{model.frame}}
\item{\code{link}}{the link function used}
\item{\code{type}}{the stochastic structure used}
\item{\code{model}}{the model frame containing all variables used in fitting}
\item{\code{y}}{the dependent variable, represented as a factor}
}
The second class, \code{strat3}, is for use with the \code{predict} method.}
\seealso{\code{\link{summary.strat}}, \code{\link{Formula}}}
\references{Jeffrey B. Lewis and Kenneth A Schultz.  2003.
\dQuote{Revealing Preferences: Empirical Estimation of a Crisis Bargaining
Game with Incomplete Information.}  \emph{Political Analysis} 11:345--367.

Curtis S. Signorino.  2003.  \dQuote{Structure and Uncertainty
in Discrete Choice Models.}  \emph{Political Analysis} 11:316--344.}
\author{Brenton Kenkel (\email{brenton.kenkel@gmail.com})}
\arguments{\item{formulas}{a list of four formulas, or a \code{Formula} object with
four right-hand sides.  See \dQuote{Details} and \dQuote{Examples}.}
\item{data}{a data frame}
\item{subset}{an optional logical vector specifying which observations from
\code{data} to use in fitting}
\item{na.action}{how to deal with \code{NA}s in \code{data}.  Defaults to
the \code{na.action} setting of \code{\link{options}}.  See
?\code{\link{na.omit}}}
\item{varformulas}{an optional list of four formulas or \code{Formula}
object with four right-hand sides.  If omitted, all error terms are assumed
to have scale parameter 1.  See \dQuote{Details}.}
\item{link}{whether to use a probit (default) or logit link structure}
\item{type}{whether to use an agent-error (\dQuote{agent}, default) or
private-information (\dQuote{private}) stochastic structure}
\item{startvals}{whether to get starting values for the optimization from
statistical backwards induction (\dQuote{sbi}, default), from a uniform
distribution (\dQuote{unif}), or to set them all to 0 (\dQuote{zero})}
\item{...}{other arguments to pass to the fitting function}
}

