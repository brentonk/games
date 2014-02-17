\name{predict.game}
\alias{predict.egame12}
\alias{predict.egame122}
\alias{predict.egame123}
\alias{predict.game}
\alias{predict.ultimatum}
\title{Predicted probabilities for strategic models}
\usage{
\method{predict}{game}(object, ...)

\method{predict}{egame12}(object, newdata, type=c("outcome", "action"), na.action = na.pass, ...)
\method{predict}{egame122}(object, newdata, type=c("outcome", "action"), na.action = na.pass, ...)
\method{predict}{egame123}(object, newdata, type=c("outcome", "action"), na.action = na.pass, ...)
\method{predict}{ultimatum}(object, newdata, na.action = na.pass, n.sim = 1000, ...)
}
\arguments{
  \item{object}{a fitted model of class \code{game}.}

  \item{newdata}{data frame of values to make the predicted
  probabilities for. If this is left empty, the original
  dataset is used.}

  \item{type}{whether to provide probabilities for outcomes
  (e.g., L, RL, or RR in \code{egame12}) or for actions
  (e.g., whether 2 moves L or R given that 1 moved R).}

  \item{na.action}{how to deal with \code{NA}s in
  \code{newdata}}

  \item{n.sim}{number of simulation draws to use per
  observation for \code{ultimatum} models (see Details).}

  \item{...}{other arguments, currently ignored.}
}
\value{
A data frame of predicted probabilities.
}
\description{
Makes predicted probabilities from a strategic model.
}
\details{
This method uses a fitted strategic model to make predictions for a new
set of data.  This is useful for cross-validating or for graphical
analysis.  For many uses, such as analyzing the marginal effect of a
particular independent variable, the function \code{\link{predProbs}} will
be more convenient.

In the \code{\link{ultimatum}} model, there is not an analytic expression
for the expected value of Player 1's offer.  Therefore, predicted values
are instead generating via simulation by drawing errors from a logistic
distribution.  The number of draws per observation can be controlled via
the \code{n.sim} argument.  For replicability, we recommend seeding the
random number generator via \code{\link{set.seed}} before using
\code{predict.ultimatum}.
}
\author{
Brenton Kenkel (\email{brenton.kenkel@gmail.com})
}
\seealso{
\code{\link{predProbs}} provides a more full-featured and
user-friendly wrapper, including plots and confidence bands.
}

