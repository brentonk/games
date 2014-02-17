\name{vuong}
\alias{clarke}
\alias{vuong}
\title{Non-nested model tests}
\usage{
vuong(model1, model2, outcome1=NULL, outcome2=NULL, level=0.05,
digits=2)

clarke(model1, model2, outcome1=NULL, outcome2=NULL, level=0.05, digits=2)
}
\arguments{
  \item{model1}{A fitted statistical model of class
  \code{"game"}, \code{"lm"}, or \code{"glm"}}

  \item{model2}{A fitted statistical model of class
  \code{"game"}, \code{"lm"}, or \code{"glm"} whose
  dependent variable is the same as that of \code{model1}}

  \item{outcome1}{Optional: if \code{model1} is of class
  \code{"game"}, specify an integer to restrict attention
  to a particular binary outcome (the corresponding column
  of \code{predict(model1)}).  For \code{\link{ultimatum}}
  models, "offer" or "accept" may also be used.  See
  "Details" below for more information on when to specify
  an outcome.  If \code{model1} is not of class
  \code{"game"} and \code{outcome1} is non-\code{NULL}, it
  will be ignored and a warning will be issued.}

  \item{outcome2}{Optional: same as \code{outcome1}, but
  corresponding to \code{model2}.}

  \item{level}{Numeric: significance level for the test.}

  \item{digits}{Integer: number of digits to print}
}
\value{
Typical use will be to run the function interactively and examine
the printed output.  The functions return an object of class
\code{"nonnest.test"}, which is a list containing: \describe{
\item{\code{stat}}{The test statistic}
\item{\code{test}}{The type of test (\code{"vuong"} or \code{"clarke"})}
\item{\code{level}}{Significance level for the test}
\item{\code{digits}}{Number of digits to print}
\item{\code{loglik1}}{Vector of observationwise log-likelihoods for
\code{model1}}
\item{\code{loglik2}}{Vector of observationwise log-likelihoods for
\code{model2}}
\item{\code{nparams}}{Integer vector containing the number of parameters
fitted in \code{model1} and \code{model2} respectively}
\item{\code{nobs}}{Number of observations of the dependent variable being
modeled}}
}
\description{
Perform Vuong's (1989) or Clarke's (2007) test for non-nested model
selection.
}
\details{
These tests are for comparing two statistical models that have the same
dependent variable, where neither model can be expressed as a special case
of the other (i.e., they are non-nested).  The null hypothesis is that the
estimated models are the same Kullback-Leibler distance from the true
model.  To adjust for potential differences in the dimensionality of the
models, the test statistic for both \code{vuong} and \code{clarke} is
corrected using the Bayesian information criterion (see Clarke 2007 for
details).

It is crucial that the dependent variable be exactly the same between the
two models being tested, including the order the observations are placed in.
The \code{vuong} and \code{clarke} functions check for such discrepancies,
and stop with an error if any is found.  Models with non-null \code{weights}
are not yet supported.

When comparing a strategic model to a (generalized) linear model, you must
take care to ensure that the dependent variable is truly the same between
models.  This is where the \code{outcome} arguments come into play.  For
example, in an \code{\link{ultimatum}} model where acceptance is observed,
the dependent variable for each observation is the vector consisting of the
offer size and an indicator for whether it was accepted.  This is not the
same as the dependent variable in a least-squares regression of offer size,
which is a scalar for each observation.  Therefore, for a proper comparison
of \code{model1} of class {"ultimatum"} and \code{model2} of class
\code{"lm"}, it is necessary to specify \code{outcome1 = "offer"}.
Similarly, consider an \code{\link{egame12}} model on the
\code{\link{war1800}} data, where player 1 chooses whether to escalate the
crisis and player 2 chooses whether to go to war.  The dependent variable
for each observation in this model is the vector of each player's choice.
By contrast, in a logistic regression where the dependent variable is
whether war occurs, the dependent variable for each observation is a
scalar.  To compare these models, it is necessary to specify \code{outcome1
= 3}.
}
\examples{
data(war1800)

## balance of power model
f1 <- esc + war ~ balanc + s_wt_re1 | 0 | balanc | balanc + s_wt_re1
m1 <- egame12(f1, data = war1800, subset = !is.na(regime1) & !is.na(regime2))

## regime type model
f2 <- esc + war ~ regime1 | 0 | regime1 + regime2 | regime1 + regime2
m2 <- egame12(f2, data = war1800)

## comparing two strategic models
vuong(model1 = m1, model2 = m2)
clarke(model1 = m1, model2 = m2)

## comparing strategic model to logit, must specify outcome1 appropriately
logit1 <- glm(war ~ balanc + s_wt_re1, data = m1$model, family=binomial)
vuong(model1 = m1, outcome1 = 3, model2 = logit1)
clarke(model1 = m1, outcome1 = 3, model2 = logit1)

logit2 <- glm(sq ~ regime1 + regime2, data = war1800, family=binomial)
vuong(model1 = m2, outcome1 = 1, model2 = logit2)
clarke(model1 = m2, outcome1 = 1, model2 = logit2)

## ultimatum model
data(data_ult)
f3 <- offer + accept ~ w1 + w2 + x1 + x2 | w1 + w2 + z1 + z2
m3 <- ultimatum(f3, maxOffer = 15, data = data_ult)
ols1 <- lm(offer ~ w1 + w2 + x1 + x2 + z1 + z2, data = data_ult)
vuong(model1 = m3, outcome1 = "offer", model2 = ols1)
clarke(model1 = m3, outcome1 = "offer", model2 = ols1)
}
\author{
Brenton Kenkel (\email{brenton.kenkel@gmail.com})
}
\references{
Quang H. Vuong.  1989.  "Likelihood Ratio Tests for Model
Selection and Non-Nested Hypotheses."  \emph{Econometrica} 57(2): 307--333.

Kevin Clarke.  2007.  "A Simple Distribution-Free Test for Nonnested
Hypotheses."  \emph{Political Analysis} 15(3): 347--363.
}

