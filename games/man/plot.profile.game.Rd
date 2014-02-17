\name{plot.profile.game}
\alias{plot.profile.game}
\title{Plot profiles of strategic model log-likelihoods}
\usage{
\method{plot}{profile.game}(x, show.pts = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{profile.game},
  typically created by running \code{\link{profile.game}}
  on a fitted \code{game} model}

  \item{show.pts}{logical: plot a point for the
  log-likelihood of each profiled model?}

  \item{...}{other arguments, currently ignored}
}
\value{
\code{x}, invisibly
}
\description{
Plot output of \code{\link{profile.game}}.
}
\details{
This method provides plots for a quick assessment of whether \code{game}
models have failed to converge to a global maximum.  For each parameter
profiled (see \code{\link{profile.game}} for details of the profiling
process), a spline interpolation of the log-likelihood profile is provided,
with an "x" marking the value at the original parameter estimate.

Sometimes the plot will seem to indicate that the original fit did not reach
the global maximum, even though \code{profile.game} did not issue the
non-covergence warning.  This is an artifact of the interpolation, which can
be confirmed by re-running \code{plot.profile.game} with the argument
\code{show.pts = TRUE}.
}
\author{
Brenton Kenkel
}
\seealso{
\code{\link{profile.game}}
}

