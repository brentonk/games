\docType{data}
\name{leblang2003}
\alias{leblang2003}
\title{Currency attacks}
\usage{
data(leblang2003)
}
\description{
Data on speculative currency attacks and devaluation decisions for 90
countries from 1985 to 1998.
}
\details{
The dataset is taken from Leblang (2003).  The unit of
observation is the country-month, and the variables are: \describe{
\item{\code{outcome}}{Whether the country faced no speculative attack,
defended its currency against an attack, or devalued in response to an
attack in the given month}
\item{\code{preelec}}{Indicator for being in the three months prior to an
election}
\item{\code{postelec}}{Indicator for being in the three months following an
election}
\item{\code{rightgov}}{Indicator for a right-wing party being in power}
\item{\code{unifgov}}{Indicator for unified government: in presidential
systems, the same party controlling the presidency and the lower house of
the legislature; in parliamentary systems, one party/coalition having a
majority of seats}
\item{\code{lreserves}}{Logged ratio of currency reserves to base money in
the previous month}
\item{\code{realinterest}}{Domestic real interest rate in the previous
month}
\item{\code{lexports}}{Logged ratio of exports to GDP in the previous month}
\item{\code{capcont}}{Indicator for capital controls in the previous year}
\item{\code{overval}}{Overvaluation of the real exchange rate}
\item{\code{creditgrow}}{Domestic credit growth in the previous month}
\item{\code{service}}{External debt service (as percentage of exports) paid
in previous month}
\item{\code{USinterest}}{U.S. domestic interest rates in the previous month}
\item{\code{contagion}}{Number of other countries experiencing speculative
attacks in the same month}
\item{\code{prioratt}}{Number of prior speculative attacks experienced by the
country}
\item{\code{nation}}{Country name}
\item{\code{month}}{Month of observation}
\item{\code{year}}{Year of observation}
}
All of the non-binary variables other than \code{nation}, \code{month},
and \code{year} are standardized to have mean 0 and unit variance.

We are grateful to David Leblang for allowing us to redistribute his data.
The original replication file is available in Stata format at
\url{http://www.faculty.virginia.edu/leblang/Leblang/Data_files/strat.dta}
(as of 2010-03-01).
}
\examples{
## replicate analysis in Leblang (2003)
data(leblang2003)
m1 <- egame12(outcome ~ capcont + lreserves + overval + creditgrow +
USinterest + service + contagion + prioratt - 1 | 1 | 1 | unifgov + lexports
+ preelec + postelec + rightgov + realinterest + capcont + lreserves, data =
leblang2003, link = "probit", type = "private")

summary(m1)
}
\references{
David Leblang.  2003.  "To Defend or Devalue: The Political
Economy of Exchange Rate Policy."  \emph{International Studies Quarterly}
47: 533--559.
}
\seealso{
\code{\link{egame12}}
}
\keyword{data}

