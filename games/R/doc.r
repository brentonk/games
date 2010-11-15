##' @include games.r
NULL

##' Dataset of militarized international disputes between 1816 and 1899.
##'
##' The dataset is taken from the Correlates of War project.  The unit of
##' observation is the dyad-year, and the variables are: \describe{
##' \item{\code{ccode1}}{Initiator's COW country code}
##' \item{\code{ccode2}}{Respondent's COW country code}
##' \item{\code{year}}{Year of dispute}
##' \item{\code{cap_1}}{Initiator's military capabilities (as percent of total
##' system capabilities)}
##' \item{\code{cap_2}}{Respondent's military capabilities (as percent of total
##' system capabilities)}
##' \item{\code{balanc}}{Balance of dyadic capabilities possessed by the
##' initiator (i.e., \code{cap_1 / (cap_1 + cap_2)})}
##' \item{\code{s_wt_re1}}{Dyadic S-score (see Signorino and Ritter 1998),
##' weighted by initiator's region}
##' \item{\code{s_wt_re2}}{Dyadic S-score, weighted by respondent's region}
##' \item{\code{dem1}}{Initiator's Polity score}
##' \item{\code{dem2}}{Respondent's Polity score}
##' \item{\code{distance}}{Distance (in miles) between initiator and respondent}
##' \item{\code{peaceyrs}}{Years since last dispute in this dyad}
##' \item{\code{midnum}}{Dispute's number in the MID data set}
##' \item{\code{revis1}}{Whether the initiator had \dQuote{revisionist} aims}
##' \item{\code{revis2}}{Whether the respondent had \dQuote{revisionist} aims}
##' \item{\code{sq}}{Indicator for status quo outcome}
##' \item{\code{capit}}{Indicator for capitulation outcome}
##' \item{\code{war}}{Indicator for war outcome}
##' \item{\code{esc}}{Indicator for escalation (i.e., either capitulation or war
##' occurs)}
##' \item{\code{regime1}}{Initiator's regime type (calculated from \code{dem1})}
##' \item{\code{regime2}}{Respondent's regime type (calculated from \code{dem2})}
##' }
##' @name war1800
##' @usage data(war1800)
##' @title 19th-century international disputes
##' @docType data
##' @references Daniel M. Jones, Stuart A. Bremer and J. David Singer.  1996.
##' \dQuote{Militarized Interstate Disputes, 1816-1992: Rationale, Coding Rules,
##' and Empirical Patterns.} \emph{Conflict Management and Peace Science}
##' 15(2): 163--213.
##' @seealso \code{\link{egame12}}
##' @keywords data
NULL

##' Simulated data for illustrating \code{\link{egame122}}.
##'
##' The variables are: \describe{
##' \item{\code{f1}, \code{f2}}{Factors with levels \dQuote{a}, \dQuote{b},
##' \dQuote{c}}
##' \item{\code{x1}--\code{x5}}{Numeric variables entering Player 1's utilities}
##' \item{\code{z1}--\code{z3}}{Numeric variables entering Player 2's utilities}
##' \item{\code{a1}}{Indicator for Player 1's move (L or R)}
##' \item{\code{a2}}{Indicator for Player 2's move (L or R)}
##' \item{\code{y}}{Factor containing outcome}
##' }
##' @name data_122
##' @usage data(data_122)
##' @title Simulated egame122 data
##' @docType data
##' @seealso \code{\link{egame122}}
##' @keywords data
NULL

##' Simulated data for illustrating \code{\link{ultimatum}}.
##'
##' The variables are: \describe{
##' \item{\code{offer}}{The offer made by Player 1}
##' \item{\code{accept}}{Whether Player 2 accepted the offer (0 for rejection, 1
##' for acceptance)}
##' \item{\code{w1}, \code{w2}}{Variables entering both players' reservation values}
##' \item{\code{x1}--\code{x4}}{Variables entering Player 1's reservation value}
##' \item{\code{z1}--\code{z4}}{Variables entering Player 2's reservation value}
##' }
##' The maximum offer size is 15.
##' @name data_ult
##' @usage data(data_ult)
##' @title Simulated ultimatum data
##' @docType data
##' @seealso \code{\link{ultimatum}}
##' @keywords data
NULL

##' Data from a trial of the ultimatum game with graduate students.
##'
##' The variables are: \describe{
##' \item{\code{offer}}{The offer made by Player 1}
##' \item{\code{accept}}{Whether Player 2 accepted the offer (0 for rejection, 1
##' for acceptance)}
##' \item{\code{gender1}}{Whether Player 1 is female}
##' \item{\code{gender2}}{Whether Player 2 is female}
##' }
##' The maximum offer size is 100.
##' @name student_offers
##' @usage data(student_offers)
##' @title Data from students playing the ultimatum game
##' @docType data
##' @seealso \code{\link{ultimatum}}
##' @keywords data
NULL
