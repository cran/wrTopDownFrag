% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreChargeCatch.R
\name{.chargeCatchingAA}
\alias{.chargeCatchingAA}
\title{Cite Charge Catching Amino-Acids}
\usage{
.chargeCatchingAA(
  chargeMode = "pos",
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{chargeMode}{(character) this value may be 'pos' (default) for the positively charged amino-acids K,R and H or, 
if this argument has any other value, than all charged amino-acids (K,R,H, S,T,N,Q, D,E, W and Y) will be considered.}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages and objects exportet to current session for debugging}

\item{callFrom}{(character) allow easier tracking of messages produced}
}
\value{
This function returns a matrix with charge-catching amino-acids and their assumed strength
}
\description{
Return a matrix with charge-catching amino-acids and their assumed strength.
So far, the strength shown/used is set rather empirically.
}
\examples{
.chargeCatchingAA()
}
\seealso{
\code{\link{fragmentSeq}}
}
