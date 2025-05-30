% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPrefFragPat.R
\name{plotPrefFragPat}
\alias{plotPrefFragPat}
\title{plot preferential fragmenation pattern
  
Plot preferential fragmenation pattern equivalent to Fig 1b of Haverland et al 2017 (J Am Soc Mass Spectrom)}
\usage{
plotPrefFragPat(
  prefPat,
  namesCex = 0.8,
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{prefPat}{(matix)}

\item{namesCex}{(numeric) expansion factor cex for display of AA-names}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages and objects exportet to current session for debugging}

\item{callFrom}{(character) allow easier tracking of messages produced}
}
\value{
This function returns a figure
}
\description{
plot preferential fragmenation pattern
  
Plot preferential fragmenation pattern equivalent to Fig 1b of Haverland et al 2017 (J Am Soc Mass Spectrom)
}
\examples{
plotPrefFragPat(.prefFragPattern())
}
\seealso{
\code{\link{scoreFragments}}
}
