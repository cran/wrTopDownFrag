% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AAfragSettings.R
\name{AAfragSettings}
\alias{AAfragSettings}
\title{Settings For AA Fragmentation}
\usage{
AAfragSettings(outTy = "all", silent = FALSE, debug = FALSE, callFrom = NULL)
}
\arguments{
\item{outTy}{(character) default "all" or any of the list-elements}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages for debugging}

\item{callFrom}{(character) allows easier tracking of messages produced}
}
\value{
This function returns a list ($knownMods, $knspecAAMods, $modChem, $neutralLossOrGain)
}
\description{
This function provides basic settings for what types of fragments may accomodate which type of modifications :  $knownMods: information about which modifications may be considered, $specAAMod: specifc AA sites (if applicable),  $specAAMod: specifc AA sites (if applicable).
For example, here 'p' codes for gain of mass for HPO3 only at S, T and Y residues.
Note: $knownMods$Nterm and $knownMods$Cterm are treated as mutually exclusive
}
\examples{
AAfragSettings()
}
\seealso{
\code{\link{makeFragments}},  \code{\link{fragmentSeq}}, \code{\link[wrProteo]{massDeFormula}}
}
