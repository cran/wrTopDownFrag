% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identifFixedModif.R
\name{.singleSpecModif}
\alias{.singleSpecModif}
\title{Add Single Specific Modifications}
\usage{
.singleSpecModif(
  pepTab,
  specModif,
  nMaxMod = 1,
  massTy = "mono",
  callFrom = NULL,
  silent = FALSE,
  debug = FALSE
)
}
\arguments{
\item{pepTab}{(matrix) matrix of fragments (cols 'no','seq','orig','ty','seqNa','beg','end','precAA','tailAA','ambig','mass')}

\item{specModif}{(list) with elements 'modOrigin' (sequence), 'modPos' (position within sequence), 'modMass' (digits, ie mass to add),
'modName' (name of modif), 'modFixed' (fixed or , logical)}

\item{nMaxMod}{(numeric) max number a given modification may occur}

\item{massTy}{(character) 'mono' or 'average'}

\item{callFrom}{(character) allow easier tracking of message(s) produced}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages and objects exportet to current session for debugging}
}
\value{
This function returns a list with $massMatch (list of exerimental peptides matching to one or more predicted), $preMa (predicted ions, including fixed modif), $pepTab (predicted neutral peptides, wo modifications), $expMa (experimental mass from input), $recalibFact (recalibration factor as from input), $docTi (time for calculations)
}
\description{
Add single specific modification to peptide/protein fragments .
}
\examples{
pep1 <- c(pe1="KPEPTI")
# The table of possible terminal fragments (for simplicity terminal only)
pepTab1 <- makeFragments(pep1, min=3, max=7, internFra=FALSE)
specModif1 <- list(modOrigin=pep1, modPos=1, modMass=579.9663, modName="p", modFixed=FALSE)
.singleSpecModif(pepTab1, specModif1 )

protP <- c(protP="PEPTIDEKR")
pep1 <- c("PTI","KPE","EPTI")
papTab1 <- cbind(no=c(7,2,6),seq=pep1, orig=rep("KPEPTI",3), origNa=rep("pe1",3), 
  ty=paste0(c("C","N","C"),"ter"), beg=c(4,1,3), end=c(6,3,4),
  mass= wrProteo::convAASeq2mass(pep1, massTy="mono"), modSpec="")

}
\seealso{
\code{\link{makeFragments}}, \code{\link{identifVarModif}}, \code{\link{identifyPepFragments}}
}
