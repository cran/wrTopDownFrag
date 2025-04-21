#' Plot Identified Fragments Relative To Their Location
#'
#' Take result from \code{fragIonMass} and display identified fragments by their location, an additional parameter (default logIntensity) is used for coloring
#' This function illustrates  the distribution of identfied peptidesa and thus common break-points.
#' 
#' @param fraL (list) result from \code{fragIonMass}, must conatain elements ''
#' @param extrCol (character) 1st should be aa seq of initial proteins (used for dimensionong graph and separating multiple input proteins), 2nd & 3rd start- and ed-site for drawing;, 4th the column to use for coloring), 5th for protein name in title of figure
#' @param useLog (logical) take values for coloring (4th element of  'extrCol') as log10 
#' @param useCol (character) custom colors
#' @param specLayout (character) custom layout
#' @param useTi (character) custom title
#' @param subTit (character) custom sub-title
#' @param footer (character) custom footer
#' @param batchFig (logical) reduce text content for multiple figues on page
#' @param legCex (numeric, length=1) expansion factor
#' @param legOffS (numeric) legend-offset (passed to \code{\link[wrGraph]{legendHist}})
#' @param legBorder (logical) legend-border (passed to \code{\link[wrGraph]{legendHist}})
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) additional messages for debugging
#' @return This function returns a figure
#' @seealso \code{\link{identifFixedModif}}
#' @examples
#' protP <- c(protP="PEPTIDE")
#' obsMassX <- cbind(a=c(199.1077, 296.1605, 397.2082, 510.2922,625.3192),
#'   b=c(227.1026, 324.1554, 425.2031, 538.2871, 653.3141),
#'   x=c(729.2937, 600.2511, 503.1984, 402.1507, 289.0666),
#'   y=c(703.3145, 574.2719, 477.2191, 376.1714, 263.0874))
#' rownames(obsMassX) <- c("E","P","T","I","D")      # all 1 & 7 ions not included
#' identP1 <- identifFixedModif(prot=protP,expMass=as.numeric(obsMassX), minFragSize=2, 
#'   maxFragSize=7, modTy=list(basMod=c("b","y")))     #
#' 
#' @export
plotFragmLoc <- function(fraL, extrCol=NULL, useLog=FALSE, useCol=NULL, specLayout=NULL, useTi=NULL, subTit=NULL, footer=NULL, batchFig=FALSE, legCex=0.5, legOffS=NULL, legBorder=NULL, silent=FALSE,callFrom=NULL,debug=FALSE){   #pdfName=NULL,fraBest,
  ## plot location of sorted fragments (based on 3rd&4th position of columns specified via 'extrCol' )
  ## 'fraL' .. fragment-list
  ## 'useLog' .. use Abundance data as log10 for color gradient (ie column from 8th element 'extrCol' from in 'fraL')
  ## 'useCol' .. custom color-scheme (ie gradient)
  ## 'specLayout'  .. optional custom layout, if NA automatic layout, if NA previous (or default) layout will be used
  ## 'useTi' .. custom header
  ## 'subTit' .. note at right bottom of image (if null, use name of argument 'fraL')
  ## 'fraBest' should have columns 'ind','rat'
  ## 'fraL' long 'reference' table, should have columns 'beg','end'
  ## version 23mar19  (not implemeted : show spec modif)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotFragmLoc")
  argNa <- deparse(substitute(fraL))
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(is.null(extrCol)) extrCol <- c("orig","beg","end","logInt","origNa")   #not any more "rowNa","orig","ty","modif",modSpec"
  if(is.character(extrCol)) {extrColNa <- extrCol; extrCol <- match(extrCol, colnames(fraL))} else extrColNa <- colnames(fraL)[extrCol]
  if(length(extrCol) < length(extrColNa) || any(is.na(extrCol))) stop(fxNa,"Problem finding columns ",wrMisc::pasteC(extrColNa[which(is.na(extrCol))])) 
  names(extrCol) <- extrColNa
  loc <- fraL[,extrCol]
  
  ## idea? opional filter by 'modif' or 'ty'
  loc <- if(length(unique(loc[,2])) >1) by(loc,loc[,1], as.matrix) else list(loc)             # split by query sequences
  legTx <- ""
  for(i in 1:length(loc)) {    # loop over all proteins tested
    if(debug && i==1) {message(fxNa, "xxPlo1a"); xxPlo1a <- list(fraL=fraL,extrCol=extrCol,i=i,loc=loc,fxNa=fxNa,useLog=useLog,useCol=useCol,specLayout=specLayout,useTi=useTi,silent=silent, extrColNa=extrColNa,batchFig=batchFig)}
    if(debug && i==2) {message(fxNa, "xxPlo1b"); xxPlo1b <- list(fraL=fraL,extrCol=extrCol,i=i,loc=loc,fxNa=fxNa,useLog=useLog,useCol=useCol,specLayout=specLayout,useTi=useTi,silent=silent, extrColNa=extrColNa,batchFig=batchFig)}
    useTit <- if(length(useTi) <1) unique(loc[[i]][,extrCol[5]])[1] else useTi 
    begEnd <- matrix(as.numeric(as.matrix(loc[[i]][,c(2:4)])), ncol=3, dimnames=list(NULL,c("beg","end","intens")))
    reOrd <- order(begEnd[,1], begEnd[,2]-begEnd[,1], decreasing=FALSE)
    begEnd <- begEnd[reOrd,]                               # sort by increasing 'beg'
    loc[[i]] <- loc[[i]][reOrd,]                           # sort by increasing 'beg'
    if(useLog) begEnd[,3] <- log10(begEnd[,3])
    if(is.null(useCol)) {
      useCol1 <- cbind(red=c(141,72,90,171, 220,253,244,255), green=c(129,153,194,221, 216,174,109,0), blue=c(194,203,185,164, 83,97,67,0))   # make gradient    
      useCol1 <- grDevices::rgb(red=useCol1[,1], green=useCol1[,2], blue=useCol1[,3], maxColorValue=255)
    } else useCol1 <- useCol
    if(sum(!is.na(begEnd[,3])) >1) {
      ## make (log) intensity values unique  - this could be done more efficiently !!
      chUni <- duplicated(begEnd[,3])   # need unique values for .cutToNgrp() !!!
    if(debug) {message(fxNa, "xxPlo1d"); xxPlo1d <- list(begEnd=begEnd,fraL=fraL,extrCol=extrCol,i=i,loc=loc,fxNa=fxNa,useLog=useLog,useCol=useCol,specLayout=specLayout,useTi=useTi,silent=silent, extrColNa=extrColNa,batchFig=batchFig)}
      cutInt <- wrMisc::cutToNgrp(begEnd[,3],useCol1)
      lineCol <- if(length(useCol1) >=nrow(fraL)) useCol1 else useCol1[as.integer(cutInt$grouped)]  
    }
    xLab <- if(batchFig) "fragment location" else "aa number"
    if(debug && i==1) {message(fxNa, "xxPlo2a"); xxPlo2a <- list(begEnd=begEnd,legTx=legTx,loc=loc,fraL=fraL,useCol=useCol,useCol1=useCol1,lineCol=lineCol,i=i,extrCol=extrCol,extrColNa=extrColNa,useTi=useTi,useLog=useLog,specLayout=specLayout,useTi=useTi,useCol1=useCol1,silent=silent)} # intGrp=intGrp,ra=ra,fragLe=fragLe,intens=intens,
    if(debug && i==2) {message(fxNa, "xxPlo2b"); xxPlo2b <- list(begEnd=begEnd,legTx=legTx,loc=loc,fraL=fraL,useCol=useCol,useCol1=useCol1,lineCol=lineCol,i=i,extrCol=extrCol,extrColNa=extrColNa,useTi=useTi,useLog=useLog,specLayout=specLayout,useTi=useTi,useCol1=useCol1,silent=silent)} # intGrp=intGrp,ra=ra,fragLe=fragLe,intens=intens,
    ## main plot
    graphics::plot(c(1,nchar(as.character(loc[[i]][1,1]))), c(1,nrow(loc[[i]])), type="n", main=useTit, xlab="aa-position", ylab="peptide number", las=1)
    if(nrow(loc[[i]]) >2) graphics::segments(0,1,nchar(as.character(loc[[i]][1,1])), nrow(loc[[i]]), lty=2, col=grDevices::grey(0.6))                      # diag increasing
    subTit2 <- if(batchFig) paste("colored by",extrColNa[4]) else paste("fragments colored by",extrColNa[4]," (blue for weak, red for high)")
    graphics::mtext(subTit2, cex=0.7 -0.1*batchFig)
    lWidth <- sort(c(0.2, round(10/sqrt(nrow(begEnd)),2),2))[2]
    graphics::segments(x0=begEnd[,1], y0=1:nrow(loc[[i]]), x1=begEnd[,2], y1=1:nrow(loc[[i]]), col=lineCol, lwd=lWidth)           #segments(x0, y0, x1=x0, y1=y0,  
    wrGraph::legendHist(begEnd[,3], useCol1, cex=legCex, offS=legOffS, border=legBorder)     
    }}
 
