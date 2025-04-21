#' Draw simplified (deconvoluted) spectrum of mgf type and highlight peaks with matches found to theoretical data
#'
#' Draw simplified (deconvoluted) spectrum of mgf type and highlight matches found
#'
#' @param lst (list) result of identificationi with element \code{$obsMass} (the column 'obsMass' should be m/z values, the column 'sc.logInt' intensity values) and $identif
#' @param basInp (numeric) alternative/custom entry of observaed masses
#' @param backgrCol (character) color of background
#' @param replNames (character) replace terms when showing names : the first character-string will be replaced by the second character-string
#' @param lwd (numeric) line width
#' @param col (character) custom colors for different types of ions/modifications
#' @param tit (character) custom title
#' @param xLim (numeric length=2) custom x axis margins
#' @param yLab (character) custom y axis label
#' @param replNames (character)
#' @param linPlot (logical) re-transform y-axis from log2 to linear scale 
#' @param listNa (character) list-elements of 'lst' to use/extract
#' @param useColN (character) columns names tu use from input ie lst$identif & lst$overview
#' @param cex (numeric) expansion factor for x- and y-label
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a mgf-like figure
#' @seealso \code{\link{makeFragments}}, \code{\link{identifVarModif}}, \code{\link{identifFixedModif}}, \code{\link{identifyPepFragments}}
#' @examples
#' set.seed(2025)
#' @export
plotMgfLike <- function(lst, basInp=NULL, backgrCol=NULL, replNames=c("by","i"), lwd=1, col=NULL, tit=NULL, xLim=NULL, yLab=NULL, linPlot=FALSE,
  listNa=c("identif","overview","obsMass"), useColN=c("obsMass","logInt","origNa","mod","modSpec","fraNa","orig", "nIdentif","minMassDec","maxMass"), 
  cex=1, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## 'col' .. custom colors for different types of ions/modifications 
  ## idea: indicate full length prot 'precursor' MH+ (at fixed modif)
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotMgfLike")
  intensDataIsLin <- TRUE
  chNa <- listNa %in% names(lst)
  if(any(!chNa)) stop(fxNa,"Cannot find fields ",wrMisc::pasteC(listNa,quoteC="'",lastCol=" and/or ")," in input !")
  chNa <- useColN[1:5] %in% colnames(lst[[listNa[1]]])
  if(any(!chNa)) stop(fxNa,"Cannot find columns ",wrMisc::pasteC(useColN[1:5],quoteC="'",lastCol=" and/or ")," in input !")
   if(is.null(backgrCol))  backgrCol <- grDevices::grey(0.85)
  x <- if(is.null(basInp)) lst[[listNa[3]]] else basInp                                                                 # orig measures : mass & logInt
  x2 <- lst[[listNa[1]]][,useColN[1:2]]                                       # identified pep (mass & logInt)
  if(!is.matrix(x2)) x2 <- matrix(as.numeric(as.matrix(x2)), ncol=2, dimnames=dimnames(x2))   # conv data.fr -> matrix
  x3 <- as.matrix(lst[[listNa[1]]][,useColN[3:6]])                            # supl info from identif
  if(length(x) >0 && length(x2) >0) {
    if(debug) {message("..xx1a"); xx1a <- list(x=x,x2=x2,x3=x3,lst=lst,basInp=basInp,listNa=listNa,useColN=useColN,replNames=replNames,backgrCol=backgrCol,lwd=lwd,col=col,intensDataIsLin=intensDataIsLin,fxNa=fxNa)}
    if(!is.matrix(x)) x <- matrix(as.numeric(as.matrix(x[,1:2])), ncol=2, dimnames=dimnames(x))   # conv data.fr -> matrix
    matchIr <- wrMisc::naOmit(match(x2[,1],x[,1]))                               # try matching measures m/z
    matchI <- wrMisc::naOmit(match(x[,1],x2[,1]))                               # try matching measures m/z
    if(length(matchI) <1) message(fxNa,"Problem fetching (",nrow(lst[[listNa[1]]]),") identified peaks -none left !!")
    x2 <- if(length(matchI) >1) x2[matchI,] else matrix(x2[matchI,], ncol=2, dimnames=list(rownames(x2)[matchI],colnames(x2)))
    x3 <- if(length(matchI) >1) x3[matchI,] else matrix(x3[matchI,], ncol=2, dimnames=list(rownames(x3)[matchI],colnames(x3)))
    ## need better check for log-data !!!
    if(diff(range(x[,2])) >10 | intensDataIsLin) x[,2] <- log2(x[,2])
  } else {
    if(!silent) message(fxNa,"Can't filter to match experimental input and identified hits")
  }
  modTy <- as.character(lst[[listNa[1]]][,useColN[4]])    #"modif"])
  ## if multiple proteins, need to fuse with prot ID
  if(length(unique(lst[[listNa[1]]][,useColN[3]])) >1) modTy <- paste(as.character(lst[[listNa[1]]][,useColN[3]]),modTy,sep="_")
  uniModTy <- unique(modTy)
  if(length(uniModTy) >20) message(fxNa,"Types of modifications (",length(uniModTy),") : too many for efficient display")
  useColor <- if(length(col) <1) {
    if(length(uniModTy) <9)  1 +(1:length(uniModTy)) else grDevices::rainbow(1.1*(length(uniModTy)))[-1]
    } else col
  ## main plotting
  if(is.data.frame(x)) x <- matrix(as.numeric(as.matrix(x[,1:2])), ncol=2, dimnames=dimnames(x))
  if(is.data.frame(x2)) x2 <- as.matrix(as.numeric(as.matrix(x2[,1:2])), ncol=2, dimnames=dimnames(x2))
  if(is.data.frame(x3)) x3 <- as.matrix(x3)
  if(linPlot) { x[,2] <- 2^x[,2]   ## retransform to intensity linear
    chNu <- is.finite(x[,2])
    if(any(!chNu)) x[which(!chNu),2] <- 0
    x2[,2] <- 2^x2[,2]
    chNu <- is.finite(x2[,2])
    if(any(!chNu)) x2[which(!chNu),2] <- 0
    }
   if(debug) {message(".. xx2"); xx2 <- list(x=x,x2=x2,x3=x3,lst=lst,listNa=listNa,useColN=useColN,useColor=useColor,modTy=modTy,uniModTy=uniModTy)}
  graphics::plot.new()
  if(is.character(x)) message(fxNa," -> x is character ")
  if(is.character(x2)) message(fxNa," -> x2 is character ")
  ploLim <- cbind(xlim=c(1, max(x2[,1], na.rm=TRUE)), ylim=c(1, max(x2[,2], na.rm=TRUE)))
  if(max(x[,1]) > ploLim[2,1]) ploLim[2,1] <- min(ploLim[2,1]*1.15, 1 +max(x[,1]))
  if(max(x[,2]) > ploLim[2,2]) ploLim[2,2] <- min(ploLim[2,2]*1.15, 1 +max(x[,2]))
  if(length(xLim) >1) ploLim[,1] <- xLim 
  graphics::plot.window(xlim=ploLim[,1], ylim=ploLim[,2])
  if(is.null(tit)) tit <- paste(wrMisc::pasteC(names(lst$pep))," : ",nrow(x2)," peaks identified")
  graphics::title(tit)
  if(nrow(lst[[listNa[2]]]) >1) {   # make grey bars for sub-runs
    yCoo <- seq(max(x2[,2])*1.02, max(x2[,2])*0.9, length.out=nrow(lst[[listNa[2]]])+1)
    graphics::mtext(paste(nrow(lst[[listNa[2]]])," identification-subsets shown as grey boxes + no of identifications"), side=3, line=0.6, cex=0.65*cex)      #  
    for(j in 1:nrow(lst[[listNa[2]]])) {
      graphics::rect(lst[[listNa[2]]][j,useColN[9]],yCoo[j+1],lst[[listNa[2]]][j,useColN[10]],yCoo[j], col=grDevices::grey(0.98-(j/40)), border=grDevices::grey(0.85-(j/40), alpha=0.2) )}
    graphics::text(round(rowMeans(lst[[listNa[2]]][,c(useColN[9],useColN[10])])), y=yCoo[-nrow(lst[[listNa[2]]])-1]+diff(yCoo)/2, labels=lst[[listNa[2]]][,useColN[8]], cex=0.8)}
  graphics::lines(x=range(x[,1], na.rm=TRUE), y=c(0,0), col=1)
  graphics::segments(x[,1], 0, x[,1] ,x[,2], lwd=lwd, col=backgrCol)
  legTx <- if(length(replNames) >0) sub(replNames[1], replNames[2], sort(unique(modTy))) else sort(unique(modTy))
  graphics::legend("topright", legend=paste("",legTx), text.col=1, col=useColor, cex=0.8, xjust=0.5, yjust=0.5, lty=1, seg.len=0.6, lwd=4)
  useColor <- useColor[as.numeric(as.factor(modTy))]
  graphics::segments(x2[,1],0,x2[,1],x2[,2], lwd=lwd, col=useColor)
  if(TRUE) graphics::points(x2[,1], x2[,2], cex=0.8)                        # enhace readability, color dot according to score,according to input ?
  graphics::axis(1, line=0.05)         # x-axis 
  graphics::axis(2, line=-2,las=1)     # y-axis
  if(length(yLab) <1) yLab <- colnames(x)[2]
  if(nchar(yLab) <1) yLab <- "Intensity"
  graphics::mtext(yLab, side=2, line=2.2, cex=1.1*cex)      #
  graphics::mtext("(deconvoluted)  m/z", side=1, line=2.2, cex=1.1*cex)   #"mono-isotopic mass"
}
    
