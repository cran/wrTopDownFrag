#' plot preferential fragmenation pattern
#'   
#' Plot preferential fragmenation pattern equivalent to Fig 1b of Haverland et al 2017 (J Am Soc Mass Spectrom)  
#'  
#' @param prefPat (matix) 
#' @param namesCex (numeric) expansion factor cex for display of AA-names
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages and objects exportet to current session for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a figure   
#' @seealso  \code{\link{scoreFragments}}  
#' @examples
#' plotPrefFragPat(.prefFragPattern())
#' @export
plotPrefFragPat <- function(prefPat, namesCex=0.8, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## plot preferential fragmenation pattern equivalent to Fig 1b of Haverland et al 2017 (J Am Soc Mass Spectrom)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotPrefFragPat")
  AAs <- sort(names(wrProteo::AAmass()[1:20]))
  dat <- matrix(0, nrow=length(AAs), ncol=length(AAs), dimnames=list(paste(AAs,"x",sep="."),paste("x",AAs,sep=".")))
  ind <- cbind(Ct=match(prefPat[,1],AAs), Nt=match(prefPat[,2],AAs),prefPat[,3])
  for(i in unique(ind[,1])) {j <- which(ind[,1]==i); dat[i,ind[j,2]] <- ind[j,3]}
  colGra <- c(RColorBrewer::brewer.pal(9,"Blues")[7:3], grDevices::gray(0.87), RColorBrewer::brewer.pal(7,"Reds")[2:6])
  legP <- seq(0,1,length.out=ncol(dat))
  legL <- legP[-1] -diff(legP[1:2])/2
  tit <- "Preferential Fragmentation Deconstructed by Residue Pair"
  graphics::image(t(dat[nrow(dat):1,]),col=colGra,xaxt="n",yaxt="n",xlab="X|x'  (cut before X)",ylab="x|X'  (cut after x)",main=tit)
  graphics::abline(h=legL,col=grDevices::grey(0.7)); graphics::abline(v=legL, col=grDevices::grey(0.7)); 
  graphics::mtext(at=legP+legL[1]/3,AAs, side=1, cex=namesCex, adj=1, las=1, line=0.4)            # x labels
  graphics::mtext(at=legP+legL[1]/5,rev(AAs), side=2, cex=namesCex, adj=0.5,las=2, line=0.7)     # y labels
  graphics::mtext("blue -> grey -> red .. increased prefered fragmentation", side=3, cex=0.8, line=0.1) 
}
   
