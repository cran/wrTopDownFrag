#' Plot the number of theoretical random fragments
#'
#' This simple function allows plotting the expected number of theoretical fragments from random fragmentation of peptides/proteins (in mass spectrometry).
#' Here, only the pure fragmentation without any variable fragmentation is considered, all fragment-sizes are included (ie, no gating). 
#' For simplicity, possible (variable) modifications like loss of neutrals, etc, are not considered.  
#' 
#' @param x (integer) length (in amino-acids) of input peptides/proteins to be considered  
#' @param tit (character) custom title 
#' @param xlab (character) custom x-axis label 
#' @param ylab (character) custom y-axis label
#' @param col (character or integer) cutsom colors
#' @param log (character) define which axis should be log (use "xy" for drawing both x- and y-axis as log-scale)
#' @param mark (matrix) first column for text and second column for where it should be stated along the top border of the figure (x-coordinate) 
#' @param cexMark (numeric) cex expansion-factor for text from argument \code{mark}
#' @return figure only
#' @seealso \code{\link{AAfragSettings}}
#' @examples
#' marks <- data.frame(name=c("Ubiquitin\n76aa", "Glutamate dehydrogenase 1\n501aa"),
#'   length=c(76,501))
#' plotNTheor(x=20:750, log="", mark=marks)
#' @export
## here simple function to plot the number of theoretical fragments, assume just b- & y- fragments
plotNTheor <- function(x,tit="Number of term and intern fragm",xlab="Number of aa",ylab="",col=2:3,log="",mark=NULL,cexMark=0.75) {
  ## plot number of theoretical fragments
  nTerm <- function(x) (x-1)*2                                   # terminal fragments  
  nInte <- function(x)  sapply(x, function(y) sum((y-2):1))      # internal fragments
  graphics::plot(x,nTerm(x) + nInte(x),type="l", las=1, xlab=xlab, ylab="", log=log, main=tit, col=col[2])
  graphics::mtext("Number of Fragments", side=2, line=5)
  graphics::lines(x, nTerm(x), lty=2, col=col[1])
  if(length(mark) >0) {graphics::abline(v=as.numeric(mark[,2]), lty=2, col=grDevices::grey(0.8))
    graphics::mtext(mark[,1], at=as.numeric(mark[,2]), side=3, line=-1.6 -(1:nrow(mark)), cex=cexMark)}
}
  
