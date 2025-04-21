#' Scoring Of Identifications (For Multi-Protein Queries)
#'
#' Make scoring for multiple protein queries: individual components :sameSite,contiguous,prefFragSite,logPeakHeight + combined (sum of scales 0->1)
#'
#' @param resTab (matrix or data.frame) identification results  (will use columns 'beg','end','orig','obsInd')
#' @param fragmInp (numeric vector or matrix) experimental m/z values, may include suppl col(s) to be considered for score (se argument 'j')
#' @param j (integer) which column of fragmInp has m/z values, the following column is assumed as peak-intensity
#' @param useCol (character) colnames from resTab to be used, 1st posiion must be present as column of 'resTab' and must represent name of original, used for splitting by proteins (eg original input protein sequence), will be passed to \code{scoreFragments}  
#' @param prefFragPat (matrix) for preferential fragmentation rules (see \code{.prefFragPattern()})
#' @param contigTermFragWe (numeric, length=1) weight to add for terminal fragments (since they cannot match other fragments beyond the protein limits)
#' @param returnCombined (logical)
#' @param figDraw (logical) make additional figure 
#' @param silent (logical) suppress messages 
#' @param callFrom (character) allow easier tracking of message(s) produced 
#' @param debug (logical) for bug-tracking: more/enhanced messages and intermediate objects written in global name-space 
#' @return This function returns a list with matrix $scaled (combined and indiv rescaled scores) and $raw (matching lines of 'resTab')
#' @seealso \code{\link{scoreFragments}}, \code{\link{identifyPepFragments}}
#' @examples
#' tab2 <- matrix(c("20","2","13","11","3","10","4", "PT","PE","EP","DE","PEP","IDE","PEPT", 
#'   rep(c("PEPTIDE","protP"),each=7), c("inter","Nter","Cter")[c(1,2,1,3,2,3,2)], 
#'   c(3,1,2,6,1,5,1, 4,2,3,7,3,7,4), "E",NA,"P","I",NA,"T",NA, "I","P","T",NA,"T",NA,"I", 
#'   c(1,6,6,20,7,19,8), c(-0.094312,-0.14707,-0.14707,0.08641,0.0084762,-0.10965,0.057087), 
#'   rep(2,7)), nrow=7, dimnames=list(NULL, c("predInd","seq","orig","origNa","ty","beg","end",
#'   "precAA","tailAA","obsInd","ppmToPred","mass")))
#' tab2 <- cbind(tab2, seqNa=paste0(tab2[,"origNa"],".",tab2[,"beg"],"-",tab2[,"end"]),Abundance=1)
#' rownames(tab2) <- paste0(tab2[,"origNa"],".", tab2[,"beg"],"", tab2[,"end"])
#' obsMassX <- cbind(a=c(199.1077,296.1605,397.2082,510.2922,625.3192),
#'   b=c(227.1026,324.1554,425.2031,538.2871,653.3141),
#'   x=c(729.2937,600.2511,503.1984,402.1507,289.0666),
#'   y=c(703.3145,574.2719,477.2191,376.1714,263.0874))
#' 
#' (outF <- scoreFragments(tab2, fragmInp=cbind(as.numeric(obsMassX), Abundance=1)))
#' (out <- scoreProteinFragments(tab2, fragmInp=cbind(as.numeric(obsMassX), Abundance=1)))
#' 
#' @export
scoreProteinFragments <- function(resTab, fragmInp=NULL, j=2, useCol=c("orig","precAA","tailAA","beg","end","ppmToPred","obsInd","predInd"),
  prefFragPat=NULL, contigTermFragWe=0.5, returnCombined=TRUE, figDraw=TRUE, silent=FALSE, callFrom=NULL,debug=FALSE) {
  ## make scoring (for multiple protein queries): individual components :sameSite,contiguous,prefFragSite,logPeakHeight + combined (sum of scales 0->1)
  ## return list with matrix $scaled (combined and indiv rescaled scores) and $raw (matching lines of 'resTab')
  ## 'resTab' .. matrix or data.frame of results  (will use columns 'beg','end','orig','obsInd')
  ## 'fragmInp' .. initial m/z input including suppl col(s) to be considered for score
  ## 'j' .. which column of fragmInp has m/z values, the following column is assumed as peak-intensity
  ## 'useCol' .. colnames from resTab to be used
  ## 'prefFragPat' .. custom preferential fragmentation pattern, if NULL use default from .prefFragPattern()
  ## 'figDraw' .. hist of combined score
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="scoreProteinFragments")
  chCol <- match(useCol[1:min(7,length(useCol))],colnames(resTab))
  if(any(is.na(chCol))) {if(all(is.na(chCol))) stop(fxNa,"Can't find ANY columns of 'resTab'") else {
    message(fxNa,"Trouble ahead ?  Can't find columns ",wrMisc::pasteC(useCol[which(is.na(chCol))],quoteC="'"))
    }} else if(!silent) message(fxNa,"All column-names of 'resTab' found.. OK")
  if(debug) {message(fxNa, "xxscoreProteinFragments1"); xxscoreProteinFragments1 <- list(resTab=resTab,fragmInp=fragmInp,useCol=useCol,prefFragPat=prefFragPat,returnCombined=returnCombined,j=j,chCol=chCol)}
  chNOri <- unique(resTab[,useCol[1]])                         # "orig"
  resTab <- cbind(resTab, index=1:nrow(resTab))
  ## split results by input protein
  resByPr <- by(resTab, resTab[,useCol[1]], as.data.frame)     # "orig"  
  ## BASE/MAIN SCORING
  ## sometimes problem  scoreFragments() not attaching intensity !
  scores <- lapply(resByPr, scoreFragments, fragmInp=fragmInp, suplTakeLog=TRUE, j=j, useResCol=c(useCol[1],"seq",useCol[-1],"Abundance"),
    prefFragPat=prefFragPat, contigTermFragWe=contigTermFragWe, figDraw=FALSE, silent=silent,debug=debug,callFrom=fxNa)
  if(length(dim(fragmInp)) <1) fragmInp <- as.matrix(fragmInp)
  if(debug) {message(fxNa, "xxscoreProteinFragments2"); xxscoreProteinFragments2 <- list(scores=scores,resByPr=resByPr,resTab=resTab,fragmInp=fragmInp,useCol=useCol,prefFragPat=prefFragPat,returnCombined=returnCombined)}  
  #sapply(scores[[1]],length)
  if(length(scores) >1) {
    ## multiple proteins, need to attach/combine lists
    out <- matrix(nrow=sum(sapply(scores,function(x) nrow(x$scaled))),ncol=ncol(scores[[1]]$scaled)+1)
    sta <- 0
    for(i in 1:length(scores)) {
      out[sta+(1:nrow(scores[[i]]$scaled)),] <- cbind(scores[[i]]$scaled,logInt= if(length(fragmInp) >0) scores[[i]]$orig[,"logintens"])
      sta <- nrow(scores[[i]]$scaled) }
  } else out <- cbind(scores[[1]]$scaled, logInt=if(length(fragmInp) >0) scores[[1]]$orig[,"logintens"])
  rownames(out) <- unlist(sapply(scores, function(x) rownames(x$orig)))
  if(is.null(colnames(out))) colnames(out) <- c(colnames(scores[[1]]$scaled), if(ncol(fragmInp) >1) "logInt")   #[(1:ncol(fragmInp))+ncol(scores[[1]]$scaled)]
  colnames(out) <- sub("\\.$","",colnames(out))
  if(debug) {message(fxNa, "xxscoreProteinFragments3"); xxscoreProteinFragments3 <- list(out=out,scores=scores,resTab=resTab,fragmInp=fragmInp,useCol=useCol,prefFragPat=prefFragPat,returnCombined=returnCombined)}
  if(figDraw) { 
    col1 <- grDevices::rgb(red=c(141,110,72,90,171, 220,253,244,255), green=c(129,133,153,194,221, 216,174,109,0), blue=c(194,198,203,185,164, 83,97,67,0), maxColorValue=255) # 9 steps : pale purple,blue, pale green, red
    his <- graphics::hist(out[,"sco4"],breaks=seq(0,max(out[,"sco4"], na.rm=TRUE) +0.1, by=0.05), main="combined score4 ", xlab="score", col=col1[c(rep(1:8,each=2),rep(9,4))])     #[1:length(his$counts)]
    graphics::mtext("using sameSite,contiguous,prefFragSite,logPeakHeight",cex=0.8,side=3,line=0.2)
    graphics::legend("topright", legend=paste(">",his$breaks[(1:9)*2-1]),text.col=1,col=col1,cex=0.8,xjust=0.5,yjust=0.5,lty=1,seg.len=0.6,lwd=4)
    }  #
  resTab <- .chColNa(resTab, colNa="index", callFrom=fxNa)
  out <- .chColNa(out, colNa="index", rmCol=colnames(resTab), callFrom=fxNa)
  if(debug) {message(fxNa, "xxscoreProteinFragments4"); xxscoreProteinFragments4 <- list(out=out,scores=scores,resTab=resTab,fragmInp=fragmInp,useCol=useCol,prefFragPat=prefFragPat,returnCombined=returnCombined)}
  if(returnCombined) { 
    if(is.null(rownames(resTab))) {
      rownames(resTab) <- paste0(resTab[,useCol[1]],".",resTab[,useCol[4]],"-",resTab[,5])
      if(!silent) message(fxNa,"Rownames for identified peptides missing, trying to reconstruct") }
    resTab <- data.frame(resTab, rowNa=rownames(resTab), fraNa=sub("-",".",rownames(resTab)), stringsAsFactors=FALSE)
    out <- data.frame(out, fraNa=rownames(out), stringsAsFactors=FALSE)
    colnames(out)[ncol(out)] <- "fraNa"
    ## make sure the correct version of the fragment names is used for merge()
    if(!"fraNa" %in% colnames(resTab)) {
       chNa <- wrMisc::naOmit(match(c("rowNa","index"), colnames(resTab)))
       if(length(chNa) >0) colnames(resTab)[chNa[1]] <- "fraNa" else stop("cannot find colnames to match scores and identifed peptides")
    }
    
    out <- merge(resTab, out, by="fraNa")
    tm2 <- c(grep("^rowNa",colnames(out))[-1],grep("^seqNa",colnames(out)), grep("^index",colnames(out)))
    out <- out[,-1*tm2]
     if(debug) {message(fxNa, "xxscoreProteinFragments5"); xxscoreProteinFragments5 <- list(out=out,scores=scores,resTab=resTab,fragmInp=fragmInp,useCol=useCol,prefFragPat=prefFragPat,returnCombined=returnCombined)}
    ## want 'fraNa' as 1st colum
    rowNCol <- match(c("fraNa","rowNa"), colnames(out))
    out <- if(any(is.na(rowNCol))) out[,c(wrMisc::naOmit(rowNCol),(1:ncol(out))[-1*wrMisc::naOmit(rowNCol)])] else {
      out[,c(rowNCol[1],(1:ncol(out))[-1*rowNCol])]}
    colnames(out)[1] <- "fraNa"   
    ## real example : prefFrag =3 levels; chargeCatch =11 levels; complemFra=81 levels
    ## on PCA prefFrag is quite orthogonal to other factors, complemFrag is redundant (& smaller) to chargeCatch
    } else out <- out[order(as.integer(out[,"index"])),]
  out }


#' Check Column Names from Matrix Or data.frame
#'
#' Check matrix or data.frame for containing columns with specified name and optionally remove other columns (by their names)
#'
#' @param x (matrix or data.frame) 
#' @param colNa (character) 
#' @param rmCol (character)
#' @param callFrom (character) allow easier tracking of message(s) produced 
#' @return This function returns a matrix or data.frame with adjusted columns
#' @seealso \code{\link{scoreFragments}}
#' @examples
#' ma1 <- matrix(1:6, nrow=2, dimnames=list(NULL, c("index","aa","bb")))
#' .chColNa(ma1)
#' .chColNa(ma1, colNa="zz", rmCol="aa") 
#' @export
.chColNa <- function(x, colNa="index", rmCol=NULL, callFrom=NULL) {     # check if 'x' contains column named colNa and correct to single instance
  xNa <- deparse(substitute(x))
  colNa <- colNa[1]                                         # must be of length=1
  chInd <- colnames(x)==colNa
  if(all(!chInd)) message(callFrom,"Trouble ahead, no column named '",colNa,"' found in '",xNa,"'")
  if(sum(chInd) >1) x <- x[,-1*(which(chInd)[-1])]
  if(length(rmCol) >0) {if(colNa %in% rmCol) rmCol <- rmCol[-which(rmCol==colNa)]}
  if(length(rmCol) >0) {    
    chCol <- colnames(x) %in% rmCol    
    if(any(chCol)) { dimNa <- dimnames(x)
      x <- x[,-which(chCol)]
      if(length(dim(x)) <2) x <- matrix(x,nrow=length(dimNa[[1]]), dimnames=list(dimNa[[1]], dimNa[[2]][-which(chCol)]))}
  }
  x }
  
