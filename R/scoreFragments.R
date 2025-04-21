#' Scoring For Single Protein : Individual Components
#'
#' Make scoring for single protein : individual components :sameSite,contiguous,prefFragSite,logPeakHeight + combined (sum of scales 0->1)
#'
#' @param resTab (matrix or data.frame) matrix or data.frame of results for SINGLE protein (will use columns 'beg','end','orig','obsInd')
#' @param fragmInp (matrix) experimental m/z values including suppl col(s) to be considered for score, its intensity column/value will be used for 'logInt' in output
#' @param suplTakeLog (logical) if suppl info should be used as log2: if T all supplemental data columns ('fragmInp') will be taken as log2
#' @param j (integer) which column of fragmInp has m/z values, the following column is assumed as peak-intensity
#' @param useResCol (character) column-names from resTab to be used
#' @param prefFragPat (matrix) for preferential fragmentation rules (see \code{.prefFragPattern()})
#' @param contigTermFragWe (numeric, length=1) weight to add for terminal fragments at 'sc.complemFra' (since they cannot match other fragments beyond the protein limits)
#' @param figDraw (logical) make additional figure 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages and objects exportet to current session for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with matrix $scaled (combined and indiv rescaled scores) and $raw (matching lines of 'resTab'; 'index' refers to predictedIndex) 
#' @seealso \code{\link{identifyPepFragments}}, \code{\link{scoreProteinFragments}}
#' @examples
#' tab2 <- matrix(c("20","2","13","11","3","10","4", "PT","PE","EP","DE","PEP","IDE","PEPT", 
#'   rep(c("PEPTIDE","protP"),each=7), c("inter","Nter","Cter")[c(1,2,1,3,2,3,2)], 
#'   c(3,1,2,6,1,5,1, 4,2,3,7,3,7,4), "E",NA,"P","I",NA,"T",NA, "I","P","T",NA,"T",NA,"I", 
#'   c(1,6,6,20,7,19,8), c(-0.094312,-0.14707,-0.14707,0.08641,0.0084762,-0.10965,0.057087), 
#'   rep(2,7)), nrow=7, dimnames=list(NULL,c("predInd","seq","orig","origNa","ty","beg","end",
#'   "precAA","tailAA","obsInd","ppmToPred","mass")))
#' tab2 <- cbind(tab2, seqNa=paste0(tab2[,"origNa"],".",tab2[,"beg"],"-",tab2[,"end"]),Abundance=1)
#' rownames(tab2) <- paste0(tab2[,"origNa"],".", tab2[,"beg"],"", tab2[,"end"])
#' obsMassX <- cbind(a=c(199.1077,296.1605,397.2082,510.2922,625.3192),
#'   b=c(227.1026,324.1554,425.2031,538.2871,653.3141),
#'   x=c(729.2937,600.2511,503.1984,402.1507,289.0666),
#'   y=c(703.3145,574.2719,477.2191,376.1714,263.0874))
#' 
#' (outF <- scoreFragments(tab2, fragmInp=cbind(as.numeric(obsMassX), Abundance=1)))
#' 
#' @export
scoreFragments <- function(resTab, fragmInp, suplTakeLog=TRUE, j=NULL, useResCol=c("orig","seq","precAA","tailAA","beg","end","ppmToPred","obsInd","predInd","Abundance"),
  prefFragPat=NULL, contigTermFragWe=0.5, figDraw=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## make scoring : individual components :sameSite,contiguous,prefFragSite,logPeakHeight + combined (sum of scales 0->1)
  ## return list with matrix $scaled (combined and indiv rescaled scores) and $raw (matching lines of 'resTab')
  ## 'resTab' .. matrix or data.frame of results for SINGLE protein (will use columns 'beg','end','orig','obsInd')
  ## 'fragmInp' .. initial m/z input including suppl col(s) to be considered for score
  ## 'suplTakeLog'.. if suppl info should be used as log2: if T all supplemental data columns ('fragmInp') will be taken as log2
  ## 'j' .. which column of fragmInp has m/z values, the following column is assumed as peak-intensity
  ## 'prefFragPat' .. custom preferential fragmentation pattern, if NULL use default from .prefFragPattern()
  ## 
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="scoreFragments")
  minNoSiteAcc <- 2       # mimimum accumulation of sites; could be chosen more dynamic depending on size of protein(s) and no of hits found
  joinStat <- contrPerFrag <- prefSco <- NULL           # initialize
  chCol <- match(useResCol, colnames(resTab))
  if(!silent && any(is.na(chCol))) message(fxNa,"Can't find column(s) ",wrMisc::pasteC( useResCol[which(is.na(chCol))],quoteC="'")," in 'resTab'")
  nAA <- nchar(as.character(resTab[,useResCol[2]]))                              #"orig"
  
   if(debug) {message(fxNa, "xxscoreFragments0a"); xxscoreFragments0a <- list(resTab=resTab,fragmInp=fragmInp,useResCol=useResCol,nAA=nAA,minNoSiteAcc=minNoSiteAcc,prefFragPat=prefFragPat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA

  datBegEnd <- matrix(as.integer(as.matrix(resTab[,useResCol[6:7]])), ncol=2, dimnames=list(rownames(resTab),useResCol[6:7]))

   if(debug) {message(fxNa, "xxscoreFragments0b"); xxscoreFragments0b <- list(resTab=resTab,datBegEnd=datBegEnd,fragmInp=fragmInp,useResCol=useResCol,nAA=nAA,minNoSiteAcc=minNoSiteAcc,prefFragPat=prefFragPat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA
  ## COMPLEM.FRA
  ## reassemble puzzle based on fragmentation ==> sc.complemFra
  ## complementary fragments (a+b=c)  ==> sc.complemFra  
  loc <- as.matrix(resTab)
  loc2 <- matrix(as.integer(loc[,c("beg","end")]), ncol=2, dimnames=list(loc[,"seqNa"],c("beg","end")))
  complemFra <- log2(1 +countChildrenParent(loc2, silent=silent, callFrom=fxNa))/2            # extremes may go higher than 1.0; here median=0, max=1.16
  names(complemFra) <- rownames(resTab)
  
  ## PREF.FRAG
  ## preferred breakage sites  ==> sc.prefFrag
  prefSco <- scorePrefFrag(resTab, useCol=useResCol[c(2:4)], prefFragPat=prefFragPat, silent=silent, debug=debug)
  if(debug) {message(fxNa, "xxscoreFragments1"); xxscoreFragments1 <- list(resTab=resTab,datBegEnd=datBegEnd,complemFra=complemFra,contrPerFrag=contrPerFrag,fragmInp=fragmInp,useResCol=useResCol,nAA=nAA,minNoSiteAcc=minNoSiteAcc,prefFragPat=prefFragPat,joinStat=joinStat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA
  ## SAME.SITE
  ## same fragmentation sites (should not count true start and end of protein ?? - so far counted)  ==> sc.sameSite
  fraRes <- wrMisc::countSameStartEnd(datBegEnd,minFreq=minNoSiteAcc) 
  fraRes[which(is.na(fraRes))] <- 0
  fragmInp <- as.matrix(fragmInp)
  if(ncol(fragmInp) <2) {fragmInp <- cbind(fragmInp, obsAbund=rep(1,nrow(fragmInp)))
    if(!silent) message(fxNa," no abundance values provided, set all to 1.0")}
  if(debug) {message(fxNa, "xxscoreFragments1b"); xxscoreFragments1b <- list(resTab=resTab,datBegEnd=datBegEnd,complemFra=complemFra,contrPerFrag=contrPerFrag,fragmInp=fragmInp,useResCol=useResCol,nAA=nAA,minNoSiteAcc=minNoSiteAcc,prefFragPat=prefFragPat,joinStat=joinStat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA

  identifRa <- range(as.numeric(as.character(resTab[,"mass"])), na.rm=TRUE)
  fraInp2 <- fragmInp[which(fragmInp[,1] >= identifRa[1] & fragmInp[,1] <= identifRa[2]),]   # filter input to range of real identifications for normalizing sameSite
  fraRe2 <- fraRes[,c("beg.n","end.n")] -2
  ch2 <- fraRe2 ==0
  if(any(ch2)) fraRe2[which(ch2)] <- 0.3
  ch2 <- fraRe2 <= 0
  if(any(ch2)) fraRe2[which(ch2)] <- 0

  if(debug) {message(fxNa, "xxscoreFragments1c\n"); xxscoreFragments1c <- list(resTab=resTab,datBegEnd=datBegEnd,complemFra=complemFra,contrPerFrag=contrPerFrag,fragmInp=fragmInp,useResCol=useResCol,nAA=nAA,minNoSiteAcc=minNoSiteAcc,prefFragPat=prefFragPat,joinStat=joinStat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA

  sameSiteOld <- 25*rowSums(fraRes[,c("beg.n","end.n")]^1.1)/nrow(fraInp2)       # extremes may go higher than 1.0 (here 1.26), almost usable
  sameSite <- 50*rowSums(log(1 +fraRe2^1.4))/nrow(resTab)
  sameSiteLin <- 20*rowSums(fraRe2^1.1)/nrow(resTab)                           # instance of 3 counts at both sides score at 0.06 ie too low
  if(debug) {message(fxNa, "xxscoreFragments1d"); xxscoreFragments1d <- list(resTab=resTab,datBegEnd=datBegEnd,complemFra=complemFra,sameSite=sameSite,fraRes=fraRes,fraRe2=fraRe2,fraInp2=fraInp2  ,contrPerFrag=contrPerFrag,fragmInp=fragmInp,useResCol=useResCol,nAA=nAA,minNoSiteAcc=minNoSiteAcc,prefFragPat=prefFragPat,joinStat=joinStat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA
  
  fraRes <- cbind(fraRes, accSite=sameSite)           
  fraRes <- fraRes[match(rownames(resTab), rownames(fraRes)),]                    # need to re-adjust order !!
  #not used#fraRes <- cbind(fraRes, complemFra/2)                                # contrPerFrag already in order of resTab
  ## CHARGE.CATCH
  ## charge Catch  ==> sc.chargeCatch
  chargeCatch <- scoreChargeCatch(resTab, scale01=FALSE, silent=silent, callFrom=callFrom)
  chCh <- chargeCatch >0
  if(any(chCh)) chargeCatch[which(chCh)] <- 1                                    # reduce counting to simple 0 or 1
  ## COMPLEMENTARY FRAGMENTS
  ## complementary fragments (a+b=c)  ==> sc.complemFra
  comSco <- cbind(sameSite=fraRes[,"accSite"], prefFrag=prefSco, chargeCatch=chargeCatch, complemFra=complemFra)             #/nrow(resTab)
  comSco <- cbind(index=as.integer(as.character(resTab[,useResCol[9]])), comSco, ppm=as.numeric(as.character(resTab[,useResCol[7]]))) 
  ## try adding intensity data  ==> logInt
  resTab <- as.matrix(resTab)        # so far data.frame
  if(!is.null(fragmInp)) {           ## find column with experim intensities and rename 'obsInt' for adding to results
    ## check formatting of 'fragmInp' and rename columns if necessary
    if(length(dim(fragmInp)) <1) fragmInp <- as.matrix(fragmInp)   # can't hold intensity info if only 1 column
    if(!useResCol[10] %in% colnames(fragmInp)) {                   # check for absence of column with experim intensities in experim values
      ## no column names 'Abundance', start looking for column with intensity data
      if(is.null(colnames(fragmInp))) colnames(fragmInp)[2] <- "obsAbund" else {
        chAbSyn <- c("Abundance","Abund","Intensity","Intens","Int","Peak.Height","PeakHeight","Height") %in% colnames(fragmInp)
        if(any(chAbSyn)) colnames(fragmInp)[which(chAbSyn)[1]] <- "obsAbund" else {
          ch2 <- lapply(c("Abundance","Abund","Intensity","Intens","Int","Peak.Height","PeakHeight","Height"),grep,colnames(fragmInp))
          chAny <- sapply(ch2, length) >0
          if(any(chAny)) {
            if(!silent) message(fxNa,"Using column '",colnames(fragmInp)[ch2[[which(chAny)[1]]]],"' as observed intensities")
            colnames(fragmInp)[ch2[[which(chAny)[1]]]] <- "obsAbund"} else {    
          chNa <- colnames(fragmInp) %in% "" 
          if(any(chNa)) { colnames(fragmInp)[which(chNa)[1]] <- "obsAbund"
            if(!silent) message(fxNa,"Abundance data most likely missing, using column without name as abundance")
          } else {
            if(!silent) message(" +++\n",fxNa,"TROUBLE AHEAD ... Can't find in argument 'fragmInp' column named '",useResCol[8],"' !\n +++")
          }
    } }}
          
    } else colnames(fragmInp)[which(colnames(fragmInp)==useResCol[10])] <- "obsAbund" 
    if(debug) {message(fxNa, "xxscoreFragments1e"); xxscoreFragments1e <- list(resTab=resTab,comSco=comSco,datBegEnd=datBegEnd,sameSite=sameSite,complemFra=complemFra, fraRes=fraRes,fragmInp=fragmInp,useResCol=useResCol,prefFragPat=prefFragPat,joinStat=joinStat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA
    intens <- wrMisc::convToNum(fragmInp[as.integer(resTab[,"obsInd"]),"obsAbund"])
    comSco <- cbind(comSco, logintens=if(suplTakeLog) log2(intens) else intens)
    }
  ## check for entire columns with all NA or all 0 (and replace by 0.01 if so to avoid NaN or Inf), eg if no charge-catching AAs found 
  ch0 <- colSums(comSco==0) ==nrow(comSco) | colSums(is.na(comSco))==nrow(comSco)  
  if(any(ch0)) comSco[,which(ch0)] <- 0.01
  ## COMBINE (old: NORMALIZE scores to similar scales)
  maxSc <- apply(comSco[,c("prefFrag","chargeCatch","complemFra","sameSite")], 2, max, na.rm=TRUE)  
  comSc2 <- cbind(prefFrag=comSco[,"prefFrag"], chargeCatch=comSco[,"chargeCatch"], complemFra=comSco[,"complemFra"],sameSite=comSco[,"sameSite"])  	  
  colnames(comSc2) <- paste0("sc.",colnames(comSc2))
  chVarVal <- apply(comSc2, 2, function(x) length(unique(x)) >1)
  if(all(!chVarVal)) {message(fxNa,"All score values have no variation => NOT USEFUL"); sco4 <- rep(0,nrow(comSc2))
  } else sco4 <- if(sum(chVarVal)==1)comSc2[,which(chVarVal)] else rowSums(comSc2[,which(chVarVal)])/sum(chVarVal)               ## the combined score, not normalized to max 
  if(debug) {message(fxNa, "xxscoreFragments2"); xxscoreFragments2 <- list(sco4=sco4,maxSc=maxSc,comSc2=comSc2,comSco=comSco,resTab=resTab,fragmInp=fragmInp,fraRes=fraRes,useResCol=useResCol,prefFragPat=prefFragPat,joinStat=joinStat,prefSco=prefSco,suplTakeLog=suplTakeLog)}     #,nAA=nAA
  comSc2 <- cbind(sco4=sco4, comSc2, index=as.integer(as.character(resTab[,useResCol[9]])))
  # message(fxNa, " done main scoring \n")
  if(figDraw) { #pairs(comSco)
    col1 <- grDevices::rgb(red=c(141,72,90,171, 220,253,244,255), green=c(129,153,194,221, 216,174,109,0), blue=c(194,203,185,164, 83,97,67,0), maxColorValue=255) # pale purple, pale green, red
    his <- graphics::hist(comSc2[,1], breaks=seq(0,max(comSc2[,1],na.rm=TRUE) +0.1,by=0.05), main="combined score4 ", col=col1[c(rep(1:7,each=2),rep(8,6))], xlab="score")
    graphics::mtext(paste(nrow(resTab)," fragments, using prefFrag,chargeCatch,complemFrags"), cex=0.8, side=3, line=0.2)
    graphics::legend("topright", legend=paste(">",his$breaks[(1:8)*2-1]), text.col=1, col=col1, cex=0.8, xjust=0.5, yjust=0.5, lty=1,seg.len=0.6,lwd=4)
  }  
  list(orig=comSco,scaled=comSc2)} 
      
