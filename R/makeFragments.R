#' Make Terminal And Internal Fragments From Proteins  
#'
#' Makes terminal and internal fragments based on protein-sequence and present as matrix including heading and/or tailing amino-acid or theoretical molecular mass of all fragments.
#' As the number of theoretically possible fragments increases with the size of the peptide/protein treated it is recommended to adopt arguments like \code{masFragSize} to 
#' realistic values for the type of mass spectrometer used, since efficient filtering will reduce considerably the amount of memory (RAM) needed and will improve overal performance.
#' 
#' @param protTab (character or matrix) named vector of protein-seqences to fragment or 
#'   or, matrix with 1st column 'ma' with mass and 2nd column 'se' with protein-/peptide-sequence (optionally sequence(s) also as rownames or 3rd column with 'trivial' names)
#' @param minFragSize (integer) minimum number of amino-acids for being considered
#' @param maxFragSize (integer) maximum number of amino-acids for being considered
#' @param internFra (logical) toggle if internal framents will be produced or not
#' @param knownMods (character) optional custom alternative to \code{AAfragSettings(ou="all")$knownMods}  
#' @param redRedundSeq (logical) reduce redundant sequences to 1st appearance in all further treatments
#' @param prefFragPat (matrix) for preferential fragmentation rules (see also \code{.prefFragPattern})
#' @param remNonConfPrefFragm (logical) allows to remove (peptide-)fragments non conform with preferential fragmentation rules (using \code{evalIsoFragm})
#' @param ambigLab (character) text-labels for ambiguities (first for duplicated sequences second for iso-mass)
#' @param massTy (character) default 'mono' for mono-isotopic masses (alterative 'average')
#' @param specModif (list) supplemental custom fixed or variable modifications (eg Zn++ at given residue) 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return matrix with fragment sequence, mass, start- and end-position, heading and tailing AA (or NA if terminal fragment)
#' @seealso \code{\link{makeFragments}};  \code{\link{evalIsoFragm}}, from package \href{https://CRAN.R-project.org/package=wrProteo}{wrProteo} \code{\link[wrProteo]{convAASeq2mass}}, \code{\link[wrProteo]{AAmass}}, \code{\link[wrProteo]{massDeFormula}}
#' @examples
#' protP <- c(protP="PEPTIDE")
#' pepT1 <- makeFragments(protTab=protP, minFragSize=2, maxFragSize=9, internFra=TRUE)
#' tail(pepT1)
#' 
#' protP2 <- cbind(se="PEPTIDE", ma="1304.7088503626")
#' pepT2 <- makeFragments(protTab=protP2, minFragSize=2, maxFragSize=9, internFra=TRUE)
#' tail(pepT1)
#' @export
makeFragments <- function(protTab, minFragSize=6, maxFragSize=300, internFra=TRUE, knownMods=NULL, redRedundSeq=FALSE, prefFragPat=NULL,
  remNonConfPrefFragm=TRUE, ambigLab=c(duplSequence="duplSequence",isoMass="isoMass"), massTy="mono",specModif=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="makeFragments")
  docTi <- rep(Sys.time(),7)                        #
  msg <- "Expecting matrix with 3 columns (name,sequence,mass) and >=1 line"
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(length(dim(protTab)) !=2) { 
    if(is.null(names(protTab))) names(protTab) <- paste0("p", 1:length(protTab))       # for checking 
    proMass <- try(wrProteo::convAASeq2mass(protTab, massTy=massTy))
    protTab <- cbind(se=protTab, ma=if(inherits(proMass, "try-error")) NA else as.character(proMass), na=names(protTab))
  }
  
  if(is.null(rownames(protTab))) rownames(protTab) <- protTab[,2]
  if(is.null(knownMods)) knownMods <- AAfragSettings(outTy="all")$knownMods    #
  names(docTi) <- c("ini","fragmentSeq",".exNamesTyDeList","convAASeq2mass","findRepeated")
  
  docTi[2] <- Sys.time()
  prefNaCol <- if(ncol(protTab) >2) 3  else 1   # column of protTab to use for 'trivial' name
  pep2 <- apply(protTab, 1, function(x) fragmentSeq(x[1], minSize=minFragSize, maxSize=maxFragSize, internFragments=internFra, 
    separTerm=TRUE, keepRedSeqs=TRUE, prefName=x[prefNaCol], callFrom=fxNa, silent=silent))                            # takes ~18% time !
  names(pep2) <- protTab[,1]
  docTi[3] <- Sys.time()
  if(debug) {message(fxNa," .. mFr1"); mFr1 <- list(protTab=protTab,minFragSize=minFragSize,maxFragSize=maxFragSize,internFra=internFra,knownMods=knownMods,massTy=massTy,pep2=pep2)}
  ## note: default fragmentSeq may already remove redudant sequences, ie here set keepRed=T to allow adding prefix !
  ## organize into table of frags for all proteins, indic if full, Nterm ... protOrig & position ... seq
  pepTab <- .exNamesTyDeList(pep2, fullSeq=protTab[,2])                          # cols  c("seq","orig","origNa","ty","seqNa","beg","end","precAA","tailAA","mass")
  docTi[3] <- Sys.time()
  ## determine duplicated sequences, determine mass (but can't separate yet to nonredundant set, otherwise problem with preferential fragmentation sites)
  duplS1 <- duplicated(pepTab[,"seq"], fromLast=FALSE)
  if(any(duplS1)) { 
    duplS2 <- duplicated(pepTab[,"seq"], fromLast=TRUE)
    pepTab[which(duplS1 | duplS2),"ambig"] <- ambigLab[1]    # "duplSequence"    
    pepMa <- wrProteo::convAASeq2mass(pepTab[which(!duplS1),"seq"], massTy=massTy, callFrom=fxNa) - wrProteo::.atomicMasses()["e",massTy]      # also subtract 1 electron mass for making (single charge) ions
    pepTab[which(!duplS1),"mass"] <- pepMa
    ## propagate mass ..
    pepTab[which(duplS1),"mass"] <- pepTab[which(duplS2)[match(pepTab[which(duplS1),"seq"], pepTab[which(duplS2),"seq"])],"mass"]  # reduce search space of match 
  } else {    
    pepTab[,"mass"] <- wrProteo::convAASeq2mass(pepTab[,"seq"], massTy=massTy, callFrom=fxNa) - wrProteo::.atomicMasses()["e",massTy]          # also subtract 1 electron mass for making (single charge) ions  
  }
  if(debug) {message(fxNa," .. mFr2"); mFr2 <- list(protTab=protTab,minFragSize=minFragSize,maxFragSize=maxFragSize,internFra=internFra,knownMods=knownMods,massTy=massTy,pep2=pep2,pepTab=pepTab)}
  ## determine precAA AND tailAA
  heaPo <- as.integer(pepTab[,"beg"])
  chHe <- heaPo >1
  if(any(chHe)) {chHe <- which(chHe); pepTab[chHe,"precAA"] <- substr(pepTab[chHe,"orig"], heaPo[chHe]-1, heaPo[chHe]-1)}
  taiPo <- as.numeric(pepTab[,"end"])
  chTa <- taiPo < nchar(pepTab[,"orig"])
  if(any(chTa)) {chTa <- which(chTa); pepTab[chTa,"tailAA"] <- substr(pepTab[chTa,"orig"], taiPo[chTa]+1, taiPo[chTa]+1)}  
  if(debug) {message(fxNa," .. mFr3"); mFr3 <- list(protTab=protTab,minFragSize=minFragSize,maxFragSize=maxFragSize,internFra=internFra,knownMods=knownMods,massTy=massTy,
    pep2=pep2,pepTab=pepTab,taoPo=taiPo,chTa=chTa,prefFragPat=prefFragPat)}
  ##  find iso-fragments (later: choose preferential cleavage sites xD.xx, xE.xx, xx.Px, need heading&tailing AA)
  pepTab <- pepTab[order(as.numeric(pepTab[,"mass"])),]
  chMa <- duplicated(pepTab[,"mass"], fromLast=FALSE)           # wo 1st instance
  if(any(chMa)) {                                        # redundant iso-masses exist, look for preferential cleavage
    chM2 <- chMa | duplicated(pepTab[,"mass"], fromLast=TRUE)           # all iso masses
    chM3 <- is.na(pepTab[which(chM2),"ambig"])                   # check for lines to mark as isoMass (sam masse but not yet marked as 'duplSequence')
    if(any(chM3)) pepTab[which(chM2)[which(chM3)],"ambig"] <- ambigLab[2]  #"isoMass"
    chM5 <- which(chM2 & pepTab[,"ty"] =="inter")                # consider for preferent cleavage (same mass & internal fragm), still too much due to duplicated seq
    if(debug) {message(fxNa," .. mFr4"); mFr4 <- list(protTab=protTab,minFragSize=minFragSize,maxFragSize=maxFragSize,internFra=internFra,knownMods=knownMods,massTy=massTy,
      pep2=pep2,pepTab=pepTab,taoPo=taiPo,chTa=chTa,chM2=chM2,chM5=chM5)}
    if(length(chM5) >0) {
      pTa <- pepTab[chM5,c("no","origNa","seq","precAA","tailAA","mass","beg")] #c("no","origNa","seq","precAA","tailAA","beg","end","mass")
      pTa <- wrMisc::sortBy2CategorAnd1IntCol(pTa, categCol=c("origNa","mass"),numCol="beg", findNeighb=TRUE, decreasing=FALSE, callFrom=fxNa)            # ad col "neiGr"
      chNA <- is.na(pTa[,"neiGr"])
      if(debug) {message(fxNa," .. mFr5"); mFr5 <- list(protTab=protTab,minFragSize=minFragSize,maxFragSize=maxFragSize,internFra=internFra,knownMods=knownMods,massTy=massTy,
      pep2=pep2,pepTab=pepTab,taoPo=taiPo,chTa=chTa,chM2=chM2,chM5=chM5, pTa=pTa,chNA=chNA,remNonConfPrefFragm=remNonConfPrefFragm,prefFragPat=prefFragPat)}
      ## remove fragments not expected due to preferential fragmentation sites
      if(!all(chNA) && remNonConfPrefFragm) {                     # lines to inspect exist
        if(any(chNA)) pTa <- pTa[which(!chNA),]
        badLi <- as.integer(unlist(by(pTa, pTa[,"neiGr"], function(y) {y <- as.matrix(y); if(nrow(y) >1) evalIsoFragm(y, prefFragPat=prefFragPat, callFrom=fxNa)})))
        if(length(badLi) >0) { pepTab <- pepTab[-1*match(badLi,pepTab[,"no"]),]
          if(!silent) message(fxNa,"Due to preferential fragmentation sites discard ",length(badLi)," fragments,  ",nrow(pepTab)," remain")}}}} 
  pepTab <- cbind(pepTab, modSpec=rep("",nrow(pepTab)))       # needed for documenting specific modifications in .singleSpecModif()
  pepTab }
   

#' Reorganize List Of Peptide Fragments To Matrix 
#'
#' This function allows reorganiziong a list (of lists) of peotide fragments into matrix
#' 
#' @param x (list) list of lists with charcter vectors of sequences with names that can be parsed eg 'x.1-7' to extract 'beg'&'end' otherwise ALL output will be NA (+message form extractLast2numericParts())
#' @param subLiNames (character) 
#' @param inclNo (logical) add 1st col with number
#' @param fullSeq (character) to reinject full sequence which may not be used in names of 'x' and not be in x[[1]][["full"]]
#' @param outCol (character) columns to create in output 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) for bug-tracking: more/enhanced messages
#' @return This function returns matrix with fragment sequence, mass, start- and end-position, heading and tailing AA (or NA if terminal fragment)
#' @seealso \code{\link{makeFragments}};  \code{\link{evalIsoFragm}}, from package \href{https://CRAN.R-project.org/package=wrProteo}{wrProteo} \code{\link[wrProteo]{convAASeq2mass}}, \code{\link[wrProteo]{AAmass}}, \code{\link[wrProteo]{massDeFormula}}
#' @examples
#' prot1 <- c(protP="KEPTIDE", pro2="MPRATE")
#' ## fragment all target proteins
#' pep3 <- lapply(prot1, fragmentSeq, minSize=3, maxSize=5, internFragments=FALSE, 
#'   separTerm=TRUE, keepRedSeqs=TRUE)
#' pepTab <- .exNamesTyDeList(pep3, fullSeq=prot1) 
#' @export
.exNamesTyDeList <- function(x,subLiNames=c("full","Nter","Cter","inter"), inclNo=TRUE, fullSeq=NULL,
  outCol=c("seq","orig","origNa","ty","seqNa","beg","end","precAA","tailAA","ambig","mass"), silent=FALSE, callFrom=NULL, debug=FALSE) {     #"ambigTy"
  ## function to extract all information from pep2 (list of lists with peptides) & organipredMae by groups in output as matrix (all full, all Nter,...)
  ## 'x' .. list of lists with charcter vectors of sequences with names that can be parsed eg 'x.1-7' to extract 'beg'&'end' otherwise ALL output will be NA (+message form extractLast2numericParts())
  ## 'inclNo' .. add 1st col with number  
  ## 'fullSeq' .. to reinject full sequence which may not be used in names of 'x' and not be in x[[1]][["full"]]
  ## NOTE : cols 'precAA' & 'tailAA' won't be filled since orig (parent) sequence for fragments not known with input of function
  ##   similar for cols 'ambig' & 'mass'
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".exNamesTyDeList")
  chLe <- sapply(x,length) <1
  if(any(chLe)) {if(all(chLe)) stop("'x' is empty !") else x <- x[which(!chLe)]}
  out <- matrix(NA, nrow=sum(sapply(x,function(y) sum(sapply(y,length)))), ncol=length(outCol), dimnames=list(NULL,outCol))
  iniNa <- names(x)
  names(x) <- NULL
  fullSequ <- unlist(sapply(x,function(y) y$full))
  if(length(fullSequ) < length(x)) fullSequ <- if(is.null(fullSeq)) iniNa else fullSeq      # try to find read full sequence, otherwise return to names of 'x'
  predMa <- 1
  subLiNames <- subLiNames[subLiNames %in% unlist(lapply(x,names))]
  for(i in 1:length(subLiNames)) {
    tm <- lapply(x,function(x) x[[subLiNames[i]]])
    chLe <- sapply(tm,length) >0 
    if(any(chLe)) { 
      if(any(!chLe)) tm <- tm[which(chLe)]
      if(length(tm) >0) { 
        protNa <- rep(iniNa[which(chLe)], sapply(tm, length))
        fullSe <- rep(fullSequ[which(chLe)], sapply(tm, length))
        tm <- unlist(tm)
        tm <- cbind(seq=tm, orig=fullSe, origNa=protNa, ty=rep(subLiNames[i],length(tm)), seqNa=names(tm), wrMisc::extractLast2numericParts(names(tm)))
        colnames(tm)[ncol(tm)+(-1:0)] <- c("beg","end")
        addCol <- (!outCol %in% colnames(tm))
        if(sum(addCol) >0) tm <- cbind(tm, matrix(NA, nrow=nrow(tm), ncol=sum(addCol), dimnames=list(NULL,outCol[which(addCol)])))
        if(!identical(colnames(tm),outCol)) {tm <- tm[,wrMisc::naOmit(match(outCol, colnames(tm)))]
          if(!silent) message(fxNa,": reduce cols to match argument 'outCol'")}
        out[predMa:(predMa +nrow(tm)-1),1:ncol(tm)] <- tm
        predMa <- predMa + nrow(tm) }}}
  if(inclNo) cbind(no=as.character(1:nrow(out)),out) else out }
       
