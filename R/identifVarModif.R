#' Idenitfy Variable Modifications
#'
#' Take result from \code{identifFixedModif} and search for variable modifications (only on identified fixed modif), ie 2nd step for identif of var modifs.
#' To reduce the complexity of the search space, only peptide fragments identified with fixed identifiactions will be considered for possible variable modifications.
#'
#' @details The main matching results are in output$massMatch : This list has one entry for each predicted mass where some matches were found. 
#' Thus, the names of the list-elements design the index from argument \code{expMass}. 
#' Each list-element contains a numeric vector giving the difference observed to predicted, the names design the unique predicted peptide index/number from output$preMa[,"no"]
#' 
#' @param zz (list) min input, result from \code{identifFixedModif}, must conatain elements 'nmassMatch','preMa','pepTab','recalibFact','recalibData'
#' @param modTy (character) type of fixed and variable modifications
#' @param expMa (matrix) experimental m/z values
#' @param maxMod (integer) maximum number of residue modifications to be consiered in fragments (values >1 will increase complexity and RAM consumption)  
#' @param identMeas (character) comparison type (used in findCloseMatch(), default ="ppm"), used with limit 'limitIdent'
#' @param knownMods (character) optional custom alternative to \code{AAfragSettings(ou="all")$knownMods}  
#' @param limitIdent (integer) limit applied to 'identMeas' 
#' @param filtAmbiguous (logical) toggle to remove all ambiguous identifications
#' @param indexStart (integer) for keeping correct index at iterative use
#' @param recalibFact (numeric, length=1) 
#' @param suplPepTab (matrix) predicted fragments (incl fixed and var modifs) to include to search (allowong to ensure overlap to include hits close to prev search)
#' @param massTy (character) 'mono' or 'average'
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return list with $massMatch (list of exerimental peptides matching to one or more predicted), $preMa (predicted ions, including fixed and variable modif),
#'  $pepTab (predicted neutral peptides, wo modifications), $expMa (experimental mass from input), $recalibFact (recalibration factor as from input), $docTi (time for calculations)
#' @seealso \code{\link{makeFragments}}, \code{\link{identifFixedModif}}, \code{\link{identifyPepFragments}}
#' @examples
#' protP <- c(protP="PEPTIDE")
#' obsMassX <- cbind(a=c(199.1077,296.1605,397.2082,510.2922,625.3192),
#'   b=c(227.1026,324.1554,425.2031,538.2871,653.3141),
#'   x=c(729.2937,600.2511,503.1984,402.1507,289.0666),
#'   y=c(703.3145,574.2719,477.2191,376.1714,263.0874))
#' rownames(obsMassX) <- c("E","P","T","I","D")      # all 1 & 7 ions not included
#' identP10 <- identifFixedModif(prot=protP,expMass=as.numeric(obsMassX),minFragSize=2, 
#'   maxFragSize=7,modTy=list(basMod=c("b","y")))     # looks ok
#' identP10v <- identifVarModif(identP10,list(varMod="h"), as.numeric(obsMassX),2)
#' identP10v$massMatch                    # list of matches
#'
#' @export
identifVarModif <- function(zz, modTy, expMa, maxMod, identMeas="ppm", knownMods=NULL, limitIdent=5, filtAmbiguous=FALSE, indexStart=1,
  recalibFact=NULL, suplPepTab=NULL, massTy="mono", silent=FALSE, callFrom=NULL, debug=TRUE){
  ## take result from  identifFixedModif() ie 2nd step for identif of var modifs
  ## return list with $massMatch for matches
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="identifVarModif")
  docTi <- rep(NA,5)
  names(docTi) <- c("ini_identifVarModif","countPotModifAAs","addMassModif","finUniqCheck","findCloseMatch")
  docTi[1] <- Sys.time()                                                            # makeFragments() consumes 95-99% of time !!
  recalibFact <- if("recalibFact" %in% names(zz)) as.numeric(zz$recalibFact) else 0
  ## adjust index in case of looping over pep-lengths
  if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange00"); xxidentifVarInAASizeRange00 <- list(zz=zz,preMa=zz$preMa,modTy=modTy,maxMod=maxMod,expMa=expMa,filtAmbiguous=filtAmbiguous,recalibFact=recalibFact,suplPepTab=suplPepTab)}   # cou=cou,
  modTv <- checkModTy(modTy, knownMods=knownMods, silent=silent, callFrom=fxNa)   # after checking modTv is already set as list of 3 , count based on $varMo2
  modTv$basMod <- "" 
  ## now address variable modifications ..
  chmodTv <- if(length(modTy$varMod) <1) FALSE else modTy$varMod !=""
  if(any(chmodTv) && length(zz$massMatch) <1) { 
    if(!silent) message(fxNa,"No hits from 1st pass or no var modif selected, can't run var modif")
    chmodTv[which(chmodTv)] <- FALSE }
  if(any(chmodTv) && length(zz$massMatch) >0) {                                                  # (if any) variable modifications : modTy not empty, 2nd pass
        ## problem : names $massMatch have heading 'x', suplPepTab 
    if(!silent) message(fxNa,"Enter 2nd pass (for variable modif)")
    if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange0"); xxidentifVarInAASizeRange0 <- list(zz=zz,modTy=modTy,modTv=modTv,maxMod=maxMod,expMa=expMa,chmodTv=chmodTv,filtAmbiguous=filtAmbiguous,preMa=zz$preMa,suplPepTab=suplPepTab)}   # cou=cou,
    nPepCand <- nrow(zz$preMa)                     # number of peptides in 1st pass (preMa)
    ### need to select those (unmodified) peptides where matches found -> check for their variable modifications in 2nd pass comparing to experimental
    usePred <- as.numeric(sub("^y","",sub("^x","", names(zz$massMatch))))
    usePred <- match(usePred, zz$preMa[,"no"])
    chNa <- is.na(usePred)
    if(all(chNa)) message("\n",fxNa," PROBLEM : none of previous massMatch indexes/hits found in $zz$preMa !!!")
    if(any(chNa)) { usePred <- usePred[which(!chNa)] 
      if(!silent) message(fxNa,"NOTE:   ",sum(chNa)," previous massMatch indexes/hits NOT found in $zz$preMa")}   
    zz$pepTab <- cbind(zz$pepTab,modif="")
    rownames(zz$pepTab) <- zz$pepTab[,"seqNa"]                     # had no rownames so far
    ## prepare 'pepTab'-like matrix using preMa already including fixed modifs
    pepTab2 <- zz$preMa[usePred,wrMisc::naOmit(match(c(colnames(zz$pepTab),"finMass","mod"),colnames(zz$preMa)))]     # use/ready for 2nd pass
  	if(length(dim(pepTab2)) <2) pepTab2 <- matrix(pepTab2, nrow=length(usePred), dimnames=list(rownames(zz$preMa)[usePred],names(pepTab2)))
    colnames(pepTab2)[match("finMass",colnames(pepTab2))] <- "mass"
    ## note : pepTab2[,"origNa"] does hold duplicated names at this level
    pepT2Na <- if(identical(as.character(zz$preMa[,"no"]),zz$pepTab[,"no"])) zz$preMa[,"seqNa"] else paste0(pepTab2[,"origNa"],".",pepTab2[,"beg"],"-",pepTab2[,"end"])    # zz$pepTab : no rownames, zz$preMa rownames with fixedMod & specMod
    chSeqNa <- colnames(pepTab2) %in% "seqNa"
    if(any(chSeqNa)) pepTab2[,which(chSeqNa)] <- pepT2Na else pepTab2 <- cbind(pepTab2,seqNa=pepT2Na)        
    couV <- countPotModifAAs(pepTab=pepTab2, maxMod=maxMod, modTy=modTv, silent=silent,debug=debug,callFrom=fxNa)   # max number of modifs gets corrected here
    ## potential problem : peptides shared betw different input proteins may appear multiple times after countPotModifAAs() ...
    if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange1"); xxidentifVarInAASizeRange1 <- list(zz=zz,modTv=modTv,couV=couV,pepTab2=pepTab2,maxMod=maxMod,expMa=expMa,chmodTv=chmodTv,filtAmbiguous=filtAmbiguous,modTy=modTy,suplPepTab=suplPepTab)}   # cou=cou,
    docTi[2] <- Sys.time()                                                            # after countPotModifAAs()
    lastIndex <- if("no" %in% colnames(zz$preMa)) max(as.integer(zz$preMa[,"no"])) else nrow(zz$preMa)
    ## make table with mass-modifications
    preMaV <- addMassModif(cou=couV$cou, pepTab=pepTab2, lastIndex=lastIndex, combTerm=couV$combTerm, modTy=modTv,silent=silent,debug=debug,callFrom=fxNa)    # list of new masses
    ## checks 9sep19: preMaV$pepTab  is already appended and has corerect masses !    
    docTi[3] <- Sys.time()                                                            # after addMassModif()
    ## correct (all) var modif peptides beeing one H too heavy at final testing
    #mass already ok#preMaV$pepTab[,"mass"] <- as.numeric(preMaV$pepTab[,"mass"]) -wrProteo::.atomicMasses()["H",massTy] 
    if(!silent) message(fxNa," addMassModif/2 done ")
    preMaV$basMod <- NULL   # should not contain var modifs
    ## now check if combining full matrix of fixed modif with selected var modif for final comparison needed
    if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange1b"); xxidentifVarInAASizeRange1b <- list(preMaV=preMaV,zz=zz,modTv=modTv,couV=couV,pepTab2=pepTab2,maxMod=maxMod,expMa=expMa,chmodTv=chmodTv,filtAmbiguous=filtAmbiguous,modTy=modTy,suplPepTab=suplPepTab,lastIndex=lastIndex)}     
    
    ## why try to append more peptides from pepTab ??
    preMaV <- preMaV$pepTab
    
    docTi[4] <- Sys.time()                                                            # after  finUniqCheck()           
    if(!silent) message(fxNa,"Ready to 2nd pass testing ",nrow(preMaV)," (including ",sum(!is.na(preMaV[,"ambig"]))," ambiguous : ",
      wrMisc::pasteC(unique(wrMisc::naOmit(preMaV[,"ambig"]))),") predicted against ",length(expMa)," input masses")
    ## (idea:) option include filter to remove all so far ambiguous ??
    if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange2"); xxidentifVarInAASizeRange2 <- list(preMaV=preMaV,zz=zz,pepTab2=pepTab2,modTy=modTy,modTv=modTv,usePred=usePred,expMa=expMa,identMeas=identMeas,filtAmbiguous=filtAmbiguous,suplPepTab=suplPepTab,couV=couV,maxMod=maxMod,recalibFact=recalibFact) }
    ## add 'suplPepTab' if overlapping search, last 'complete' table of fragments before 2nd run of matching 
    if(length(suplPepTab) >0) { 
      colnames(suplPepTab)[which(colnames(suplPepTab)=="mass")] <- "mass"
      suplPepTab <- cbind(suplPepTab, seqNa=rownames(suplPepTab), modif=NA, ppmClos=NA,bdifClos=NA,bindClos=NA)
      zz$preMa <- rbind(zz$preMa, suplPepTab[,match(colnames(zz$preMa),colnames(suplPepTab))])  }
    if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange3"); xxidentifVarInAASizeRange3 <- list(preMaV=preMaV,zz=zz,pepTab2=pepTab2,modTy=modTy,modTv=modTv,usePred=usePred,
      suplPepTab=suplPepTab,expMa=expMa,identMeas=identMeas,filtAmbiguous=filtAmbiguous,suplPepTab=suplPepTab,couV=couV,maxMod=maxMod,recalibFact=recalibFact) }
    ## now (new) run of mass comparison with initial + PTM (ie varMod) enriched
    ## recalibration will be applied here to experimental masses
    predMa <- as.numeric(preMaV[,"mass"])
    names(predMa) <- preMaV[,"no"]
    expMass <- as.numeric(expMa) +recalibFact
    names(expMass) <- if(length(names(expMass)) >0) names(expMass) else 1:length(expMa) 
    preMaRa <- range(predMa,na.rm=TRUE) +c(-1,1)              
    chExpM <- expMass > preMaRa[1] & expMass < preMaRa[2]
    if(any(!chExpM)) {expMass <- if(all(!chExpM)) NULL else expMass[which(chExpM)]}
    if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange3b"); xxidentifVarInAASizeRange3b <- list(preMaV=preMaV,zz=zz,pepTab2=pepTab2,modTy=modTy,modTv=modTv,usePred=usePred,
      predMa=predMa,expMass=expMass,chExpM=chExpM,limitIdent=limitIdent,identMeas=identMeas,
      suplPepTab=suplPepTab,expMa=expMa,filtAmbiguous=filtAmbiguous,suplPepTab=suplPepTab,couV=couV,maxMod=maxMod,recalibFact=recalibFact) }
    massMatch1 <- wrMisc::findCloseMatch(x=predMa, y=expMass, compTy=identMeas, limit=limitIdent, sortMatch=FALSE,callFrom=fxNa)      ## FINAL selection of suitable hits (within limits)
    if(length(massMatch1) >0) names(massMatch1) <- sub("^x","", names(massMatch1))                #not needed any more
    if(debug) {message(fxNa, " .. xxidentifVarInAASizeRange4");  xxidentifVarInAASizeRange4 <- list(massMatch1=massMatch1,preMaV=preMaV,zz=zz,predMa=predMa,pepTab=zz$pepTab,pepTab2=pepTab2,modTy=modTy,modTv=modTv,usePred=usePred,expMa=expMa,identMeas=identMeas,limitIdent=limitIdent,suplPepTab=suplPepTab)}    #toPreMa=toPreMa,origMassInf=origMassInf,cou=cou
    } else {preMaV <- NULL; massMatch1 <- zz$massMatch}                 # end of variable modifications     
  docTi[5] <- Sys.time()                      # (final) after findCloseMatch()
  ## note: zz$pepTab remains unchanged (since wo var modif)
  list(massMatch=massMatch1, preMa=preMaV, pepTab=zz$pepTab, expMa=expMa, recalibFact=recalibFact, docTi=docTi) }
       
