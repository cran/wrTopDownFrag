#' Identify terminal and internal protein/peptide-fragments as matches to experimental MS-peaks
#'
#' Function for predicting internal and terminal peptide-fragments and compare them with experimental monoisotopic masses.    
#' The accuracy of results is given in ppm and a false discovery rate (FDR) for the identification is estimated.
#' The identifed fragments are also checked for preferential break sites, a score including this and other parameters is given with the results.
#' 
#' @param expMass (matrix or data.frame) 
#' @param pep (character) protein/peptide sequences to be used for fragmentation
#' @param modTy (list) defining fixed and variable modifications
#' @param minFragSize (integer) min length in AA of peptides to be considered (please see you spectrometers characteristics)
#' @param maxFragSize (integer) max length in AA of peptides to be considered (please see you spectrometers characteristics)
#' @param identMeas (character) comparison type (used in findCloseMatch(), default ="ppm"), used with limit 'limitIdent'
#' @param limitIdent (integer) limit applied to 'identMeas' 
#' @param internFra (logical) switch from including all internal fragments to terminal fragments only (if F)
#' @param specModif (list) optional custom single-site modifications (eg ions bound), will be processed using \code{.singleSpecModif}
#' @param massTy (list) list of modifications/fragmentation-type(s) to consider, organipredMae as 'basMod' (any occurance) and 'varMod' (optional aoccurance), 'modPos' (position of modif, integer), 'modMass' (mass to be added), 'modName' (name), 'modFixed' (fixed or variable modif, logical)
#' @param chargeCatchFilter (logical) filter (upfront) to consider only peptides containing AAs capable of catching extra charges (K, R, H, defined via \code{.chargeCatchingAA()})  
#' @param corMutShift (numeric) (numeric) vector of decoy-type possible mass shifts (eg from load("C:/E/projects/MassSpec/fragmIdentif/corMutShift.RData"))
#' @param nProc (integer)  number of preocessors to use
#' @param parallDefault (logical)  if 'parallDefault'=F no multiprocessor parameters set for BiocParallel
#' @param multParam (list)
#' @param maxMod (integer) maximum number of residue modifications to be consiered in fragments (values >1 will increase complexity and RAM consumption)  
#' @param recalibrate (logical) recalibrate based on region with highest density of experim values 
#' @param filtAmbiguous (logical)
#' @param prefFragmPat (matrix) optional custum preferential fragmentation pattern (otherwise \code{.prefFragPattern()} will be used)
#' @param sortOutputByMass (logical)
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return matrix of idenitfied ions
#' @seealso \code{\link{makeFragments}}, \code{\link{identifVarModif}}, \code{\link{identifFixedModif}}, \code{\link[wrMisc]{findCloseMatch}}, \code{scoreProteinFragments}                                                                                                         
#' @examples
#' protP <- c(protP="PEPTIDE")
#' obsMassX <- cbind(a=c(199.1077,296.1605,397.2082,510.2922,625.3192),
#'   b=c(227.1026,324.1554,425.2031,538.2871,653.3141),
#'   x=c(729.2937,600.2511,503.1984,402.1507,289.0666),
#'   y=c(703.3145,574.2719,477.2191,376.1714,263.0874))
#' rownames(obsMassX) <- c("E","P","T","I","D")      # all 1 & 7 ions not included
#' modTy1 <- list(basMod=c("b","y"), varMod=c("p","o","q"))
#' frag1 <- identifyPepFragments(ex=as.numeric(obsMassX), pe=protP, modTy=modTy1, 
#'   minFragSize=2, chargeCatchFilter=FALSE)
#' (frag1b <- if(length(unlist(frag1$identif)) >0) modifFragmTabOutput(frag1))
#' @export
identifyPepFragments <- function(expMass, pep, modTy=NULL, minFragSize=6, maxFragSize=75, identMeas="ppm", limitIdent=5, internFra=TRUE, specModif=NULL, massTy="mono",chargeCatchFilter=TRUE, # ambigPref="",
  corMutShift=NULL, nProc=1, parallDefault=TRUE, multParam=NULL, maxMod=c(p=3,h=1,k=1,o=1,m=1,n=1,u=1,r=1,s=1), recalibrate=TRUE, filtAmbiguous=FALSE, prefFragmPat=NULL,
  sortOutputByMass=FALSE, silent=FALSE, callFrom=NULL, debug=FALSE){                                                                                                               
  ##  version '3'
  ## Identify protein/peptide fragments by monoisotopic m/predMa mass given in 'expMass'
  ## return matrix for all (potentially identified proteins/peptides)
  ##  details : predict mass of fragmented peptides (from 'pep', character vector) wo variable modifications, for all good hits include variable modifications
  ## 'expMass'.. experimental mass
  ## 'pep' .. protein/peptide sequences to be used for fragmentation
  ## 'minFragSize','maxFragSize'  ..
  ## ''  ..
  ## 'identMeas' .. comparison type (used in wrMisc::findCloseMatch(), default ="ppm")
  ## 'limitIdent' .. threshold for identification (used in findCloseMatch(), default =5)
  ## 'modTy' .. list of modifications/fragmentation-type(s) to consider, organipredMae as 'basMod' (any occurance) and 'varMod' (optional aoccurance)
  ##   'p', 'q' .. gain/loss Phosphorlyation; loss of H2PO4 (S,T)  (S,T,Y+PO4); set 'maxMod["p"]' to 0 for no dephosphorylation
  ##   'i' .. loss of CO (terminal)   
  ##   'd' ... terminal loss of H20 (not internal ions)  
  ##   'h' ... loss of water (S,T,E,D)   ??also in (E,D)??   (HL: less frequent in S & T)
  ##   'n','k' .. n: loss of N-terminal ammonia (at Q); 'k': full ammonia loss (Q,K,R,N), ie if 'k' then 'n' can be omitted
  ##   'u' .. loss of Urea (C-termial R only)
  ##   'o', 'r' .. o: oxidation at M;  r: oxidation at C
  ##   's' .. acetylation at S
  ##   'm' .. methylation (at R,K)
  ##  what is 'f', 'j' ?? (-H20, -CO  both at intern fragm)  ?
  ##
  ## 'internFra' .. switch from including all internal fragments to terminal fragments only (if F)
  ## 'massTy' .. to work in residue weight (ie -H20, otherwise molecule weight)
  ## 'ambigPref' .. prefix to add to names of non-unique (ie ambiguous) masses
  ## 'corMutShift' .. (numeric) vector of decoy-type possible mass shifts (eg from load("C:/E/projects/MassSpec/fragmIdentif/corMutShift.RData"))
  ## 'nProc','parallDefault' .. number of preocessors to use; if 'parallDefault'=F no multiprocessor parameters set for BiocParallel()
  ## 'ppm' .. if numeric & length=1, for additional check/column 'ambig.ppm' showing if values in column 'finMass' separate within given range (in ppm units)
  ## 'maxMod' .. max number of residue specific modifications
  ##   not used so far : (e f g) j l t u   (s prepared)
  ## 'specModif' .. special modification at single AA location : provide list with 'modOrigin' (seq of prot to be modif),
  ##     'modPos' (position of modif, integer), 'modMass' (mass to be added), 'modName' (name), 'modFixed' (fixed or variable modif, logical)
  ## 'recalibrate' ..(logical) to recalibrate based on region with highest density of experim values (at min intensity ??)
  ## 'prefFragmPat' .. optional custum preferential fragmentation pattern (otherwise .prefFragPattern() will be used)
  ## note/idea - speed improvements : rather treat as tree of sequences  (instead of indiv =independent)! (eg for .convAASeq2mass)
  ## note : so far no speed improvement when going 3 procs with .parCombinateAllAndSum()
  ## note : sequence parts identical between various 'pep' : will be marked in column 'ambig' as 'duplSequence', but so far not specifying which 'pep' (need to evolve get1stOfRepeatedByCol())
  #potIsoFragm=c("D","E","P")                           #
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="identifyPepFragments")
  message("++ Note of caution : This function might not find all matches, the next version  of this package will fix potential problems ++")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  maxP <- 2             # max Phospho considered
  lossOfP <- TRUE          # also include loss of P
  #aaSegments  loopOverIdent
  overl <- 0          # overlap to prev aa length
  sinMa <- 20         # max range of le to pred in single run
  nRandRepeat <- 25                       # number of time FDR estimation gets repeated
  docTi <- c(Sys.time(), rep(NA,5))               # ini ;
  iPF0g <- NULL;
  names(docTi) <- c("ini","main identif","formatMatrix","scoring")
  docTi2 <- matrix(ncol=14, nrow=1, dimnames=list(1,c("ini_identifFixedModif","makeFragments","countPotModifAAs","addMassModif","finUniqCheck","findCloseMatch","recalib",
     "ini_identifVarModif","countPotModifAAs","addMassModif","finUniqCheck","findCloseMatch", "randData","FDR")))   # 7+5+2
  if(length(expMass) <1) stop(fxNa,"Argument 'expMass' seems empty -> nothing to compare with !!")
  if(length(pep) <1) stop(fxNa,"Argument 'pep' seems empty -> no sequence to fragment !!")
  ## NOTE : new modifications must be entered in knownMods, specAAmods (for molec compos); 'f','m','c'  used for designating full-length/Nterm/Cterm
  knownMods <- AAfragSettings(outTy="all")
  specAAMod <- knownMods$specAAMod
  modChem <- knownMods$modChem
  knownMods <- knownMods$knownMods
  # check modChem ?
  .countFirstChar <- function(char,x) as.numeric(substring(x, 1, 1) %in% char)              # count occurances as  1st character of 'x'
  .countLastChar <- function(char,x) as.numeric(substring(x, nchar(x), nchar(x)) %in% char) # count occurances as last character of 'x'
  if(debug) silent <- FALSE
  setSeedFDR <- NULL
 ## checking
  chModNa <- all(names(specAAMod) %in% unlist(knownMods))                            # just to check
  if(any(!chModNa)) stop(fxNa,"Modification types with unknown corresponding AA found !")
  AAmods <- list()                                   #
  if(!is.list(modTy)) {modTy <- list(basMod=modTy)
    if(!silent) message(fxNa,"'modTy' not list, suppose 'basMod' only")}
  modTy <- checkModTy(modTy, knownMods=knownMods, silent=silent, callFrom=fxNa)
  if(debug) {message(fxNa, " .. iPF0"); iPF0 <- list(modTy=modTy,minFragSize=minFragSize,expMass=expMass,
    internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep)}
  ## prepare experimental input
  if(length(dim(expMass)) <2) {
    naCh <- is.na(expMass)
    if(!silent) message(fxNa,sum(naCh)," out of ",length(expMass)," NAs omitted from 'expMass'")
    if(any(naCh)) expMass <- expMass[which(!naCh)]                                                             # remove NAs from 1st col
    if(length(names(expMass)) <1) names(expMass) <- 1:length(expMass)
    origMassInf <- as.matrix(expMass)    
  } else {                                                                      # split supplemental information (eg peak heights) from m/predMa measures
    naCh <- is.na(expMass[,1])
    if(debug) message(fxNa, "iPF0a .. dim expMass",nrow(expMass)," ",ncol(expMass))
    if(any(naCh)) {                                                             # remove NAs from 1st col
      if(!silent) message(fxNa,sum(naCh)," out of ",nrow(expMass)," NA-lines omitted from 'expMass'")
      expMass <- expMass[which(!naCh),]}
    origMassInf <- expMass
    rowNa <- rownames(expMass)
     ## idea: make rowNa unique ?
    if(length(rowNa) <1) rowNa <- 1:nrow(expMass)                               # if no names -> index
    expMass <- as.numeric(expMass[,1])
    names(expMass) <- rowNa
  }   
  ## check for duplicated experimental values (now numeric vector)
  chDuEx <- duplicated(expMass, fromLast=FALSE)
  if(any(chDuEx)) { if(!silent) message(fxNa,"Note: Making experimental data non-redundant (",sum(chDuEx)," replicated/redundant, keep only 1st)")
    expMass <- expMass[which(!chDuEx)] 
    origMassInf <- origMassInf[which(!chDuEx), , drop=FALSE] }   
  ## prepare pep, ...
  iniPepLe <- length(pep)
  pep1 <- gsub("[[:digit:]]|[[:lower:]]|[[:blank:]]|[[:punct:]]","", pep)                  # keep only caps letters of peptide input
  if(any(nchar(pep1) <1)) stop("No (single letter caps code-) AA found in 'pep'")
  if(is.null(names(pep1))) names(pep1) <- pep1 else {          # if names not given for 'pep' use sequence
    noNa <- nchar(names(pep1)) <1
    if(any(noNa)) names(pep1)[which(noNa)] <- pep1[which(noNa)]}
  if(debug) {message(fxNa, " .. iPF0b"); iPF0b <- list(modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,expMass=expMass,
    internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf,specModif=specModif,
    maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous)}
  ##
  ## run LOOP over different mass-ranges : easier to loop over aa-bandwidth
  runs <- 0:floor((min(max(nchar(pep)), maxFragSize) -minFragSize) /(sinMa-overl))
  startL <- minFragSize +(sinMa*(runs) -overl*(runs)) 
  endL <- startL +sinMa -1
  if(any(endL > maxFragSize)) endL[which(endL > maxFragSize)] <- maxFragSize
  runs <- cbind(startL, endL, startR=c(startL[1],endL[-length(endL)]+1), recalibrate=0, recalibFact=0, nPred=NA, prevMaxIndex=NA, curMaxInd=NA,
    minMass=NA ,maxMass=NA, minMassDec=NA, nUniqPred=0, nIdentif=NA, nRandHit=NA, FDR=NA, nSupl=NA)
  chRu <- runs[,3] >= runs[,2]
  if(any(chRu)) {runs <- if(sum(!chRu) >1) runs[which(!chRu),] else matrix(runs[which(!chRu),], nrow=1, dimnames=list(NULL,colnames(runs)))}
  #docTi[1] <- Sys.time()
  if(nrow(runs) >1) docTi2 <- matrix(NA, ncol=ncol(docTi2), nrow=nrow(runs), dimnames=list(NULL,colnames(docTi2)))
  ## check which run will be used for recalibration
  if(recalibrate) {  
    ## need to figure out which loop run 1st for recalibrating (should be run with high density of experim values)
    ## rougthly cut experim values by estimated borders of loop-segments (aa-length at avermass for resutant peptides based on human frq of aa usage)  
    if(nrow(runs) >1) {
      avMolWe <- 118.9                        # aver mass for aa
      avMolWe <- c(runs[1,1], runs[,2]) *avMolWe
      chReCalRu <- table(as.numeric(cut(expMass, breaks=avMolWe)))
      recalRun <- which(chReCalRu==max(chReCalRu, na.rm=TRUE))
      if(!silent) message(fxNa," RECALIBRATION was set to subset ",paste(runs[recalRun,1:2],collapse="-")," aa length (set no ",recalRun,", apr ",100*round(chReCalRu[recalRun]/length(expMass),2),"% of experim values) ")
      if(recalRun != 1) {runs <- runs[c(recalRun, (1:nrow(runs))),]               # repeat run for recalibration as 1st
        docTi2 <- rbind(docTi2, NA)}
      runs[1,"recalibrate"] <- 1}         
  } 
  if(debug) {message(fxNa, " .. iPF0c"); iPF0c <- list(modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,expMass=expMass,specModif=specModif,maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous,prefFragmPat=prefFragmPat,overl=overl,sinMa=sinMa,overl=overl,docTi2=docTi2,
    runs=runs,recalibrate=recalibrate,internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf)}
  ## loop over start&end-points as aa-length
  massMatch1 <- preMaFull <- NULL                           # collect preMa-lines where hits have been found
  suplPepTab <- NULL                          # in case of overlapping peptide-length loops: matrix for re-injecting pepTab to next round
  fullMatch <- list()                         # collect mass-matches
  massChar <- list(n=NULL, perc=NULL)
  identCount <- array(dim=c(nrow(runs),0,2))  # initialipredMae
  for(ru in 1:nrow(runs)) {  # run no
    if(!silent) message("\n+++++\n",fxNa,"  start run ",ru)
    recal <- if(runs[ru,"recalibrate"] >0) TRUE else runs[ru,"recalibFact"] 
    runs[ru,"prevMaxIndex"] <- indexLast <- if(ru <2) 0 else {if(runs[ru-1,1] >runs[ru,1]) 0 else max(as.integer(massMatch1$preMa[,1]),na.rm=TRUE) }        #runs[ru-1,"prevMaxIndex"]
    if(ru >1) if(runs[ru,"startL"] <= runs[ru-1,"startL"]) suplPepTab <- NULL   
    massMatch1 <- identifFixedModif(prot=pep1, minFragSize=runs[ru,"startR"], maxFragSize=runs[ru,"endL"],suplPepTab=NULL, indexStart=indexLast+1, internFra=internFra,chargeCatchFilter=chargeCatchFilter,
      modTy=modTy, massTy=massTy, maxMod=maxMod, expMass=expMass, specModif=specModif, knownMods=knownMods, identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous, 
      recalibrate=recal, prefFragPat=prefFragmPat, silent=silent,debug=debug,callFrom=fxNa)
       message(fxNa, " zytzr")
      allCurMz <- unique(sort(as.numeric(massMatch1$preMa[,"finMass"]))) 
    if(debug) message(fxNa," -> done identifFixedModif,  dim(runs) ",dim(runs),"  dim docTi2",wrMisc::pasteC(dim(docTi2)))
    if(debug) {message(fxNa, " .. iPF0d"); iPF0d <- list(modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,
      expMass=expMass,specModif=specModif,maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous,prefFragmPat=prefFragmPat,overl=overl,docTi2=docTi2,
      identCount=identCount,massMatch1=massMatch1, runs=runs,fullMatch=fullMatch,massChar=massChar,ru=ru,nRandRepeat=nRandRepeat,setSeedFDR=setSeedFDR,suplPepTab=suplPepTab,
      recalibrate=recalibrate,internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf)}
    docTi2[ru,1:7] <- if(length(massMatch1$massMatch) >0) massMatch1$docTi else 0               #,Sys.time()
    if(runs[ru,"recalibrate"] >0) runs[,"recalibFact"] <- massMatch1$recalibFact    # if recalibrating: set recalib fact for all runs
    if(c(runs[-1,"startL"] < runs[-nrow(runs),"startL"],FALSE)[ru]){                      # only recalib (don't use results)
      runs[ru,c("prevMaxIndex","nSupl")] <- 0
      if(!silent) message(fxNa,"Run no ",ru," only used for recalibration")        
    } else {   
      ## prepare for var modif part
      ## now 2nd part for var modif (if any)
      runs[ru,"nSupl"] <- if(length(suplPepTab) >0) nrow(suplPepTab) else 0                  # supl peptides from prev run that will be added
      ## if no initial matches : don't run identifVarModif()
      
      if(length(massMatch1$massMatch) >0 && length(modTy$varMod) >0) {
        massMatch1 <- identifVarModif(massMatch1, expMa=expMass, modTy=modTy, massTy=massTy, identMeas=identMeas, limitIdent=limitIdent,
          maxMod=maxMod, knownMods=knownMods, suplPepTab=suplPepTab, silent=silent,debug=debug,callFrom=fxNa)   
         if(debug) message(fxNa,"xxmassMatch1X .. done var Modif : new massMatch length=",length(massMatch1$massMatch)) 
        allCurMz <- unique(sort(c(allCurMz,as.numeric(massMatch1$preMa[,"mass"])))) 
        if(!silent) message(fxNa,"Run=",ru," : at varible modif : found ",length(massMatch1$massMatch)," matches to experimental masses") 
        
         if(debug) {message(fxNa, " .. iPF0d2"); iPF0d2 <- list(massMatch1=massMatch1,modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,
          expMass=expMass,specModif=specModif,maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous,prefFragmPat=prefFragmPat,overl=overl,docTi2=docTi2,
          identCount=identCount,massMatch1=massMatch1, runs=runs,fullMatch=fullMatch,massChar=massChar,ru=ru,nRandRepeat=nRandRepeat,setSeedFDR=setSeedFDR,suplPepTab=suplPepTab,
          recalibrate=recalibrate,internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf)}
      } else {
        if(!silent) message(fxNa,"Nothing to do for variable modifications ...")
        if(length(massMatch1$preMa) >0) colnames(massMatch1$preMa)[which(colnames(massMatch1$preMa)=="finMass")] <- "mass" 
      }
      if(length(massMatch1$preMa) >0) {docTi2[ru,8:12] <- if(length(massMatch1$docTi) <1 || length(massMatch1$docTi) >5) 0 else massMatch1$docTi }              #,Sys.time()
      ## check & append massMatch1 for export  (remove prev suplPepTab from current massMatch1$massMatch )  
      if(length(massMatch1$massMatch) >0) { 
        massMatch1$massMatch <- wrMisc::fuseCommonListElem(massMatch1$massMatch, removeDuplicates=TRUE)     # why do some peptide IDs appear twice ?? (after var modif pep shoud have differnt no !)
        ## add to fullMatch
        if(length(fullMatch) >0) {newEl <- length(fullMatch) +(1:length(massMatch1$massMatch))
          fullMatch[newEl] <- massMatch1$massMatch
          names(fullMatch)[newEl] <- names(massMatch1$massMatch)
          if(overl >0) fullMatch <- wrMisc::fuseCommonListElem(fullMatch)
        } else fullMatch <- massMatch1$massMatch 
        ## append results in preMa to form preMaFull
        useMaLi <- as.integer(names(massMatch1$massMatch))                       # index of predicted
        expMaNa <- sub("y","", unlist(lapply(massMatch1$massMatch, names)))                    # index of experim with hits
        usePrLi <- rep(wrMisc::naOmit(match(useMaLi, massMatch1$preMa[,"no"])), sapply(massMatch1$massMatch,length))
        tmp <- cbind(if(length(usePrLi) >1) massMatch1$preMa[usePrLi,] else {
          matrix(massMatch1$preMa[usePrLi,], nrow=length(usePrLi), dimnames=list(rownames(massMatch1$preMa)[usePrLi], colnames(massMatch1$preMa)))},
          runNo=ru, obsInd=expMaNa )
        preMaFull <- if(length(preMaFull) >0) rbind(preMaFull,tmp) else tmp
        if(!silent) message(fxNa," ru=",ru," addding ",nrow(tmp)," prediced peptides to 'preMaFull'")
      }
      if(length(massMatch1$preMa) >0) {  
        runs[ru,c("minMass","maxMass")] <- range(as.numeric(massMatch1$preMa[,"mass"]),na.rm=TRUE)
        runs[ru,"nIdentif"] <- length(massMatch1$massMatch)
        runs[ru,"nPred"] <- nrow(massMatch1$preMa)
        runs[ru,"curMaxInd"] <- max(as.integer(massMatch1$preMa[,"no"]))
      
        if(debug) {message(fxNa, " .. run",ru," .. iPF0e"); iPF0e <- list(modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,expMass=expMass,
          specModif=specModif,maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous,prefFragmPat=prefFragmPat,overl=overl,runs=runs,identCount=identCount,
          massMatch1=massMatch1,fullMatch=fullMatch,massChar=massChar,ru=ru,nRandRepeat=nRandRepeat,setSeedFDR=setSeedFDR,suplPepTab=suplPepTab,docTi2=docTi2,        #massMatchLi=massMatchLi, 
          recalibrate=recalibrate,internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf,
          preMaFull=preMaFull,sinMa=sinMa)}
          
        ## prepare/'export' overlapping peptides for next round (inject to var modif search)
        if(ru < nrow(runs) && overl >0) {
          chOv <- as.integer(massMatch1$pepTab[,"end"]) -as.integer(massMatch1$pepTab[,"beg"]) > as.integer(runs[ru,"endL"]) -overl-1   # pick longest peptides (wo var mod) ..length at max range
          if(debug) message(fxNa," .. run",ru," next run suplPepTab : add pep with length > ",as.integer(runs[ru,"endL"]) -overl )
          suplPepTab <- if(any(chOv)) massMatch1$pepTab[which(chOv),-ncol(massMatch1$pepTab)] else NULL} else suplPepTab <- NULL
        ## continue for FDR estimation  
        ## prepare for 'decoys' : experim m/predMa within range of predictions
        mpredMaMarg <- c(-0.5,0.5) + as.numeric(runs[ru,c("minMass","maxMass")]) #range(expMass,na.rm=TRUE)
        chRef <- expMass >= mpredMaMarg[1] & expMass <= mpredMaMarg[2]
        expMassRu <- if(any(chRef)) expMass[which(chRef)] else signif(mean(runs[ru,c("minMass","maxMass")],na.rm=TRUE),4)
         
        xChar1 <- c(n=length(expMassRu), minV= as.numeric(runs[ru,"minMass"]), maxV=as.numeric(runs[ru,"maxMass"]))  # not perfect of overmapping runs
        names(xChar1)[1] <- "n"                             # realy need this name for 1st position
        if(debug) {message(fxNa, " .. run",ru," .. iPF0f"); iPF0f <- list(modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,expMass=expMass,
          specModif=specModif,maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous,prefFragmPat=prefFragmPat,suplPepTab=suplPepTab,identCount=identCount,
          massMatch1=massMatch1,runs=runs,fullMatch=fullMatch,massChar=massChar,ru=ru,nRandRepeat=nRandRepeat,setSeedFDR=setSeedFDR,xChar1=xChar1,overl=overl,ru=ru,
          xChar1=xChar1,nRandRepeat=nRandRepeat,identMeas=identMeas,limitIdent=limitIdent,
          runs=runs,recalibrate=recalibrate,internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf)} #massMatchLi=massMatchLi,
        ## make decoy-like random mass (using preMa -> new massMatch1 )
        evRan <- function(expMa,predMa) {   # evaluate (by counting only) random data 'predMa' against given values (real predicted): run&count instances of closeMatch
          nThresh <- length(wrMisc::countCloseToLimits(list(0), limitIdent))              # need to know in advance how many steps of FDR will be tested 
          if(length(expMa) <2) return(rep(0,nThresh)) else {                             ## don't run neither if only 1 experim val within suitable range ...
            massMatchR <- wrMisc::findCloseMatch(x=predMa, y=expMa, compTy=identMeas, limit=limitIdent, sortMatch=FALSE, maxFitShort="10%")  
            out <- wrMisc::countCloseToLimits(if(length(massMatchR) >0) massMatchR else NA, limitIdent=limitIdent)
            if(all(is.na(out))) out <- rep(0,nThresh)
            out }}
        if(debug) message(fxNa,"Run=",ru," paramerters for 'decoy' random mass ",wrMisc::pasteC(xChar1))            ## 15mar19 BUG ???
        if(debug) {message(fxNa, ".. iPF0f2"); iPF0f2 <- list(modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,expMass=expMass,
          specModif=specModif,maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous,prefFragmPat=prefFragmPat,suplPepTab=suplPepTab,identCount=identCount,
          massMatch1=massMatch1,runs=runs,fullMatch=fullMatch,massChar=massChar,ru=ru,nRandRepeat=nRandRepeat,setSeedFDR=setSeedFDR,xChar1=xChar1,overl=overl,ru=ru,
          xChar1=xChar1,nRandRepeat=nRandRepeat,identMeas=identMeas,limitIdent=limitIdent,docTi2=docTi2,
          runs=runs,recalibrate=recalibrate,internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf)}
  
        if(as.integer(runs[ru,"nIdentif"]) >1) {        # some values were identified, now proceed to FDR estimation
          if(is.na(xChar1[2])) xChar1[2] <- floor(as.numeric(runs[ru,"minMass"]))
          ranVa <- matrix(randMassByStochastic(xChar1, nRepeat=nRandRepeat, negAvoid=FALSE, setSeed=setSeedFDR, callFrom=fxNa,silent=silent),ncol=nRandRepeat)    
          docTi2[ru,13] <- Sys.time()
          runs[ru,"minMassDec"] <- min(ranVa,na.rm=TRUE)   # useful to document ???
          
          ranCount <- t(apply(ranVa, 2, evRan, allCurMz))                ## TAKES TIME !!
          if(length(ranCount) >0 && length(dim(ranCount)) <2) ranCount <- matrix(ranCount, nrow=1, dimnames=list(NULL,names(ranCount)))
          if(length(ranCount)==ncol(ranVa) && all(ranCount==0)) {
            colNa <- names(wrMisc::countCloseToLimits(list(0), limitIdent)) 
            ranCount <- matrix(0, nrow=ncol(ranVa), ncol=length(colNa), dimnames=list(NULL,colNa))}
          if(is.null(colnames(ranCount))) colnames(ranCount) <- names(wrMisc::countCloseToLimits(list(0), limitIdent))
          chNA <- is.na(ranCount)
          if(any(chNA)) ranCount[which(chNA)] <- 0
          if(ru==1 | dim(identCount)[2]==0) identCount <- array(NA,dim=c(nrow(runs),ncol(ranCount),2), dimnames=list(paste("ru",1:nrow(runs),sep="_"), colnames(ranCount),c("ident","decoy"))) 
          identCount[ru,,2] <- if(ncol(ranCount)==dim(identCount)[2]) colMeans(ranCount) else rowMeans(ranCount)
          runs[ru,c("nUniqPred","nRandHit")] <- c(length(allCurMz),round(identCount[ru,ncol(ranCount),2],2))
          identCount[ru,,1] <- wrMisc::countCloseToLimits(massMatch1$massMatch, limitIdent=limitIdent)    # count number of hits at different ppm
        } else { identCount[ru,,] <- 0; ranVa <- NULL; ranCount <- NULL }
        docTi2[ru,14] <- Sys.time()
        if(debug) {message(fxNa, " .. iPF0g"); iPF0g <- list(modTy=modTy,minFragSize=minFragSize, maxFragSize=maxFragSize,internFra=internFra,fxNa=fxNa,expMass=expMass,
          specModif=specModif,maxMod=maxMod,identMeas=identMeas,limitIdent=limitIdent,filtAmbiguous=filtAmbiguous,prefFragmPat=prefFragmPat, preMaFull=preMaFull,
          xChar1=xChar1, ranVa=ranVa,ranCount=ranCount,suplPepTab=suplPepTab,docTi2=docTi2,docTi=docTi,allCurMz=allCurMz,indexLast=indexLast,
          massMatch1=massMatch1,runs=runs,fullMatch=fullMatch,massChar=massChar,ru=ru,nRandRepeat=nRandRepeat,setSeedFDR=setSeedFDR,identCount=identCount,ranCount=ranCount,         
          runs=runs,recalibrate=recalibrate,internFra=internFra,knownMods=knownMods,specAAMod=specAAMod,modChem=modChem,pep=pep,pep1=pep1,massTy=massTy,origMassInf=origMassInf)} 
          if(ru==1 && debug) {message(fxNa, " .. iPF0g_1"); iPF0g_1 <- iPF0g}    #
          if(ru==2 && debug) {message(fxNa, " .. iPF0g_2"); iPF0g_2 <- iPF0g}    #
          if(ru==3 && debug) {message(fxNa, " .. iPF0g_3"); iPF0g_3 <- iPF0g}    #
          if(ru==4 && debug) {message(fxNa, " .. iPF0g_4"); iPF0g_4 <- iPF0g}    #
          if(debug){message(fxNa, fxNa," .. finished run-no",ru," (incl FDR estimation for this run) +++ \n")}
        } }
    }          # end of loop ru
  docTi[2] <- Sys.time()
  chdocTi <- is.na(docTi2)
  if(any(chdocTi)) docTi2[which(chdocTi)] <- 0
  docTi2 <- rowSums(apply(docTi2,1,diff))[-7]      # 12 measures
  chdocTi <- docTi2 <0
  if(any(chdocTi)) docTi2[which(chdocTi)] <- 0     # time consumed must be >=0
  if(debug){message(fxNa, ".. xxidentifyPepFragments5c"); xxidentifyPepFragments5c <- list(runs=runs,massMatch1=massMatch1,fullMatch=fullMatch,preMaFull=preMaFull,modTy=modTy,expMass=expMass,origMassInf=origMassInf,identCount=identCount,identMeas=identMeas,pep=pep)}   #,cou=cou,preMaV=preMaV,toPreMa=toPreMa,pepTab=pepTab,pepTab2=pepTab2,modTv=modTv, usePred=usePred,
  ## exploit FDR & plot 
  if(any(identCount >0)) {  
    limPpm <- as.numeric(sub("lim_","", colnames(identCount)  ))
    FDR <- identCount[,,2]/(identCount[,,1] +identCount[,,2])
    if(length(dim(FDR)) <2) FDR <- matrix(FDR, nrow=1, dimnames=dimnames(identCount)[1:2])
    chNA <- is.na(FDR)
    if(any(chNA)) FDR[which(chNA)] <- 0
    if(all(FDR==0)) {
      if(!silent) message(fxNa,"All FDR estimations =0, probably due too short protein or too small sample tested, FDR plot won't be informative")     
    } else {
      graphics::plot( limPpm, if(length(dim(FDR)) >1) FDR[1,] else FDR,type="n", main=paste(wrMisc::pasteC(names(pep))," : FDR estimation"),xlab=identMeas,ylim=c(min(FDR,na.rm=TRUE),min(max(FDR,na.rm=TRUE),0.6)))
      graphics::abline(h=c(0.1,0.2,0.3), col=grDevices::grey(0.7), lty=3)
      if(nrow(runs) >1) {for(k in 1:(nrow(runs))) graphics::lines(limPpm,FDR[k,], type="b", col=k)} else graphics::lines(limPpm, FDR, type="b", col=1)
      graphics::legend("topleft", paste(runs[,1],"-",runs[,2],"aa"), text.col=1:nrow(runs), col=1:nrow(runs), lty=1, lwd=1, seg.len=1.2, cex=0.8, xjust=0, yjust=0.5) 
    }
     #memo# range(as.integer(preMaFull[,1]))                    # name/index of predicted  
     #memo# range(as.numeric(names(fullMatch)))                 # name/index of predicted, expect to be unique
     #memo# range(as.integer(unlist(sapply(fullMatch,names))))  # name/index of experim
    colnames(preMaFull)[1] <- "predInd"
    useLi <- match(rep(names(fullMatch), sapply(fullMatch,length)), preMaFull[,1])
    out <- if(length(useLi) >1) preMaFull[useLi,] else matrix(preMaFull[useLi,],
      ncol=ncol(preMaFull), dimnames=list(rownames(preMaFull)[useLi], colnames(preMaFull))) # get those where peptide in fullMatch
    out <- cbind(out, obsMass=expMass[as.integer(unlist(sapply(fullMatch, names)))])       # fuse experimMass
    out <- cbind(out, diff=as.numeric(out[,"obsMass"]) -as.numeric(out[,"mass"]),
      ppmToPred=wrMisc::XYToDiffPpm(as.numeric(out[,"obsMass"]), as.numeric(out[,"mass"]), nSign=5, callFrom=fxNa))       
    if(debug){message(fxNa, ".. xxidentifyPepFragments5d"); xxidentifyPepFragments5d <- list(out=out,runs=runs,FDR=FDR,massMatch1=massMatch1,fullMatch=fullMatch,preMaFull=preMaFull,modTy=modTy,expMass=expMass,origMassInf=origMassInf,identCount=identCount,identMeas=identMeas,prefFragmPat=prefFragmPat,limPpm=limPpm)}   #,cou=cou,preMaV=preMaV,toPreMa=toPreMa,pepTab=pepTab,pepTab2=pepTab2,modTv=modTv, usePred=usePred,
    finFDR <- FDR[,ncol(FDR)]
    if(any(!is.finite(finFDR))) finFDR[which(!is.finite(finFDR))] <- NA
    if(any(finFDR >1)) finFDR[which(finFDR >1)] <- 1    
    runs[,"FDR"] <- round(finFDR,4)
    out <- cbind(out, FDR=finFDR[as.integer(out[,"runNo"])])
    docTi[3] <- Sys.time()  
    ## REMOVE redundant identifications with loss of neutral group
    neutralLoss <- AAfragSettings(outTy="neutralLossOrGain")                       # c("b","d","f","h","n","k","z")             ## see AAfragSettings()$modChem
    ## remove all neutral loss prefixes, thus they end up with same prefix to fuse to pepSeq &n pepName to choose line with lowest 'diff' 
    pepNa <- paste0(out[,"origNa"],"__",out[,"seq"], gsub(paste(neutralLoss,collapse="|"),"",out[,"mod"]))
    dupIdent <- duplicated(pepNa,fromLast=TRUE)
    if(any(dupIdent))  {
      dupIden2 <- duplicated(pepNa, fromLast=FALSE)
      noRepLi <- !dupIdent & !dupIden2     
      out0 <- if(length(noRepLi) >0) out[which(noRepLi),] else out[NULL,]
      repLi <- which(!noRepLi)
      out <- cbind(rownames(out)[repLi],out[repLi,])
      out <- by(out, pepNa[repLi], function(x) {x <- as.matrix(x); abV <- abs(as.numeric(x[,"diff"])); mi <- which(abV==min(abV,na.rm=TRUE)); x[mi,] })
      if(!silent) message(fxNa,"Removing ",length(repLi)-nrow(out)," redundant neutral-loss identifications") 
      out <- matrix(unlist(out), nrow=length(out), byrow=TRUE)
      rownames(out) <- out[,1]
      out <- rbind(out0, out[,-1])
      out <- out[order(out[,"origNa"], out[,"beg"], out[,"end"], out[,"mass"]),]
    }
    if(debug){message(fxNa, ".. xxidentifyPepFragments5e"); xxidentifyPepFragments5e <- list(out=out,runs=runs,FDR=FDR,massMatch1=massMatch1,fullMatch=fullMatch,preMaFull=preMaFull,modTy=modTy,expMass=expMass,origMassInf=origMassInf,identCount=identCount,identMeas=identMeas,prefFragmPat=prefFragmPat,limPpm=limPpm)}   #,cou=cou,preMaV=preMaV,toPreMa=toPreMa,pepTab=pepTab,pepTab2=pepTab2,modTv=modTv, usePred=usePred,
    ## MAIN SCORING
    if(nrow(out) >1) out <- scoreProteinFragments(out, fragmInp=origMassInf, j=0, useCol=c("orig","seq","precAA","tailAA","beg","end","ppmToPred","obsInd","predInd","Abundance"),
      prefFragPat=prefFragmPat, returnCombined=TRUE, figDraw=TRUE, silent=silent, debug=debug, callFrom=fxNa) 
  } else {
    FDR <- NULL
    docTi[3] <- Sys.time()
    out <- matrix(NA, nrow=0, ncol=11, dimnames=list(NULL,c("fraNa","seq","orig","origNa","precAA","tailAA","beg","end","ppmToPred","obsInd","predInd")))
  } 
  docTi[4] <- Sys.time()  
  if(debug){message(fxNa, ".. xxidentifyPepFragments5f"); xxidentifyPepFragments5f <- list(out=out,runs=runs,FDR=FDR,massMatch1=massMatch1,fullMatch=fullMatch,preMaFull=preMaFull,modTy=modTy,expMass=expMass,origMassInf=origMassInf,identCount=identCount,identMeas=identMeas,prefFragmPat=prefFragmPat,pep=pep)}   #,cou=cou,preMaV=preMaV,toPreMa=toPreMa,pepTab=pepTab,pepTab2=pepTab2,modTv=modTv, usePred=usePred,
  docT4 <- docT3 <- c(total=docTi[4]-docTi[1], mainIdentf=docTi[2]-docTi[1], scoring=docTi[4]-docTi[3])
  chdocTi <- docT3 > 200
  if(any(chdocTi)) docT4[which(chdocTi)] <- docT4[which(chdocTi)]/60
  docT4 <- cbind(round(docT4,1)," sec ")
  if(any(chdocTi)) docT4[which(chdocTi),2] <- " min "
  if(!silent) message(fxNa,"Found ",nrow(out),"  out of ",length(expMass)," experimental values in ",if(all(is.na(runs[,"curMaxInd"]))) 0 else max(as.integer(runs[,"curMaxInd"]),na.rm=TRUE)," theoretical fragments in ",docT4[1,],"(",docT4[2,],"main identification, ",docT4[3,],"scoring)")
  list(identif=out, overview=runs, identCount=identCount, FDR=FDR, obsMass=origMassInf, pep=pep, modTy=modTy, specModif=specModif) }
    
