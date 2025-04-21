                                   #' Identify Fixed Modifications
#'
#' Identify peptide/protein fragments based on experimental m/z values 'expMass' for given range of aa-length.
#' Internally all possible fragments will be predicted and their mass compared to the experimental values (argument \code{expMass}).
#'
#' @details The main matching results are in output$massMatch : This list has one entry for each predicted mass where some matches were found. 
#' Thus, the names of the list-elements design the index from argument \code{expMass}. 
#' Each list-element contains a numeric vector giving the difference observed to predicted, the names design the unique predicted peptide index/number from output$preMa[,"no"]
#'
#' The main element of the output is the $massMatch -list, which is in the format of  \code{\link[wrMisc]{findCloseMatch}}.
#' Thus, the list-elements names represent the line-number of mass-predictions and the values the delta-mass and their names the position of the initial query.
#' 
#' @param prot (character) amino-acid sequene of peptide or protein
#' @param expMass (numeric) experimental masses to identify peptides from
#' @param minFragSize (integer) min number of AA residues for considering peptide fragments
#' @param maxFragSize (integer) max number of AA residues for considering peptide fragments
#' @param indexStart (integer) for starting at correct index (if not 1)
#' @param suplPepTab (matrix) additional peptides to be add to theoretical peptides
#' @param internFra (logical) decide whether internal fragments should be consiered
#' @param chargeCatchFilter (logical) by default remove all peptides not containing charge-catching (polar) AAs (K, R, H, defined via \code{.chargeCatchingAA()} )
#' @param maxMod (integer) maximum number of residue modifications to be consiered in fragments (values >1 will increase complexity and RAM consumption)  
#' @param modTy (character) type of fixed and variable modifications
#' @param specModif (list) supplemental custom fixed or variable modifications (eg Zn++ at given residue) 
#' @param knownMods (character) optional custom alternative to \code{AAfragSettings(ou="all")$knownMods}  
#' @param identMeas (character) default 'ppm' 
#' @param limitIdent (character) thershold for identification in 'identMeas' units 
#' @param filtAmbiguous (logical) allows filtering/removing ambiguous results (ie same mass peptides)
#' @param recalibrate  (logical or numeric) may be direct recalibration-factor (numeric,length=1), if 'TRUE'  fresh determination of 'recalibFact' or 'FALSE' (no action);  final recalibration-factor used exported in result as $recalibFact
#' @param massTy (character) 'mono' or 'average'
#' @param prefFragPat (numeric) pattern for preferential fragmentation (see also Haverland 2017), if \code{NULL} default will be taken (in function \code{evalIsoFragm}) from \code{.prefFragPattern()} 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) additional messages and objects exportet to current session for debugging
#' @return This function returns a list with $massMatch (list of exerimental peptides matching to one or more predicted), $preMa (predicted ions, including fixed modif), $pepTab (predicted neutral peptides, wo modifications), $expMa (experimental mass from input), $recalibFact (recalibration factor as from input), $docTi (time for calculations)
#' @seealso \code{\link{makeFragments}}, \code{\link{identifVarModif}}, \code{\link{identifyPepFragments}}, \code{\link[wrMisc]{findCloseMatch}}
#' @examples
#' pro3 <- "HLVDEPQNLIK"
#' exp3 <- c( b4=465.2451, b5=594.2877, b6=691.3404,  y7=841.4772, y6=712.4347, y5=615.3819)
#' ident3 <- identifFixedModif(prot=pro3, expMass=exp3, minFragSize=4, 
#'   maxFragSize=60, modTy=list(basMod=c("b","y")))
#' ident3$massMatch                                                                                              
#' ## as human readable table:
#' ident3$preMa[ ident3$preMa[,"no"] %in% (names(ident3$massMatch)),]
#' @export
identifFixedModif <- function(prot, expMass, minFragSize=5, maxFragSize=60, indexStart=1, suplPepTab=NULL, internFra=TRUE, chargeCatchFilter=TRUE,
  maxMod=c(p=3,h=1,k=1,o=1,m=1,n=1,u=1,r=1,s=1), modTy=NULL, specModif=NULL, knownMods=NULL, identMeas="ppm",limitIdent=5,filtAmbiguous=FALSE, 
  recalibrate=FALSE, massTy="mono", prefFragPat=NULL, silent=FALSE, debug=FALSE, callFrom=NULL){         #
  ## identify predicted mass based on 'prot' (AA-sequence) compared to 'expMass' for given range of aa-length
  ## return list, ie result of massMatch() on 'pepTab' and 'expMass'
  ##    need to consider min intensity of experim values ??
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="identifFixedModif")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(is.null(maxFragSize) || !is.numeric(maxFragSize)) maxFragSize <- 400
  if(debug) {message(fxNa," .. xxidentifFixedModif0 ")}    #
  AAmass <- wrProteo::AAmass(massTy=massTy, inPept=TRUE)
  mH20 <- wrProteo::massDeFormula("2HO", massTy=massTy, callFrom=fxNa)
  massMatch1 <- recalibFact <- NULL   # initialize (for case no peptides remain after filtering charge-cathing AAs)
  docTi <- rep(NA,7)
  names(docTi) <- c("ini_identifFixedModif","makeFragments","countPotModifAAs","addMassModif","finUniqCheck","findCloseMatch","recalib")
  docTi[1] <- Sys.time()                                                            # 
  recalibFact <- 0
  recalibrate <- isTRUE(recalibrate)
  ##                     
  ## basic mass predictions : peptides (wo modification)
  massIni <- cbind(se=prot, na=names(prot),  mass=wrProteo::convAASeq2mass(prot, seqName=TRUE, callFrom=fxNa))
  ## make table of peptides (wo considering optional modifications)
  if(debug) {message(fxNa," .. xxidentifFixedModif1 ")}   #
  pepTab <- makeFragments(protTab=massIni, minFragSize=minFragSize,maxFragSize=maxFragSize, internFra=internFra,
    knownMods=knownMods, massTy=massTy, prefFragPat=prefFragPat, silent=silent, debug=debug, callFrom=fxNa)
  if(!identical(indexStart,1)) {
    if(!silent) message(fxNa," ** Increase pepTab index from ",min(as.integer(pepTab[,"no"]))," by ", indexStart," to ",min(as.integer(pepTab[,"no"])+indexStart-1))
    pepTab[,"no"] <- as.integer(pepTab[,"no"]) +indexStart -1}
  if(nrow(pepTab) <4) message(fxNa,"NOTE : only ", nrow(pepTab)," initial fragments predicted from makeFragments !!!")
  docTi[2] <- Sys.time()                                                            # makeFragments() consumes 95-99% of time !!  
  if(debug) {message(fxNa,"iFM1"); iFM1 <- list(prot=prot, expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1,chargeCatchFilter=chargeCatchFilter,chargeCatchFilter=chargeCatchFilter)}   #
  
  ## need to document number of peptides removed (in nRemPep) !!
  ## optional filter to remove all peptides wo charge-catching AAs
  nRemPep <- 0
  if(chargeCatchFilter) {   
    chaLi <- unique(unlist(sapply(.chargeCatchingAA()[,1], grep, pepTab[,"seq"])) )
    if(length(chaLi) >0) {
      if(length(chaLi) <nrow(pepTab)) {
        if(!silent) message(fxNa,"Removing ",nrow(pepTab) -length(chaLi)," out of ",nrow(pepTab)," (initial) peptides not containing any charge-catching residues")         
        pepTab <- pepTab[chaLi,] } else warning(fxNa,"Bizzare, when checking for non-charged peptides to remove found more than peptides")
    } else {
      message(fxNa,"PROBLEM : NO peptides remaining when filtering for peptides containing one of  ",nrow(.chargeCatchingAA())," charge-catching AAs !! (keep only 1st & 2nd)")
      if(debug) {message(fxNa,"iFM1b");  iFM1b <- list(pepTab=pepTab, chaLi=chaLi,recalibFact=recalibFact, recalibData=NULL, docTi=docTi,chargeCatchFilter=chargeCatchFilter)}
      pepTab <- if(nrow(pepTab) >2) pepTab[1:2,] else pepTab} }

  ## add custom/specific single location mass modifications (eg bound ions)
  if(debug) {message(fxNa,"Ready for single location mass modifications (bound ions)  iFM2"); iFM2 <- list(prot=prot, expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1,chargeCatchFilter=chargeCatchFilter)}
  if(length(specModif) >0) pepTab <- .singleSpecModif(pepTab, specModif, callFrom=fxNa, silent=silent, debug=debug)
  if(length(suplPepTab) >0) {if(ncol(pepTab)==ncol(suplPepTab)) pepTab <- rbind(suplPepTab, pepTab) else message(fxNa,"Problem with incompatible 'suplPepTab', ignoring !")}
  if(debug) {message(fxNa,"iFM3"); iFM3 <- list(prot=prot, expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1)}
  if(nrow(pepTab) >2) {
    if(length(modTy) >0) {
      ## fixed modifications ..  (later consider variable modif only when fixed modif found in experimtal data)
      modTb <- list(basMod=checkModTy(modTy)$basMod, varMod=NULL, varMo2=NULL)                  # copy for treating fixed modif only
      cou <- countPotModifAAs(pepTab=pepTab, modTy=modTb, maxMod=maxMod, silent=silent, debug=debug, callFrom=fxNa)
      if(debug) {message(fxNa,"iFM3b"); iFM3b <- list(prot=prot, expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1,cou=cou,modTb=modTb)}
      docTi[3] <- Sys.time()                                 #
      ## transform to ions
      preMa <- addMassModif(cou=cou$cou, pepTab=pepTab, combTerm=cou$combTerm, modTy=modTb, basVarMod="basMod", silent=silent, debug=debug, callFrom=fxNa)   # wo $varMod
      chColN <- colnames(preMa[[1]])=="mass"
      if(any(chColN)) colnames(preMa[[1]])[which(colnames(preMa[[1]])=="mass")] <- "finMass"
    } else {
      docTi[3] <- Sys.time() 
      message(fxNa,"NOTE : No modifications given, using neutral peptide masses")                                                           #
      preMa <- cbind(pepTab, finMass=pepTab[,"mass"], modSpec="", mod="")
    } 
    docTi[4] <- Sys.time()                                                            #
    if(debug) {message(fxNa,"iFM4"); iFM4 <- list(prot=prot,preMa=preMa, expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1,recalibFact=recalibFact)}
    if(is.list(preMa)) preMa <- preMa$pepTab
    docTi[5] <- Sys.time()                                                            #
    ## compare 'preMa' & 'expMass' : findCloseMatch(),.compareByPPM(),.compareByDiff()
    ## now filter experim values to go within range of predicted (no sense in testing even further ..)
    expMass <- as.numeric(expMass) + recalibFact
    names(expMass) <- 1:length(expMass)                           #
    preMaRa <- range(as.numeric(preMa[,"finMass"]), na.rm=TRUE) +c(-1,1)
    chExpM <- expMass > preMaRa[1] & expMass < preMaRa[2]
     if(debug) {message(fxNa,"iFM4a")}
    if(all(!chExpM)) { if(!silent) message(fxNa,"No hits found !! (range of experimental masses not compatible)")
      return(list(massMatch=list(), preMa=preMa, pepTab=pepTab, recalibFact=recalibFact, recalibData=NULL, docTi=docTi))
    } else {      
      if(!silent) message(fxNa,sum(chExpM)," out of ",length(chExpM)," experim masses in range of ",nrow(preMa)," predicted (max ",signif(preMaRa[2]-0.5,3),")")  
      if(any(!chExpM)) expMass <- expMass[which(chExpM)]
      ## 1st run of identification (and later address/identify fixed modifications) :
      ## last 'complete' table of fragments before 1st run of matching : xxFrag5$preMa
      ## massMatch1 is simple list with index-names of matches close enough & mass values
      predMa <- as.numeric(preMa[,"finMass"])
      names(predMa) <-  preMa[,"no"]
       if(debug) {message(fxNa," ..  iFM5"); iFM5 <- list(prot=prot,preMa=preMa,predMa=predMa, limitIdent=limitIdent,identMeas=identMeas,expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1)}
      ## also limit expMass to stay within range of predicted (and adjust names)
      massMatch1 <- wrMisc::findCloseMatch(x=predMa, y=expMass, compTy=identMeas, limit=limitIdent, sortMatch=FALSE, callFrom=fxNa, silent=silent) # 'sortMatch' as F for not inversung
      ## 'x' corrsp to preMa line, name of match to expMass
      if(!silent) message(fxNa," 1st pass: compare ",nrow(preMa)," predicted (incl ",sum(!is.na(preMa[,"ambig"]))," ambiguous : ",
        wrMisc::pasteC(unique(wrMisc::naOmit(preMa[,"ambig"]))),"\n     with  ",length(expMass)," input (measured) masses, found  " ,length(massMatch1),
        " groups of matches to experimental masses (total=",if(length(massMatch1) >0) sum(sapply(massMatch1,length)) else 0,"))")
      if(length(massMatch1) <1) message(fxNa,"\n ** NO matches found !! **\n    to masses like ",paste0(utils::head(preMa[,"finMass"]),collapse=" "),"...\n")
      if(debug) {message(fxNa," .. iFM6"); iFM6 <- list(prot=prot,preMa=preMa,predMa=predMa, expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1)}
      docTi[6] <- Sys.time()                                                            #
      ## problem : so far need also to export full pepTab for 2nd round search  -> gains in RAM diminish ... ==> need to add test for var modif
      ##
      ## RECALIBRATE
      dif <- NULL
      if(length(massMatch1) <21) {
        if(recalibrate) {recalibrate <- FALSE
          if(!silent) message(fxNa," ",length(massMatch1),"matches are insufficient for determining calibration factor") }
        recalibFact <- 0; dif <- NULL}
      if(recalibrate) {
        chLe1 <- sapply(massMatch1, length)==1
        if(sum(chLe1) > length(chLe1)/1.5 && sum(chLe1) >20) {                  # > 66.7% single hit and >20 pep
          msg <- c(fxNa,"RECALIBRATION :  based on diff of  ",sum(chLe1)," single hits :  ")
          if("diff" %in% identMeas) {
            dif <- unlist(massMatch1[which(chLe1)])                           # data finally used for recalibration
          } else {                                                            # prec for extracting diff in case of ppm
            expMaN <- massMatch1[which(chLe1)]
            names(expMaN) <- NULL
            expMaN <- names(unlist(expMaN))                                   # names of expMass to use
            dif <- as.numeric(preMa[match(names(massMatch1[which(chLe1)]), preMa[,"no"]), "finMass"]) -expMass[expMaN]
          }
        } else {
          chLe2 <- sapply(massMatch1,length) <4    # now also consider up to 3 matches
          if(sum(chLe2) < length(chLe2)/10) {recalibFact <- 0; dif <- 0; 
            if(!silent) message("\n++++",fxNa,"TROUBLE finding right data for recalibration ? (too few data below 4 matches)\n")
          } else {
          if(sum(chLe2) <15) { chLe2 <- sapply(massMatch1, length) <11
            if(!silent) message(fxNa,"Opening recalibration to all sets of with up to 10 hits")
          }
          msg <- c(fxNa,"RECALIBRATION :  based on median diff of  ",sum(chLe2)," muti-hits (below 4 hits):  ")
          if("diff" %in% identMeas) {
            dif <- sapply(massMatch1[which(chLe2)], function(x) {xz <- abs(x); x[which(xz==min(xz))]})
          } else { 
            match3 <- match3N <- massMatch1[which(chLe2)]
            names(match3N) <- NULL
            match3N <- names(unlist(match3N))        
            dif <- preMa[match(rep(names(match3), sapply(match3,length)), preMa[,"no"]), "finMass"] -expMass[match3N]}
        } }
        recalibFact <- wrMisc::stableMode(dif, histLike=TRUE, silent=silent, callFrom=fxNa)               #median(c(unlist(massMatch1[which(chLe1)]),zz))
        if(!silent) message(msg, signif(recalibFact,3))    }
      docTi[7] <- Sys.time()
      if(debug) {message(" ..  iFM7"); iFM7 <- list(prot=prot,preMa=preMa, expMass=expMass,pepTab=pepTab,massIni=massIni, massMatch1=massMatch1,recalibFact=recalibFact)}
      list(massMatch=massMatch1, preMa=preMa, pepTab=pepTab, recalibFact=recalibFact, recalibData=dif, docTi=docTi) }
    } else { if(!silent) message(fxNa,"TROUBLE ahead, pepTab is empty !")
      list(massMatch=massMatch1, preMa=NULL, pepTab=pepTab, recalibFact=0, recalibData=NULL, docTi=docTi)}
  }  


#' Add Single Specific Modifications
#'
#' Add single specific modification to peptide/protein fragments .
#'
#' @param pepTab (matrix) matrix of fragments (cols 'no','seq','orig','ty','seqNa','beg','end','precAA','tailAA','ambig','mass')
#' @param specModif (list) with elements 'modOrigin' (sequence), 'modPos' (position within sequence), 'modMass' (digits, ie mass to add),
#'   'modName' (name of modif), 'modFixed' (fixed or , logical) 
#' @param nMaxMod (numeric) max number a given modification may occur 
#' @param massTy (character) 'mono' or 'average'
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) additional messages and objects exportet to current session for debugging
#' @return This function returns a list with $massMatch (list of exerimental peptides matching to one or more predicted), $preMa (predicted ions, including fixed modif), $pepTab (predicted neutral peptides, wo modifications), $expMa (experimental mass from input), $recalibFact (recalibration factor as from input), $docTi (time for calculations)
#' @seealso \code{\link{makeFragments}}, \code{\link{identifVarModif}}, \code{\link{identifyPepFragments}}
#' @examples
#' pep1 <- c(pe1="KPEPTI")
#' # The table of possible terminal fragments (for simplicity terminal only)
#' pepTab1 <- makeFragments(pep1, min=3, max=7, internFra=FALSE)
#' specModif1 <- list(modOrigin=pep1, modPos=1, modMass=579.9663, modName="p", modFixed=FALSE)
#' .singleSpecModif(pepTab1, specModif1 )
#' 
#' protP <- c(protP="PEPTIDEKR")
#' pep1 <- c("PTI","KPE","EPTI")
#' papTab1 <- cbind(no=c(7,2,6),seq=pep1, orig=rep("KPEPTI",3), origNa=rep("pe1",3), 
#'   ty=paste0(c("C","N","C"),"ter"), beg=c(4,1,3), end=c(6,3,4),
#'   mass= wrProteo::convAASeq2mass(pep1, massTy="mono"), modSpec="")
#' 
#' @export
.singleSpecModif <- function(pepTab, specModif, nMaxMod=1, massTy="mono", callFrom=NULL, silent=FALSE, debug=FALSE) {  #massIni
  ## 'pepTab'  matrix of fragments (cols 'no','seq','orig','ty','seqNa','beg','end','precAA','tailAA','ambig','mass')
  ##    note : at this level pepTab is typically for neutral peptides, thus cumparsion to ions will seem 1 H+ too low (or 1 e- too high)
  ## 'specModif' .. list with elements 'modOrigin' (sequence), 'modPos' (position within sequence), 'modMass' (digits, ie mass to add),
  ##   'modName' (name of modif), 'modFixed' (fixed or , logical)
  ## 'nMaxMod' .. (numeric) max number a given modification may occur
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".singleSpecModif")
  if(debug) silent <- FALSE
  .extrSpcFeat <- function(x, spc, no=1) if(spc %in% names(x)) x[[which(names(x)==spc)]] else x[[no]]
  naChe <- modOrigin <- modPos <- modMass <- modName <- modFixed <- NULL
  if(debug) {message(fxNa," .. xxsingleSpecModif00 ")}  
  if(length(specModif) >0) { if(!is.list(specModif)) specModif <- as.list(specModif)
    if(length(specModif) <5) specModif <- NULL else {
      nEl <- sapply(specModif, length)
      if(any(nEl > 1) && any(nEl==1)) {for(i in which(nEl==1)) specModif[[i]] <- rep(specModif[[i]], max(nEl,na.rm=TRUE))
        if(!silent) message(fxNa," .. augmenting ",wrMisc::pasteC(names(specModif)[i])," to length ",max(nEl,na.rm=TRUE)) }
      if(debug) {message(fxNa," .. xxsingleSpecModif0 ")}
      modOrigin <- .extrSpcFeat(specModif, spc="modOrigin", no=1)
      modPos <- .extrSpcFeat(specModif, spc="modPos", no=2)
      modMass <- as.numeric(.extrSpcFeat(specModif, spc="modMass", no=3))
      modName <- .extrSpcFeat(specModif, spc="modName", no=4)
      modFixed <- as.logical(.extrSpcFeat(specModif, spc="modFixed", no=5))}
      msg <- " .. inconsistent length of args to 'specModif', ignoring 'specModif'"
      if(length(unique(sapply(specModif[1:5], length))) >1) {if(!silent) message(fxNa,msg); specModif <- NULL}
      nEl2 <- c(modOrigin=length(modOrigin), modPos=length(modPos), modMass=length(modMass), modName=length(modName), modFixed=length(modFixed))
      if(length(unique(nEl2)) >1) {
        unexpLe <- table(nEl2)
        unexpLe <- which(nEl2==as.numeric(names(unexpLe)[which(unexpLe==min(unexpLe,na.rm=TRUE))]))
        message(fxNa,"Problem : 'specModif' part ",wrMisc::pasteC(names(unexpLe))," of length ",unexpLe," won't fit to rest !!") }
      if(debug) {message(fxNa," .. xxsingleSpecModif1 ")}
      if(!is.numeric(modMass)) {if(!silent) message(fxNa," .. invalid 'modMass'"); specModif <- NULL}
      if(!is.logical(modFixed)) {if(!silent) message(fxNa," .. invalid 'modFixed'"); specModif <- NULL}
  }
  oblNames <- c("orig","beg","end","mass","no","seqNa","ty","modSpec")         # check 'pepTab' for these obligatory names
  chNames <- match(oblNames,colnames(pepTab))
  if(any(is.na(chNames))) stop(fxNa,"Can't find obligatory colnames ",wrMisc::pasteC(oblNames[is.na(chNames)],quoteC="'")," in input 'pepTab'")
  ## main
  if(length(modMass) >0) {    # which protein/peptide [modOrigin], AAposition [modPos], deltaMass [modMass], name [modName], fixedModif [modFixed], (short name)
    if(debug) {message(fxNa," .. xxsingleSpecModif1a")}    
    ## 1st step: try to match names of proteins 'modOrigin' to those from pepTab (ie to which proteins modifications apply)
    modOrigInd <- lapply(modOrigin, function(x) which(pepTab[,"origNa"]==x))       #now index to pepTab  #which(pepTab[,"origNa"]==modOrigin)
    naChe <- sapply(modOrigInd, function(x) (all(is.na(x)) | length(x) <1))      # find NA or empty list (no matches found)
    if(all(naChe)) {             ## can't find ANY of names given in modOrigin .. more elaborate search, possibly only half of name given
      pepTaNa <- unique(pepTab[,"origNa"])                             # want to find match these with modOrigin 
      pepTaNa <- matrix(unlist(sapply(pepTaNa, strsplit, "\\.")), ncol=2, byrow=TRUE)    # will cause problem if any of pepTab[,"origNa"] do NOT contain '\\.' !!
      chNa2 <- apply(pepTaNa, 1, function(x) sapply(modOrigin, function(y) y %in% x))                   # each col represents one of unique pepTab[,"origNa"]
      if(length(dim(chNa2)) <2) chNa2 <- matrix(chNa2, nrow=nrow(pepTaNa), dimnames=list(modOrigin,NULL))
      colnames(chNa2) <- paste0(pepTaNa[,1],pepTaNa[,2],sep=".")
      ## want to find modOrigin in pepTab[,"origNa"]
      ## oppose unique(rownames(chNa2)) {from specModif} to colnames(chNa2) [from pepTab)]
      for(i in 1:length(unique(modOrigin))) {
        ## check for full match
        j <- unique(modOrigin)[i]
        if(any(chNa2[j,])) {modOrigin[which(modOrigin==j)] <- rep(colnames(chNa2)[which(chNa2[j,])], sum(modOrigin==j))
        } else {        
          trimI <- wrMisc::.trimLeft(unique(modOrigin))[i]
          gr2 <- grep(trimI,colnames(chNa2))
          if(!silent) message(fxNa,"  matching '",j,"' by grep to ",wrMisc::pasteC(colnames(chNa2)[gr2], quoteC="'"))
          if(length(gr2) >0) modOrigin[which(rownames(chNa2) %in% j)] <- rep(colnames(chNa2)[gr2], sum(rownames(chNa2) %in% j))        
      } }      
      modOrigInd <- lapply(modOrigin, grep, pepTab[,"origNa"])                       # now index to pepTab 
      naChe <- sapply(modOrigInd, function(x) (all(is.na(x)) | length(x) <1))        # reset: find NA or empty list (no matches found)
      naChe <- if(all(sapply(modOrigInd, length)) >0) sapply(modOrigInd, function(x) (all(is.na(x)) | length(x) <1)) else TRUE      # reset: find NA or empty list (no matches found)
    }
    if(debug) {message(fxNa," .. xxsingleSpecModif1b "); xxsingleSpecModif1b <- list(pepTab=pepTab,specModif=specModif,nMaxMod=nMaxMod,massTy=massTy,naChe=naChe,
      naChe=naChe,modOrigin=modOrigin,modPos=modPos,modMass=modMass,modName=modName,modFixed=modFixed)}
    ##
    ## APPLY MODIFICATIONS to the correspoding fragments (of proteins concerned) : loop along valid (nonNA) modifications
    nIter <- 1                                                    # not really used any more !
    if(any(!naChe)) {
      tmp <- sapply(modName, rep, nMaxMod)                          # prepare repeated abbreviation for spec modif for search by grep (specModif$modName)
      tmp <- if(length(dim(tmp)) <2) as.matrix(tmp, nrow=nMaxMod) else t(tmp)
      maxModNa <- apply(tmp, 1, paste, collapse="\\.")                          # max concatenated repeat-names of each modif
      for(i in which(!naChe)) {                                   # run loop along modifications ..
       if(debug) message(fxNa,"checking for modif ",) 
        modPep <- which(as.integer(pepTab[,"beg"]) <= modPos[i] & as.integer(pepTab[,"end"]) >= modPos[i] & pepTab[,"origNa"]==modOrigin[i])  #index rel to pepTab
        names(modPep) <- pepTab[modPep,"no"]
        if(debug && i==which(!naChe)[1]) {message(fxNa," .. xxsingleSpecModif2a")}
        if(debug && i==which(!naChe)[2]) {message(fxNa," .. xxsingleSpecModif2b")}
        if(length(modPep) >0) {                   # ie current modif indeed occurs somewhere.. 
          chMaxMod <- grep(paste0("\\.", maxModNa[i],"\\.",sep=""), pepTab[modPep,"seqNa"])            # check if current modif already performed (eg prev loop) 
          if(length(chMaxMod) > 1) {              # modification exists, examine max no
            if(nMaxMod ==1) modPep <- modPep[-1*chMaxMod] else {
              ## need to count no modifs already present
              nRep <- (nchar(pepTab[modPep,"seqNa"]) -nchar(gsub(paste0("\\.", maxModNa[i],"\\.",sep=""),"", pepTab[modPep,"seqNa"])))/(2+nchar(maxModNa[i]))
              chRep <- nRep > nMaxMod             # not yet at max no
              modPep <- if(any(chRep)) modPep[which(chRep)] else NULL }}}
        if(length(modPep) >0) {                   # valid cases of current modif         
          ## now duplicate concerned lines of pepTab to new lines at end - if variable modif (ie keep orig at end)
          nPepT0 <- nrow(pepTab)                                                                  # initail no of peptides (if only fixed modif)
          if(!modFixed[i]) { 
            pepTab <- rbind(pepTab, pepTab[modPep,])                               # varModif: add new lines for keeping orig peptides            
            if(!silent) message(fxNa,modName[i],"/",modOrigin[i],"/",modPos[i]," variable modif: augmenting by ",length(modPep)," to ",nrow(pepTab)," peptides")
            pepTab[nPepT0 +(1:length(modPep)),"no"] <- max(nPepT0, as.integer(pepTab[,"no"]), na.rm=TRUE) +(1:length(modPep))   # increase no
          } else if(!silent) message(fxNa,modName[i],"/",modOrigin[i],"/",modPos[i]," fixed modif ",length(modPep)," peptides OK within range for modif")   
          chNo <- duplicated(pepTab[,"no"])
          if(any(chNo) && !silent) message(fxNa," +++ DUPLICATED 'no' in i=",i)
          pepTab[modPep,"mass"] <- as.numeric(pepTab[modPep,"mass"]) + modMass[i]                 # add the current modification mass
          ## new col in pepTab initallly filled with ""#
          pepTab[modPep,"modSpec"] <- paste0(sub("NA","", pepTab[modPep,"modSpec"]), modName[i])    # add modif name 
          pepTab[modPep,"seqNa"] <- paste0(gsub("[[:digit:]]+-[[:digit:]]+$","", pepTab[modPep,"seqNa"]),    # add modif name to 'seqNa' with '.'
            modName[i],".",pepTab[modPep,"beg"],"-", pepTab[modPep,"end"],sep="")
          nIter <- nIter +1 } 
       } }
    pepTab <- pepTab[order(pepTab[,"orig"], as.integer(pepTab[,"beg"]), as.integer(pepTab[,"no"]), as.integer(pepTab[,"end"])),]               # re-order (1st by orig then by 'beg' & 'end')
    chNo <- duplicated(pepTab[,"no"])
    if(any(chNo)) message(fxNa,"Trouble ahead: Problem with non-unique fragment numbers") } 
  pepTab }
   
