#' Add modifications to peptide mass
#'
#' Adjust/add mass for modifications from 'modTy' to all peptides in 'pepTab' based on count 'cou' of occurances of modifications :
#' Either fixed or variable modifications will be added to the mass of initial peptides from argument \code{papTab}.
#' Terminal ionization (like 'b' or 'y' -fragments) is treated as fixed modification and the resulting masses will correspond to standard mono-protonated ions.
#' Since variable and fixed modification types can't be run in a single instance, the function has to get calles twice, it is recommended to always start with the fixed modfications,
#' In the case of fixed modifications (like defining 'b' or 'y' fragments) neutral peptide masses should be given to add the corresponding mass-shift (and to obtain mono-protonated ions).
#' In case of variable modifications (like 'd' or 'p'), the corresponding ions from the fixed modifications should get furnished to add the corresponding mass-shift,
#' the masses resulting from the initial fixed modifications run can be used.  
#' Note, that transforming a neutral precursor M into MH+ is also considered a modification. 
#' The results are also correct with obligatory fragments that can't occur the same time (eg x & y ions can't be same time, need to make add'l lines...).
#' This function has a multiprocessor mode, with small data-sets (like the toy example below) there is typcally no gain in performance. 
#' @param cou (list) list of matrixes with counts for number of modifications per peptide
#' @param pepTab (matrix) table with peptide properties
#' @param combTerm (matrix) table with separate rows for $basMod that are exclusive (ie can't be accumulated, eg x & y ions)
#' @param modTy (character) list of modification types to be considered
#' @param lastIndex (integer) index-1 (ie last index from prev matrix) from which new peptide-variants should start from 
#' @param modChem (character) optional modifications
#' @param basVarMod (character) toggle if fixed ('basMod') or variable ('varMod') modificatons should be calculated
#' @param massTy (character) default 'mono' 
#' @param knownMods (list) optional custom definition whoch modification is N-term, etc (see \code{\link{AAfragSettings}} 
#' @param nProc (integer) number of processors in case of multi-processor use (requires Bioconductor package \code{BiocParallel}) 
#' @param parallDefault (logical) for use of other/previously set \code{register(bpstart())} in case \code{.parCombinateAllAndSum} is called
#' @param silent (logical) suppress messages
#' @param debug (logical) for bug-tracking: more/enhanced messages and intermediate objects written in global name-space 
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return list of $pepTab (table of peptide as single charge positive ions), $abc ('representative' list of all combinations to add). Main result in $pepTab
#' @seealso \code{\link[wrMisc]{convToNum}}
#' @examples
#' pep1 <- c(pe1="KPEPTI")
#' # The table of possible terminal fragments (for simplicity terminal only)
#' pepTab1 <- makeFragments(pep1, min=3, max=7, internFra=FALSE)
#' # Which fragment may be subject to how many modification (including ionization by H+)
#' cou1 <- countPotModifAAs(pepTab=pepTab1, modTy=list(basMod=c("b","y")))
#' # Add modifications (here: ionize all pepptides by H+)
#' preMa1 <- addMassModif(cou=cou1$cou, pepTab=pepTab1, combTerm=cou1$combTerm, 
#'   modTy=list(basMod=c("b","y")), basVarMod="basMod")
#' preMa1
#' 
#' ## Example including variable modifications
#' modT3 <- list(basMod=c("b","y"),varMod=c("p","h","d"))
#' cou3 <- countPotModifAAs(pepTab=pepTab1, modTy=modT3)
#' ## Now we re-use/inject the results for the fixed modificatons
#' preMa3 <- addMassModif(cou=cou3$cou, pepTab=preMa1$pepTab, combTerm=cou1$combTerm, 
#'   modTy=modT3, basVarMod="varMod")
#' head(preMa3$pepTab,12)
#' @export
addMassModif <- function(cou, pepTab, combTerm, modTy, lastIndex=NULL, modChem=NULL, basVarMod="basMod", massTy="mono", knownMods=NULL, nProc=1, parallDefault=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL){
  ## adjust/add mass for modifications from 'modTy' to all peptides in 'pepTab' based on count 'cou' :
  ## also OK woth obligatory fragments that can't occur the same time (eg x & y ions can't be same time, need to make add'l lines...)
  ##  return list of $pepTab,$basMod,$varMod,$abc ('representative' list of all combinations to add)
  ##  main result in $varMod (inlcudes basMod, has add'l col 'isoMasMod' for warning) if empty (no var modifs done) use $basMod
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="addMassModif")
  if(is.null(modChem)) modChem <- AAfragSettings(outTy="all")$modChem
  msg <- " 'basVarMod' should be either 'basMod' or 'varMod' (length=1) "
  if(length(basVarMod) >1) {message(fxNa,"truncating"); basVarMod <- basVarMod[1]}
  if(!any(c("basMod","varMod") %in% basVarMod)) stop(msg)
  modTy <- checkModTy(modTy,knownMods=knownMods,silent=silent,callFrom=fxNa)
  abc <- NULL
  dataOK <- FALSE  
  if(length(cou) <1 | length(pepTab) <1) {
    if(!silent) message(fxNa," 'cou' and/or 'pepTab' is/are empty - nothing to do !")
  } else if(nrow(pepTab) >0) {
    dataOK <- TRUE
    if(length(lastIndex) <1) lastIndex <- if("no" %in% colnames(pepTab)) max(pepTab[,"no"], na.rm=TRUE) else nrow(pepTab)
  }
  if(dataOK) {
    ## finish defining mass-modifications :
    thisIsBasMod <- any(nchar(modTy$basMod) >0) & basVarMod %in% c("basMod","all")
    uniqCo1 <- wrMisc::firstOfRepLines(cou[[if(thisIsBasMod) "basMod" else "varMod"]], outTy="all")                      # unique combination schemes (for easy repeating)
    if(debug) {message(" .. xxaddMassModif0 \n")}  
    couX <- cou[[if(thisIsBasMod) "basMod" else "varMod"]][uniqCo1$ind,]                                                 # changed 26oct17
    if(!is.matrix(couX)) couX <- matrix(couX, nrow=length(uniqCo1$ind), dimnames=list(
      rownames(cou[[if(thisIsBasMod) "basMod" else "varMod"]])[uniqCo1$ind], colnames(cou[[if(thisIsBasMod) "basMod" else "varMod"]]) ))
    if(any(couX >0)) {
      ## remove cols/modifs not encountered
      chCol <- !(colSums(couX) <1 & colnames(couX) %in% names(AAfragSettings("specAAMod"))) 
      if(all(!chCol)) {
        couX <- 0
        abc <- NULL
        if(!silent) message(fxNa," nothing to do for mass-modifications")
      } else if(any(!chCol)) {couX <- if(sum(chCol) >1) couX[,which(chCol)] else {
         matrix(couX[,which(chCol)], ncol=1, dimnames=list(rownames(couX), colnames(couX[which(chCol)]))) }}}
    if(debug) {message(" .. xxaddMassModif0b\n")}
    ## complete/finish for non-AAspec ie terminal modifications
    if(any(combTerm >1)) {                              # continue (add more lines) for obligatory modifs that can't occur the same time (eg y & z ions)
      ii <- if(nrow(combTerm) >2) 2 else 1
      for(i in ii:nrow(combTerm)) {         
        useCol <- which(!combTerm[i,] %in% combTerm[(i-1):1,])
        if(length(useCol) >0){
          useLi <- sort(unique(unlist(lapply(combTerm[-i,useCol], grep, add2FraNames))))
          massMod2 <- wrProteo::massDeFormula(modChem[match(combTerm[i,], modChem[,1]),2], massTy=massTy, silent=TRUE, callFrom=fxNa)   # mass modifications (simple)
          couX <- cou$basMod[useLi,]
          colnames(couX) <- massMod2
          massMod2 <- rowSums(matrix(as.numeric(wrMisc::conv01toColNa(couX, pasteCol=FALSE)), nrow=nrow(couX)), na.rm=TRUE)        # mass modifications in order of output
          colnames(couX) <- combTerm[i,]
          tmX <- cbind(no=useLi, modif=wrMisc::.pasteCols(wrMisc::conv01toColNa(couX)), mass=as.numeric(pepTab[useLi,"mass"]) +massMod2)
          basMod <- rbind(basMod,tmX) }
      } }
    ## note : dephospho q  and loss of water will give same mass !  
    ## prepare mass changes, split as fixed OR var modif
    if(any(couX >0)) {
      if(any(nchar(modTy$basMod) >0) & basVarMod %in% c("basMod","all")) {
        ## this is FIXED modif !
        if(debug) {message(" .. xxaddMassModif1a - fixed modif \n")}
        if(!"mod" %in% colnames(pepTab)) pepTab <- cbind(pepTab, mod=rep("", nrow(pepTab)))
        newNa <- wrProteo::massDeFormula(modChem[match(if(is.matrix(couX)) colnames(couX) else names(couX),modChem[,1]),2], massTy=massTy, silent=TRUE, callFrom=fxNa)
        chNewN <- newNa==0
        if(any(chNewN)) names(newNa)[which(chNewN)] <- ""
        if(is.matrix(couX)) colnames(couX) <- newNa else couX <- matrix(couX, nrow=1, dimnames=list(NULL,newNa))    
        mod <- .multMatByColNa(couX)
        mod <- mod[uniqCo1$num] 
        add2FraNames <- wrMisc::.pasteCols(wrMisc::conv01toColNa(cou[[1]]))                         # only oblig modifs
        protMa <- wrProteo::.atomicMasses()["H",massTy]              # need to ionize in pos mode ...
        pepTab[,"mass"] <- as.numeric(pepTab[,"mass"]) +mod +protMa
        pepTab[,"mod"] <- paste0(pepTab[,"mod"],add2FraNames,sep="")
        ## need to re-check for iso-masses
        chAmb <- pepTab[,"ambig"]=="isoMass"
        if(any(wrMisc::naOmit(chAmb))) pepTab[which(chAmb),"ambig"] <- NA
        chMa <- duplicated(as.numeric(pepTab[,"mass"]), fromLast=FALSE)
        if(any(chMa)) {
          chM2 <- duplicated(as.numeric(pepTab[,"mass"]), fromLast=TRUE)
          pepTab[which(chMa | chM2),"ambig"] <- "isoMass" }
        ## better reconstruct full name
        supNa <- gsub(" ",".",pepTab[,"modSpec"])
        chHeadPo <- nchar(supNa) >0 & substr(supNa,1,1) != "."
        if(any(chHeadPo)) supNa[which(chHeadPo)] <- paste0(".",supNa,sep="")       # add heading '.' separartor if any special modif      
        rownames(pepTab) <- paste0(pepTab[,"origNa"],".",pepTab[,"beg"],"-",pepTab[,"end"],supNa,".",pepTab[,"mod"],sep="")
      } else {    
        ## this is VARIABLE modif !
        if(debug) {message(" .. xxaddMassModif1b - variable modif \n")}
        chMod <- sum(cou$varMod, na.rm=TRUE)
        if(chMod <1) {
          ## no sites found, nothing to do ..
          abc <- NULL
        } else {
        ## cou$varMo2 remove q (since p+q =0 modif); ==> finally use $varMo2      
        if(length(cou$varMo2) < length(cou$varMod)) {
          cou$varMo2 <- cou$varMod
          chMo <- colnames(cou$varMod) %in% "q"
          if(all(c("p","q") %in% colnames(cou$varMod))) cou$varMo2 <- if(sum(!chMo) >1) cou$varMo2[,which(-chMo)] else {
            matrix(cou$varMo2[,which(-chMo)], ncol=1, dimnames=list(rownames(cou$varMod), colnames(cou$varMod)[which(-chMo)]))}}
        uniqCo2 <- wrMisc::firstOfRepLines(cou$varMo2, outTy="all", callFrom=fxNa)              # unique combination schemes (for easy repeating)
        tm2 <- match(colnames(cou$varMo2), modChem[,1])
        massModV <- wrProteo::massDeFormula(modChem[tm2,2], massTy=massTy, silent=TRUE, callFrom=fxNa)       # variable mass modifications
        names(massModV) <- modChem[tm2,1]
        nVMod <- unique(sort(cou$varMo2))                                     # number of types of indiv modifications (across cou$varMod)
        nVMod <- nVMod[nVMod >0]
        uniqCo <- cou$varMo2[uniqCo2$ind,]                                    # max no of modifications by type/group (lines)
        ## values not exceeding max no of modifs normally already considered in countPotModifAAs() making cou
        if(length(dim(uniqCo)) <2)  uniqCo <- matrix(uniqCo, nrow=length(uniqCo2$ind), ncol=ncol(cou$varMo2), dimnames=list(rownames(cou$varMo2)[uniqCo2$ind],colnames(cou$varMo2)))
        chPa <- requireNamespace("BiocParallel", quietly=TRUE)
        isWin <- "windows" %in% .Platform$OS.type
        if(!chPa) { message(fxNa,": package 'BiocParallel' not installed, can't run parallel processing")
          nProc <- 1}
        ## main
        if(debug){
          msg <- if(nrow(uniqCo)==1) "too few peptides for multi-proc" else {  if(nProc <2) "not multi-proc" else "multi-proc"}
          message(fxNa," - ",msg) }        
        ## may take more time when nProc >1 
        abc <- if(nrow(uniqCo)==1) combinateAllAndSum(as.numeric(uniqCo), massModV, notSingle=c("q","p"), callFrom=fxNa, silent=silent) else {        # (representative) list of all combinations to add
          if(nProc <2) { apply(uniqCo,1, combinateAllAndSum, massModV, notSingle=c("q","p"), callFrom=fxNa, silent=silent)      
        } else .parCombinateAllAndSum(uniqCo, massModV, nProc=nProc, parRegDefault=parallDefault, callFrom=fxNa)}
        if(!is.list(abc)) {abc <- list(abc); names(abc) <- as.character(rownames(uniqCo))} 
        if(debug) {message(" .. xxaddMassModif3\n")}
        ## now add/dispatch abc to each peptide concerned  (uses abc !)
        ab0 <- lapply(abc, function(x) {ch0 <- x %in% 0; if(all(ch0)) NULL else {if(any(ch0)) x[which(!ch0)] else x}})   # remove ocuurances of 0 mass shift
        couY <- cou$varMo2[which(rowSums(cou$varMo2) >0),]
        if(length(dim(couY)) <2) couY <- matrix(couY, ncol=ncol(cou$varMo2), dimnames=list(rownames(cou$varMo2)[which(rowSums(cou$varMo2) >0)], colnames(cou$varMo2)))
        uniqCo3 <- wrMisc::firstOfRepLines(couY, outTy="all", callFrom=fxNa)                                # unique combination schemes (for easy repeating)
        chLe <- sapply(ab0,length) <1
        if(all(chLe)) { message(fxNa," no modifications left !!"); return(list(pepTab=pepTab, abc=abc))
        } else { if(any(chLe)) ab0 <- ab0[which(!chLe)] }
        if(debug) {message(" .. xxaddMassModif3a\n") }
        ## clean list of types of mass-changes (ab0) for repeating mass-changes , sort by alphabet names to get 'd' displayed instead of 'h'
        ab0 <- lapply(ab0, function(x) {if(length(x) >1) {x <- x[order(names(x))]; x[!duplicated(x, fromLast=FALSE)]} else x })
        ## need to introduce mass-change of modifs & names of modifs to subset of main table
        names(ab0) <- NULL                                                   # will get composed modif-names otherwise
        modNa <- names(unlist(ab0))
        if(length(modNa) <1) modNa <- unlist(sapply(abc,function(x) rownames(x)[which(x !=0)])) else {
          if(any(nchar(names(modNa)) <1)) modNa <- unlist(sapply(abc, function(x) rownames(x)[which(x !=0)]))}
        nRep <- sapply(ab0,length)[uniqCo3$num]
        repI <- rep(which(rowSums(cou$varMod) >0), nRep)
        pepTa3 <- if(length(repI)==1) matrix(pepTab[repI,], nrow=1, dimnames=list(rownames(pepTab)[repI],colnames(pepTab))) else pepTab[repI,]
        addToNa <- unlist(sapply(ab0[uniqCo3$num],names))
        if(debug) {message(" .. xxaddMassModif3b\n") }
        names(ab0) <- NULL                                                   # will get composed modif-names otherwise
        ab0 <- unlist(ab0[uniqCo3$num])
        pepTa3[,"mass"] <- as.numeric(pepTa3[,"mass"]) + ab0
        modColNo <- wrMisc::naOmit(match(c("mod","modif","modSpec"), colnames(pepTa3)))[1]   # search for comumn to use for adding modif-names
        if(any(nchar(names(ab0)) <1)) message(fxNa," Trouble ahead !?  Some variable modifications names don't appear !") 
        pepTa3[,modColNo] <- paste0(pepTa3[,modColNo], addToNa, sep="")      # add var mod name to modif column 
        addSpe <- gsub(" ", ".", pepTa3[,"modSpec"])
        chSpe <- nchar(addSpe) >0 & substr(addSpe,1,1) != "."
        if(any(chSpe)) addSpe[which(chSpe)] <- paste0(".",addSpe[which(chSpe)],sep="")   # obtain heading '.' if followed by something
        chModSpe <- grep("modSpe",colnames(pepTa3))
        if(length(chModSpe) >0) colnames(pepTa3)[chModSpe[1]] <- "mod"       # reset pepTa3 to basic colnames
        rownames(pepTa3) <- paste0(pepTa3[,"origNa"],".",pepTa3[,"beg"],"-",pepTa3[,"end"], addSpe,sep="")   # add var mod name to rownames
        if(debug) {message(" .. xxaddMassModif3e\n") }
        ## make unique (new) index for var modif
        pepTa3[,"no"] <- as.integer(lastIndex) +1 + 1:nrow(pepTa3)                          # increase index
        chIso <- pepTa3[,"ambig"] %in% "isoMass"
        if(any(chIso)) pepTa3[which(chIso),"ambig"] <- NA                    # need to re-check iso-mass once combined ...
        if(debug) {message(" .. xxaddMassModif4\n") }
        ## final fusing identif fixed modif and var modif (& re-check for 'isoMass') 
        rownames(pepTab) <- paste0(pepTab[,"seqNa"], ".", pepTab[,"modSpec"], sep="")
        pepTab <- rbind(pepTab, pepTa3)
        chInd <- duplicated(pepTab[,"no"])
        if(any(chInd)) message(fxNa," BEWARE ! Some index numbers not unique !!")
        chIso <- duplicated(pepTab[,"mass"], fromLast=FALSE)
        if(any(chIso)) {
          chIs2 <- duplicated(pepTab[,"mass"], fromLast=TRUE)
          pepTab[which(chIso | chIs2),"ambig"] <- "isoMass" }
        ## why does this add a duplicted column named "seqNa" -> remove
        if(colnames(pepTab)[ncol(pepTab)]=="seqNa") pepTab <- pepTab[,-ncol(pepTab)]   # varModif : remove redundant column 'seqNa'
    } } } }
  list(pepTab=pepTab, abc=abc)}       

#' @export
.multMatByColNa <- function(mat,sumByRow=TRUE,...) {             
  ## multiply values of 'mat' by its colnames (numeric equivalent to conv01toColNa() which repates concatenated text)
  out <- matrix(as.numeric(mat)*as.numeric(rep(colnames(mat), each=nrow(mat))), nrow=nrow(mat))
  if(sumByRow) {out <- rowSums(out); names(out) <- rownames(mat)} else names(out) <- rownames(mat)
  out }
   
