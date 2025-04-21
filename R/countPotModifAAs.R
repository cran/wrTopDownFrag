#' Make Table With Counts of Potential Modification Sites
#'
#' Makes table 'cou' with counts of (potential) modification sites based on column 'seq' in matrix 'pepTab'.
#' Note: if multiple N-or C-term modifs, then only the first is shown in resulting table 'cou'.
#' 
#' @param pepTab (matrix) peptide sequences, start and end sites, typically result from \code{\link{makeFragments}} 
#' @param modTy (list) modifications : $basMod for character vector of fixed modifications and $varMod for variable modifications. For one letter-code see AAfragSettings("modChem")
#' @param maxMod (integer) maximal number variable modifications will be considered in given fragment (may increase complexity and RAM consumption)
#' @param specAAMod (list) optional custom list showing which AA to be considered with which (one-letter) modification code (default \code{\link{AAfragSettings}})
#' @param knownMods (list) optional custom list showing which modification appears at what type of location, eg N-terminal, internal ... (default \code{\link{AAfragSettings}})
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) for bug-tracking: more/enhanced messages and intermediate objects written in global name-space  
#' @return list of matrixes $cou and $combTerm, with number of modifications per peptides (line in 'pepTab') for basMod, varMod & varMo2
#' @seealso \code{\link{AAfragSettings}}, \code{\link{makeFragments}} 
#' @examples
#' protP2 <- c(mesp="MESPEPTIDES", pepe="PEPEPEP")
#' pepTab1 <- makeFragments(protTab=protP2, minFra=6, internFr=TRUE, massTy="mono")
#' cou1 <- countPotModifAAs(pepTab=pepTab1, modTy=list(basMod=c("b","y"),
#'   varMod=c("p","h")))
#' modTy2 <- list(basMod=c("b","y","h"), varMod=c("x","p","o","q","e","j"))
#' cou2 <- countPotModifAAs(pepTab=pepTab1, modTy=modTy2)
#' @export
countPotModifAAs <- function(pepTab, modTy, maxMod=c(p=3,h=1,k=1,o=1,m=1,n=1,u=1,r=1,s=1), specAAMod=NULL, knownMods=NULL, silent=FALSE, callFrom=NULL, debug=FALSE){
  ## make table 'cou' with count of modifications based on column 'seq' in matrix 'pepTab'
  ## return list of matrixes with number of modifications per peptide (line in 'pepTab') for  basMod, varMod & varMo2
  ## note: if multiple N-or C-term modifs, then only the first is shown in resulting table 'cou'
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="countPotModifAAs")
  if(nrow(pepTab) <1) return(list(cou=NULL, combTerm=NULL)) else {
  restrMod <- c("basMod","varMod")
  if(is.null(specAAMod)) specAAMod <- AAfragSettings(outTy="all")$specAAMod
  if(is.null(knownMods)) knownMods <- AAfragSettings(outTy="all")$knownMods
  if(is.null(maxMod)) maxMod <- c(p=3,h=1,k=1,o=1,m=1,n=1,u=1,r=1,s=1)
  ## table for converting names of fragment types:
  useKnoMo <- cbind(c("Nterm","Cterm","NCterm","intern","any"), c("Nter","Cter","NCter","inter","any")) 
  modTy <- checkModTy(modTy, knownMods=useKnoMo[1:2,1], silent=silent, callFrom=fxNa)
  useModT <- which(names(modTy) %in% restrMod & sapply(modTy, function(x) {if(length(x) >0) any(nchar(x) >0) else FALSE}))
  if(debug) {message(fxNa," .. xxcountPotModifAAs0")}
  pepTSup <- cbind(protIndex=as.integer(as.factor(pepTab[,"origNa"])), isTerm=pepTab[,"ty"] %in% c("Nter","Cter","full")) 
  ## multiple obligat modif (basMod) may be exclusive (can't be on same fragm) : eg multiple Cterm modifications :
  ##  select single pair (and add results adjusted by additive factor at end)
  useGr <- as.list(useKnoMo[,2])
  names(useGr) <- useKnoMo[,1]
  termMod <- lapply(modTy[useModT], function(x) sapply(knownMods[useKnoMo[,1]], function(y) y %in% x))                  # any terminal or internal
  chMultNC <- sapply(lapply(modTy[useModT], function(x) sapply(knownMods[useKnoMo[1:4,1]], function(y) y %in% x)),sapply,sum)    #which known (exclusive) N or C-term modifs present in modTy
  if(any(unlist(termMod[[1]]))) {                                                     # terminal modifs exist .. 
    combTerm <- sapply(termMod[[1]], sum, na.rm=TRUE)
    combTerm <- matrix(0+combTerm, nrow=1, dimnames=list(NULL,names(combTerm)))         # matrix, 2nd and later lines: which modification(pairs) need to be done later !!
  } else combTerm <- matrix(nrow=0, ncol=2, dimnames=list(NULL,c("Nterm","Cterm")))
  ## reduce/adjust init testing so that no mutually exclusive modifs remain present
  if(any(chMultNC[c("Nterm","Cterm"),] >1)) {         
    modTy <- lapply(modTy, function(x)  x[c(which(x %in% knownMods[["Nterm"]])[1],which(x %in% knownMods[["Cterm"]])[1], 
       which(x %in% unlist(knownMods[-1*match(c("Nterm","Cterm","intern"),names(knownMods))] ))) ] )
    if(!silent) message(fxNa,"Avoid exclusive modifications : adjusting modTy$basMod modifications to ",modTy$basMod,"  ")
  }
  if(debug) {message(fxNa," .. xxcountPotModifAAs1"); xxcountPotModifAAs1 <- list(pepTab=pepTab,modTy=modTy,maxMod=maxMod,useModT=useModT,specAAMod=specAAMod,knownMods=knownMods)}     
  ## make table 'cou' with count of modifications :  count no of AA for dependent modifs for basMod/varMo2 (varMo2 wo dephospho -> use to create varMod later)
  ## NEW CHANGES 2oct19: add col with protein-index !
  ## 1st step : 'protIndex', 'isTerm' terminal info in pepTSup
  ## 2nd : make non-redund (?, need index of orig ?)
  ## how to integrate shared between mult prot ??? (accompagnig list -of same length-with prot indexes ?  
  ## make pep seq unique within prot ? (if no terminal/internal mixing ?)
  ##  .. before counting aa spec events per/peptide    
  cou <- lapply(modTy[useModT], function(x) if(length(x) >0) .countModif(pepTab[,"seq"], modTyp=x, specAAMod, knownMods=knownMods, silent=debug, debug=debug))           # count $basMod (ie wo $varMod)
  if(debug)  {message(fxNa," .. xxcountPotModifAAs2")}        #,chMod=chMod
  ## consider max number of optional modifications : (eg max phospho & max de-pho )
  if(length(maxMod) >0) {
    for(i in names(cou)) {
      ## complete cou : search and modify which parts contain terminal modif (not yet integrated to cou)
      chTerm <- colnames(cou[[i]]) %in%  c(knownMods$NCterm)
      ## mark internal
      chInt <- pepTab[,"ty"]=="inter"
      chMoTy <- colnames(cou[[i]]) %in%  c(knownMods$intern)
      if(any(chInt) & any(chMoTy)) {
        cou[[i]][which(chInt),which(chMoTy)] <- 1 } 
      ## mark various variants of terminal
      if(any(chTerm)) for(k in which(chTerm)) {
        ch2 <- pepTSup[,"isTerm"] >0       
        if(any(ch2)) cou[[i]][which(ch2),k] <- 1  }
      chTerm <- colnames(cou[[i]]) %in%  c(knownMods$Nterm)
      if(any(chTerm)) for(k in which(chTerm)) {
        ch2 <- pepTab[,"ty"] =="Nter"       
        if(any(ch2)) cou[[i]][which(ch2),k] <- 1  }
      chTerm <- colnames(cou[[i]]) %in%  c(knownMods$Cterm)
      if(any(chTerm)) for(k in which(chTerm)) {
        ch2 <- pepTab[,"ty"] =="Cter"       
        if(any(ch2)) cou[[i]][which(ch2),k] <- 1  }
      chTerm <- colnames(cou[[i]]) %in%  c(knownMods$spcNterm)
      if(any(chTerm)) for(k in which(chTerm)) {           # spcNterm : set non-terminal to 0
        ch2 <- pepTab[,"ty"] !="Nter"    
        if(any(ch2)) cou[[i]][which(ch2),k] <- 0  }
      chTerm <- colnames(cou[[i]]) %in%  c(knownMods$spcCterm)
      if(any(chTerm)) for(k in which(chTerm)) {           # spcCterm; set non-terminal to 0
        ch2 <- pepTab[,"ty"] !="Cter"    
        if(any(ch2)) cou[[i]][which(ch2),k] <- 0  }
      ## complete cou : correct maxMod
      chMaxM <- colnames(cou[[i]]) %in% names(maxMod)                     # see if any defind type of modif with max number (of modif) to consider present
      if(any(chMaxM)) for(k in which(chMaxM)) {                           # loop along modifs
        chLi <- cou[[i]][,k] > maxMod[colnames(cou[[i]])[k]]
        if(any(chLi)) cou[[i]][which(chLi),k] <- maxMod[colnames(cou[[i]])[k]]
        }}}          
  if(debug)  {message(fxNa," .. xxcountPotModifAAs3")}        #
  if(any(sapply(termMod[[1]][c("Nterm","Cterm","NCterm")], sum) >1)) message("  ",fxNa," NOTE : MULTIPLE terminal modifications for ","(finish by new fx)")  # make function to give row&colnames of elements >thrsh
  ## so far the speaciall AA-linked modifs are counted
  ## complete cou : search and modify which parts contain terminal modif (not yet integrated to cou) ; limit to fixed modif ?
  if(debug) {message(fxNa," .. xxcountPotModifAAs3b")}       
  ## mutually excluding modifs
  whMod <- sapply(modTy,function(x) any(nchar(x))) >0
  if(any(whMod)) whMod <- names(modTy)[which(whMod)[1]] else message(fxNa,"don't understad which group of modifications to treat !?!")
  ## make 'varMo2' (with de-phspho) if given in 'modTy$varMod'
  if(all(c("p","q") %in% modTy$varMod,"p" %in% colnames(cou$varMod))) {       # for not running de-phospho alone ! (use cou$varMo2 ONLY for combinations)
    cou$varMo2 <- cbind(cou$varMod, q=cou$varMod[,"p"])}     
  ## making $varMo2 (including potential de-phospho) is useful but increases memory charge ! avoid !?!
  if(debug) {message(fxNa," .. xxcountPotModifAAs4")} 
  list(cou=cou,combTerm=combTerm) }}
  
  
#' Count For All Proteins The Occurance Of Modification Types
#'
#' Count for all protein 'sequ' the occurance of modification types defined in list 'modTyp' (only if in names(specAAMod)).
#' 
#' @param sequ (character) peptide sequence(s)
#' @param modTyp (list) modifications : $basMod for character vector of fixed modifications and $varMod for variable modifications. For one letter-code see AAfragSettings("modChem")
#' @param specAAMod (list) optional custom list showing which AA to be considered with which (one-letter) modification code (default \code{\link{AAfragSettings}})
#' @param knownMods (list) optional custom list showing which modification appears at what type of location, eg N-terminal, internal ... (default \code{\link{AAfragSettings}})
#' @param detailedCount (logical) 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) for bug-tracking: more/enhanced messages and intermediate objects written in global name-space  
#' @return This function returns a list of matrixes $cou and $combTerm, with number of modifications per peptides (line in 'pepTab') for basMod, varMod & varMo2
#' @seealso \code{\link{AAfragSettings}}, \code{\link{makeFragments}} 
#' @examples
#' protP2 <- c(mesp="MESPEPTIDES", pepe="PEPEPEP")
#' pepTab1 <- makeFragments(protTab=protP2, minFra=6, internFr=TRUE, massTy="mono")
#' modTy2 <- list(basMod=c("b","y","h"), varMod=c("x","p","o","q","e","j"))
#' .countModif(pepTab1[,2], modTyp=c("b","y"), specAAMod=AAfragSettings(outTy="all")$specAAMod, 
#'   knownMods=AAfragSettings(outTy="all")$knownMods)
#' @export
.countModif <- function(sequ, modTyp, specAAMod, knownMods, detailedCount=FALSE, silent=FALSE, debug=FALSE, callFrom=NULL){
  ## count for all protein 'sequ' the occurance of modification types defined in list 'modTyp' (only if in names(specAAMod))
  ## 'modTyp'.. character vector
  ## 'specAAMod','knownMods' for matching modification type to single letter AA code
  ## 'detailedCount' .. if TRUE return list of matrixes with counts per 'sequ' (rows) & all elements (AA) of each modTyp
  ##  (otherwise) default return matrix with 'sequ' (rows) and sum of counts per 'modTyp' (cols)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".countModif")
  if(is.list(modTyp)) {
    if(is.list(modTyp[[1]])) message(fxNa," argument 'modTyp' not designed as list of lists (will loose/append layers)")
    modTyp <- unlist(modTyp) }
  modTyp <- modTyp[which(modTyp %in% unlist(knownMods))]
  msg <- " 'modTyp' seems empty -nothing to count !"
  if(is.null(modTyp)) stop(fxNa,msg) else {if(all(modTyp== "")) stop(fxNa,msg)}
  chModTy <- modTyp %in% names(specAAMod)
  exclTy <- if(any(!chModTy)) modTyp[which(!chModTy)] else NULL
  addBl <- FALSE
  if(all(!chModTy) || length(sequ) <1) count <- matrix(0,nrow=length(sequ), ncol=length(chModTy), dimnames=list(sequ, modTyp)) else {
    modTyp <- modTyp[which(chModTy)]
    chM <- lapply(modTyp, function(x) {z <- specAAMod[[which(names(specAAMod) %in% x)]]; if(length(z) >0) z else NULL})     # AA-letters to consider for each modTyp
    names(chM) <- modTyp
    chM <- chM[which(sapply(chM,length) >0)]
    if(length(sequ) <2) { sequ <- c(sequ,"")}
    if(sequ[length(sequ)]=="" & length(sequ) >1) addBl <- TRUE                  # remove last sequence (since empty)
    if(all(sapply(chM,length)==0)) return(matrix(0, nrow=length(sequ) -addBl, ncol=length(modTyp), dimnames=list(sequ[1:(length(sequ)-addBl)],modTyp)))
    chN <- unlist(chM)
    names(chN) <- rep(names(chN), sapply(chN, length))
    count <- sapply(chN, function(x) .countLET(sequ, x, silent=debug))
    if(!detailedCount) { tmp <- matrix(nrow=length(sequ), ncol=length(chM), dimnames=list(sequ,names(chM)))
      preCol <- 0
      for(i in 1:length(chM)) {
        tmp[,i] <- if(length(chM[[i]]) <2) count[,preCol+1] else rowSums(count[,preCol+(1:length(chM[[i]]))])
          preCol <- preCol +length(chM[[i]])}
      count <- if(length(exclTy) >0) cbind(matrix(0, nrow=length(sequ), ncol=length(exclTy), dimnames=list(sequ,exclTy)),tmp) else tmp
      } }
  if(addBl) count <- matrix(count[1,], nrow=1, dimnames=list(sequ[1],colnames(count)))
  count }
  
  
#' Count Letters
#'
#' Count how many time a given letter occurs in each element of 'seq'
#' 
#' @param sequ (character) eg peptide sequence(s)
#' @param countCh (charcter) 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) for bug-tracking: more/enhanced messages and intermediate objects written in global name-space  
#' @return This function returns a numeric vector of counts for 'countCh' (single element !) in each element of 'seq'
#' @seealso \code{\link{AAfragSettings}}, \code{\link{makeFragments}} 
#' @examples
#' protP2 <- c(mesp="MESPEPTIDES", pepe="PEPEPEP")
#' .countLET(protP2,"P")
#' @export
.countLET <- function(sequ, countCh="K", silent=FALSE, debug=FALSE, callFrom=NULL){
  ## return numeric vector of counts for 'countCh' (single element !) in each element of 'seq'
  fxNa <- wrMisc::.composeCallName(callFrom,newNa=".countLET")
  if(length(countCh) >1) {countCh <- countCh[1]
    if(!silent) message(fxNa,"Trim argument 'countCh' to length 1 !")}
  x <- grep(countCh, sequ)
  out <- rep(0, length(sequ))
  if(length(x) >0) out[x] <- nchar(sequ[x]) - nchar(gsub(countCh,"",sequ[x]))
  names(out) <- sequ
  out }
    
