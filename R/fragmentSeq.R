#' Fragment protein or peptide sequence
#'
#' Makes internal/terminal fragments of a SINGLE peptide/protein input (as single letter amino-acid code) and returns list of all possible sequences ($full, $Nter, $Cter, $inter).
#' 
#' @param sequ (character, length=1) sequence used for fragmenting, as as mono-aminoacid letter code (so that cuting will be perfomed between all the letters/characters)
#' @param minSize (integer) min number of AA residues for considering peptide fragments
#' @param maxSize (integer) max number of AA residues for considering peptide fragments
#' @param internFragments (logical) logical (return only terminal fragments if 'FALSE')
#' @param separTerm (logical) if 'TRUE', separate N-terminal, C-terminal and internal fragments in list
#' @param keepRedSeqs (logical) if 'FALSE' remove fragments with redundant content (but my be from different origin in 'sequ'); remove redundant so far only when no separation of Nterm/Cterm/intern as list
#' @param prefName (logical) alternative name for all fragments (default the sequence itself), avoid separators '.' and '-'
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return numeric vector with mass
#' @seealso \code{\link{makeFragments}};   \code{\link[wrProteo]{convAASeq2mass}}
#' @examples
#' fragmentSeq("ABCDE")
#' fragmentSeq("ABCDE", minSize=3, internFragments=FALSE)
#' fragmentSeq("ABCDE", minSize=3, internFragments=TRUE)
#' 
#' ## Run multiple peptides/proteins
#' twoPep <- cbind(c("a","ABCABCA"), c("e","EFGEFGEF"))
#' apply(twoPep, 2, function(x) fragmentSeq(x[2], mi=3, kee=FALSE, sep=TRUE, pre=x[1]))
#' 
#' ## Ubiquitin example
#' P0CG48 <- "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"  
#' system.time( fra1 <- (fragmentSeq(P0CG48, mi=5, kee=FALSE)))      # < 0.5 sec  
#'
#' @export
fragmentSeq <- function(sequ,minSize=3,maxSize=300,internFragments=TRUE,separTerm=FALSE,keepRedSeqs=TRUE,prefName=NULL,silent=FALSE,callFrom=NULL){
  ## make internal/terminal fragments as list ($full, $Nter, $Cter, $inter) of SINGLE input sequence 'sequ' and return as list
  ## needs wrMisc::firstOfRepeated()
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="fragmentSeq")
  tx <- c("argument '","minSize","' shoud be of length 1  (truncating !!)","sequ")
  if(length(sequ) >1) {if(!silent) message(fxNa,tx[c(1,4,3)]); sequ <- sequ[1]}
  if(length(minSize) <1) {minSize <- 3; if(silent) message(fxNa,"setting 'minSize' to default =3")}
  maxSize <- c(minSize, maxSize)
  minSize <- min(minSize, na.rm=TRUE)
  maxSize <- max(maxSize, na.rm=TRUE)
  if(is.null(prefName)) prefName <- if(length(unique(names(sequ)))==length(sequ)) names(sequ) else sequ
  if(nchar(sequ) < minSize & !silent) {
    message(fxNa," sequence given as 'sequ' already shorter than 'minSize'");return(NULL)}
  cut1 <- .termPepCut(sequ, mi=minSize, ma=maxSize ,sepNC=TRUE, mainName=prefName)
  ## now cut1 may be list -> force to list !!
  if(!is.list(cut1)) cut1 <- list(cut1)
  ## make internal fragments (run loop to reduce 'sequ' at both ends by 1 unit & re-run terminal fragments )
  if(internFragments) { nCha <- nchar(sequ) 
    frTo <- cbind(from=2:floor(nCha/2), to=(nCha-1):ceiling(1+nCha/2))
    frTo <- cbind(frTo,seqc=apply(frTo,1, function(x) substr(sequ,x[1],x[2])))
    chLe <- nchar(frTo[,"seqc"]) <minSize
    if(!all(chLe)) { if(any(chLe)) frTo <- frTo[which(!chLe),]
    cut1$inter <- if(nrow(frTo) >1) {unlist(apply(frTo[,c(3,1)],1, function(x) .termPepCut(x[1], mi=minSize, ma=maxSize, indexOffs=as.numeric(x[2])-1, mainName=prefName,sepNC=FALSE))) 
      } else .termPepCut(frTo[1,3], mi=minSize, ma=maxSize, indexOffs=as.numeric(frTo[1,1])-1, mainName=prefName, sepNC=FALSE) }}  
  nFrag <- sum(sapply(cut1,length))
  chRed <- unique(unlist(cut1))
  if(nFrag > length(chRed) & !silent) message(fxNa,nFrag- length(chRed)," out of ",if(nFrag >10e3)c(" ~",signif(nFrag,4)) else nFrag," fragments not unique")
  if(!keepRedSeqs) {
    uniq <- wrMisc::naOmit(match(chRed,unique(cut1$Nter)))
    if(length(cut1$Nter) >0) {
      if(length(uniq) <length(cut1$Nter)) cut1$Nter <- cut1$Nter[uniq]          # keep only non-redundant of Nter
      if(length(uniq) >0) chRed <- chRed[-1*uniq] }
    uniq <- wrMisc::naOmit(match(chRed, unique(cut1$Cter)))
    if(length(chRed) >0 & length(cut1$Cter) >0) {
      if(length(uniq) <length(cut1$Cter)) cut1$Cter <- cut1$Cter[uniq]          # keep only non-redundant of Cter
      if(length(uniq) >0) chRed <- chRed[-1*uniq] }
    uniq <- wrMisc::naOmit(match(chRed, unique(cut1$inter)))
    if(length(chRed) >0 & length(cut1$inter) >0) {
      if(length(uniq) <length(cut1$inter)) cut1$inter <- cut1$inter[uniq] }     # keep only non-redundant of inter
  }
  if(separTerm) cut1 else unlist(cut1) }   

#' @export
.termPepCut <- function(pe,mi,ma=1000,se1=".",se2="-",mainName=NULL,sepNC=FALSE,indexOffs=NULL) {
  ## make named character vector of sequential terminal fragments
  ## 'pe' .. single (!) peptide (character vector, length=1)
  ## 'mi','ma' .. min/max fragment length, should be <= length(pe) (otherwise the full length of 'pe' ALWATYS returned !)
  ## 'se1', 'se2' .. separators for adding numbers to specify partial/fragment locations
  ## 'sepNC' .. if TRUE, separate fragments from both ends as $Nter & $Cter in list
  ## 'indexOffs' .. offset to add for custom numbering in names (numeric, length=1), ie '1' will already increase by +1
  mi <- min(mi,nchar(pe))                                         # can't be shorter than 'pe'
  if(nchar(pe) <ma) ma <- nchar(pe)  
  if(is.null(mainName)) mainName <- pe
  indexOffs <- if(is.null(indexOffs)) 0 else as.numeric(indexOffs) 
  names(pe) <- paste(mainName,se1,indexOffs+1,se2,indexOffs+nchar(pe),sep="")
  if(mi==nchar(pe) & ma==nchar(pe)) return(if(sepNC) list(full=pe) else pe)
  x <- substring(pe,1,mi:min(ma,nchar(pe)-1))                               # N-term part
  chMa <- nchar(x) > ma
  if(any(chMa)) x <- x[which(!chMa)]
  y <- if(nchar(pe) > mi) substring(pe, (2:(nchar(pe)-mi+1)), nchar(pe)) else ""               # C-term part  
  chMa <- nchar(y) > ma
  if(any(chMa)) y <- y[which(!chMa)]
  basInd <- list(xL=1, xU=mi:min(nchar(pe)-1, ma), yL=max(nchar(pe)-ma+1,2):(nchar(pe)-mi+1), yU=nchar(pe))
  indexNa <- if(indexOffs==0) basInd  else lapply(basInd, function(x) x +indexOffs[1])
  names(x) <- paste(mainName,se1,indexNa[[1]],se2,indexNa[[2]],sep="")
  if(identical(y,"")) {if(sepNC) list(Nter=x) else x} else {
    names(y) <- if(identical(y,"")) "" else paste(mainName,se1,indexNa[[3]],se2,indexNa[[4]],sep="")
    fu <- if(nchar(pe) > ma) NULL else pe
    if(length(fu)>0) names(fu) <- paste(mainName,se1,1+indexOffs,se2,nchar(fu)+indexOffs,sep="")
    if(sepNC) list(full=fu, Nter=x, Cter=y) else c(fu, x, y)}
}

#' @export
.CtermPepCut <- function(pe,mi,se1=".",se2="-",mainName=NULL,indexOffs=NULL) {
  ## make named character vector of sequential terminal fragments
  ## 'pe' .. single peptide (character vector, length=1)
  ## 'mi' .. min fragment length
  ## 'se1', 'se2' .. separators for adding numbers to specify partial/fragment locations
  ## 'indexOffs' .. offset for custom numbering in names (numeric, length=1)
  mi <- min(mi,nchar(pe))                                         # can't be shorter than 'pe'
  y <- if(nchar(pe) > mi) substring(pe,(2:(nchar(pe)-mi+1)), nchar(pe)) else ""               # C-term part    
  namX <- paste(pe,se1,"1",se2,mi:nchar(pe),sep="")
  if(is.null(mainName)) mainName <- if(length(names(pe)) >0) names(pe)[1] else pe
  basInd <- list(xL=1, xU=mi:nchar(pe), yL=2:(nchar(pe) -mi +1), yU=nchar(pe))
  indexNa <- if(is.null(indexOffs)) basInd  else lapply(basInd, function(x) x +indexOffs[1])
  names(y) <- if(identical(y,"")) "" else paste(mainName,se1,indexNa[[3]],se2,indexNa[[4]],sep="")
  y }
  #

#' @export
.NtermPepCut <- function(pe,mi,se1=".",se2="-",mainName=NULL,sepNC=FALSE,indexOffs=NULL) {
  ## make named character vector of sequential terminal fragments
  ## 'pe' .. single peptide (character vector, length=1)
  ## 'mi' .. min fragment length
  ## 'se1', 'se2' .. separators for adding numbers to specify partial/fragment locations
  ## 'sepNC' .. if TRUE, separate fragments from both ends as $Nter & $Cter in list
  ## 'indexOffs' .. offset for custom numbering in names (numeric, length=1)
  ##   bugfix : won't return full lenfth query any more
  mi <- min(mi, nchar(pe))                                              # can't be shorter than 'pe'
  x <- substring(pe, 1, mi:(nchar(pe) -1))                              # N-term part
  namX <- paste(pe,se1,"1",se2,mi:(nchar(pe) -1), sep="")
  if(is.null(mainName)) mainName <- pe
  basInd <- list(xL=1, xU=mi:(nchar(pe)-1), yL=2:(nchar(pe)-mi), yU=nchar(pe))
  indexNa <- if(is.null(indexOffs)) basInd  else lapply(basInd, function(x) x+indexOffs[1])
  names(x) <- paste(mainName,se1,indexNa[[1]],se2,indexNa[[2]],sep="")
  x }
  
