#' Identify Children/Parent Settings As a+b=c
#'
#' This functions helps identifying fragments ('parent') characterized by a start- and end-position, that got split into 2  'children' fragments. 
#' So, each one of the new 'children' conserves either the start- or end-site of the parent and the the remaining ends are on consecutive positions.
#' For example if the sequence 'BCDEFG' (parent) gets split into 'BCD' (positions 1-3) and 'EFG' (positions 4-6), 
#' this will be identified as a children/parent 'family' which could be represented as 'a+b=c' case.
#' Note : At this point only settings with 2 children are considered, for more complex scenarions one may build trees using \code{\link[wrMisc]{buildTree}} (however, this function does not identify 'parents').
#' In proteomics-applications some start- and end-sites may occur multiple times, representing eg unmodified and modified versions of the same basal peptide-sequence.
#' Such duplicated start- and end-cases are handeled as allowed, a 'child' (characterized by its start- and end-position) may occur multiple times, and the 
#' corresponding redundant rownames (eg peptide sequence like 'BCD') will be conserved. However, information reflecting eg different peptide modifications must be stored separately.
#' If redudant start- and end-sites accur with different row-names, repeated start- and end-sites will display \code{NA}.  
#'
#' @param fragments (matrix or data.frame) integer values in 1st column, for start site of fragment, and in 2nd column as end-sites of fragments, rownames as IDs
#' @param output (character) choose simply returning results as counts or as list with \code{$counts} and \code{$detailIndex} (list with details showing each child1,child2 & parent)   
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This functions returns either numeric vector with cumulated counts (corresponding to rows of \code{fragments}) or list with $count and $detailIndex (list with indexes refering to non-redundant entries of all a+b=c settings identified)  
#' @seealso \code{\link[wrMisc]{simpleFragFig}} ; for building longer consecutive trees (without identification of 'parent') \code{\link[wrMisc]{buildTree}}
#' @examples
#' frag3 <- cbind(beg=c(4,2,3,7,13,13,15, 2,9,2,9), end=c(14,6,12,8,18,20,20, 8,12,12,18)) 
#' rownames(frag3) <- c("K","A","E","B","C","D","F", "H","G","I","J")
#' countChildrenParent(frag3)
#' ## example with duplicate start- and end-position positions
#' frag3c <- cbind(beg=c(4,2,3,7, 7,13, 13,13,15, 2,9,2,9,9),
#'   end=c(14,6,12,8, 8,18, 18,20,20, 8,12,12,12,18))
#' rownames(frag3c) <- c("K","A","E", "B","B", "C","C","D","F", "H","G","I","G","J")
#' countChildrenParent(frag3c, out="det")
#' 
#' @export
countChildrenParent <- function(fragments, output="count", silent=FALSE, debug=FALSE, callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="countChildrenParent")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  msg <- "expecting matrix (or data.frame) with 2 columns with integer values for start end end-sites !" 
  if(length(fragments) <1 || length(dim(fragments)) <2) stop(msg)
  if(is.data.frame(fragments)) fragments <- as.matrix(fragments)
  if(is.null(rownames(fragments))) rownames(fragments) <- 1:nrow(fragments)
  ## need to remove duplicate entries 
  chDu <- duplicated(paste0(fragments[,1],fragments[,2]), fromLast=FALSE)
  iniP <- rownames(fragments)
  if(any(chDu)) {
    fragments <- fragments[which(!chDu),] 
    if(!silent) message(fxNa, "Remove ",sum(chDu)," duplicate entries")  }
  if(!is.numeric(fragments)) { fragments <- matrix(wrMisc::convToNum(fragments,spaceRemove=TRUE,convert=NULL,
    callFrom=fxNa, silent=silent), ncol=ncol(fragments), dimnames=dimnames(fragments))}
    if(!is.numeric(fragments)) stop(msg)
    if(any(is.na(fragments))) { fragments <- fragments[which(rowSums(is.na(fragments)) >0),]
      if(!silent) message(fxNa,"Expecting 'fragments' without NAs ! (remove lines)")
      if(length(dim(fragments)) <2) fragments <- matrix(fragments, ncol=2)}
  out <- rep(0,nrow(fragments))
  names(out) <- rownames(fragments)
  ## strategy : parent was broken into 2 children, thus one child must have another continuing child (fragment)
  ## parents must duplicate either start or end-site of children, use to pre-filter data
  ## filter for duplicated start sites  (later among those passing filter for continuing children/fragments)
  chDuSt1 <- duplicated(fragments[,1], fromLast=TRUE) 
  chDuSt2 <- duplicated(fragments[,1], fromLast=FALSE) 
  uniSt <- !chDuSt1 & !chDuSt2
  names(uniSt) <- rownames(fragments)
  chCont3 <- NULL
  if(any(!uniSt)) {
    ## parents do exist, now search for 'continuing other child'
    st <- which(!uniSt)
    if(length(st) >length(unique(st))) message(fxNa," (potential) Problem with repeatedly found start fragm !")
    ## search for continuing children/fragments (ie after end) 
    chCont2 <- lapply(st, function(x) {y <- x; names(y) <- rownames(fragments)[x]; c(y, which(fragments[x,2] +1 ==fragments[,1]))})
    names(chCont2) <- st
    chLe <- sapply(chCont2,length) <2
    if(any(chLe)) chCont2 <- chCont2[which(!chLe)]
    ## search for parent fragment
    if(any(!chLe)) {chCont3 <- lapply(chCont2, function(x) wrMisc::naOmit(as.integer(c(x[1:2],
      sapply(x[-1], function(y) which(fragments[x[1],1] ==fragments[,1] & fragments[y,2] ==fragments[,2]))))))
      chL2 <- sapply(chCont3,length) <3
      chCont3 <- if(any(!chL2)) chCont3[which(!chL2)] else NULL
      ## add names (if full output has been chosen) 
      if(length(chCont3) >0 && !identical(output,"count")) for(i in 1:length(chCont3)) names(chCont3[[i]]) <- rownames(fragments)[chCont3[[i]]]
    }
  }
  if(length(chCont3) >0) {
    tab <- table(unlist(chCont3))
    out[as.integer(names(tab))] <- tab  
  }
  ## need to reintroduce removed duplicates
  if(any(chDu)) {
    tmp <- rep(NA,length(iniP))
    names(tmp) <- iniP
    out <- wrMisc::getValuesByUnique(tmp, out, silent=silent, callFrom=fxNa)
  }  
  if(identical(output,"count")) out else list(count=out, detailIndex=chCont3) }
   
