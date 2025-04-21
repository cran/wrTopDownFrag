#' Evaluate Selected Lines Of PepTab (iso-mass) For Preferential Cutting Sites
#'
#' Evaluate selected lines of pepTab (iso-mass) for preferential cutting sites. Such sites are taken by default from \code{.prefFragPattern()} simplified from a publication 
#'  by the Kelleher group (Haverland 2017, J Am Soc Mass Spectrom) or can be furnished by the user.
#' @param z (matrix) main input, must contain cols specified as seqCol and "no","tailAA","precAA"
#' @param prefFragPat (matrix) specifies preferential fragmentation (which combination of AA to consider cols cTer,nTer,score), default made by \code{.prefFragPattern()}
#' @param seqCol (character) column names for the column containing the sequence to search for preferential cutting sites 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns line ID-numbers (pepTab[,"no"]) for those below median score (ie to remove from pepTab) or NULL if nothing to remove due to preferential fragmentation
#' @seealso \code{\link{makeFragments}}
#' @examples
#' peTab <- matrix(c("9","13","14","15", "LPVIAGHEAAG","PVIAGHEAAGI","EKKPFSI","KKPFSIE", 
#'   "P","L","E","E", "I","V","E","E"),nr=4,dimnames=list(NULL,c("no","seq","precAA","tailAA")))
#' evalIsoFragm(peTab)
#' @export
evalIsoFragm <- function(z, prefFragPat=NULL, seqCol="seq", silent=FALSE, debug=FALSE, callFrom=NULL){
  ## evaluate selected lines of pepTab (iso-mass) for preferential cutting sites
  ## 'z' .. matrix, must contain cols specified as seqCol and "no","tailAA","precAA"
  ## 'prefFragPat' .. matrix specifying preferential fragmentation (which combination of AA to consider cols cTer,nTer,score)
  ## return line ID-numbers (pepTab[,"no"]) for those below median score (ie to remove from pepTab) or NULL if nothing to remove due to preferential fragmentation
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="evalIsoFragm")
  if(is.null(prefFragPat)) prefFragPat <- .prefFragPattern()
  z <- as.matrix(z)
  nAA <- nchar(z[,seqCol])
  chRep <- duplicated(nAA, fromLast=FALSE) | duplicated(nAA, fromLast=TRUE)
  if(all(!chRep)) return(NULL) else {
    if(any(!chRep)) {z <- z[which(chRep),]; nAA <- nchar(z[,seqCol])}   # remove (unexpected case of) single instance of given AA length
    chLe <- table(nAA)
    if(length(chLe) >1) {
      z <- unlist(by(z, nAA, .evalIsoFra, prefFragPat=prefFragPat, seqCol=seqCol))
    } else .evalIsoFra(z, prefFragPat=prefFragPat, seqCol=seqCol)
  }}
  

#' Evaluate Selected Lines Of PepTab
#'
#' Evaluate selected lines of pepTab of SAME AA-length AND iso-mass for preferential cutting sites.
#' 
#' @param x (matrix) main input, must contain cols specified as seqCol and "no","tailAA","precAA"
#' @param prefFragPat (matrix) specifies preferential fragmentation (which combination of AA to consider cols cTer,nTer,score), default made by \code{.prefFragPattern()}
#' @param seqCol (character) column names for the column containing the sequence to search for preferential cutting sites 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns line ID-numbers (pepTab[,"no"]) for those below median score (ie to remove from pepTab)
#' @seealso \code{\link{makeFragments}}
#' @examples
#' peTab <- matrix(c("9","13","14","15", "LPVIAGHEAAG","PVIAGHEAAGI","EKKPFSI","KKPFSIE", 
#'   "P","L","E","E", "I","V","E","E"),nr=4,dimnames=list(NULL,c("no","seq","precAA","tailAA")))
#' .evalIsoFra(peTab)
#' @export
.evalIsoFra <- function(x, prefFragPat=NULL, seqCol="seq", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## evaluate selected lines of pepTab of SAME AA-length AND iso-mass for preferential cutting sites
  ## return line ID-numbers (pepTab[,"no"]) for those below median score (ie to remove from pepTab)
  x <- as.matrix(x)
  nAA <- nchar(x[,seqCol])
  prefSi <- if(is.null(prefFragPat)) .prefFragPattern() else prefFragPat
  sc <- rep(0, nrow(x))
  out <- NULL
  locSi <- paste0(x[,"precAA"],substr(x[,seqCol],1,1)) %in% paste0(prefSi[,1],prefSi[,2])       # N-terminal fragmentation sites
  if(any(locSi)) sc[which(locSi)] <- prefSi[match(paste0(x[which(locSi),"precAA"], substr(x[which(locSi),seqCol],1,1)), paste0(prefSi[,1],prefSi[,2])),3]
  nAA <- nchar(x[,seqCol])
  locSi <- paste0(substr(x[,seqCol], nAA, nAA), x[,"tailAA"]) %in% paste0(prefSi[,1],prefSi[,2])   # C-terminal fragmentation sites
  if(any(locSi)) sc[which(locSi)] <- sc[which(locSi)] + prefSi[match(paste0(substr(x[which(locSi), seqCol],nAA,nAA), x[which(locSi), "tailAA"]), paste0(prefSi[,1],prefSi[,2])), 3]
  if(any(sc >0)) {med <- stats::median(sc, na.rm=TRUE)            # evaluate scores for N & C-term
    goodSc <-  which(sc >= if(med ==0) 0.1 else med)
    if(length(goodSc) < length(sc)) out <- as.integer(x[-1*goodSc,1])}
  out }
  
  
#' Return data.frame with pattern of perferential fragmentation sites
#'
#' Return data.frame with pattern of perferential fragmentation sites
#' Here a simplified version (elaborate see Kelleher group: Haverland 2017, J Am Soc Mass Spectrom)
#' 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a data.frame with pattern of perferential fragmentation sites 
#' @seealso \code{\link{makeFragments}}
#' @examples
#' peTab <- matrix(c("9","13","14","15", "LPVIAGHEAAG","PVIAGHEAAGI","EKKPFSI","KKPFSIE", 
#'   "P","L","E","E", "I","V","E","E"),nr=4,dimnames=list(NULL,c("no","seq","precAA","tailAA")))
#' head(.prefFragPattern())
#' @export
.prefFragPattern <- function(silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## return data.frame with pattern of perferential fragmentation sites x|y (1st & 2nd col) and score (3rd col)
  ## here a simplified version (elaborate see Kelleher group: Haverland 2017, J Am Soc Mass Spectrom)
  AA <- wrProteo::AAmass()[1:20]            # so far exclude ornithine O & selenocysteine U
  pat <- data.frame(cTer=c(rep(c("D","E"), each=length(AA)), names(AA)[c(-4,-6)]), nTer=c(rep(names(AA),2), rep("P",length(AA)-2 )), score=0.5)
  ## include K, L, V ??
  pat[which((pat[,1]=="D" | pat[,1]=="E") & pat[,2]=="P"),3] <- 1   
  pat } 
   
