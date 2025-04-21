#' Identifcation and scoring of preferential cuting sites  
#'
#' Search for preferential fragmentation sites from 'pepTab' among the 3 colums specified via 'useCol' (for full AA sequence, preceeding AA, tailing AA) 
#' and return sum of scores (from 3rd column of prefFragPat) for both ends.
#' Note : proteins must be witten as single lettre code.
#'
#' @param pepTab (matrix) peptide-fragments with lines for peptides, cols as sequence/preceedingAA/tailingAA
#' @param useCol (character) column names for peptide-sequence,preceeding and tailing AA
#' @param prefFragPat (matrix) for preferential fragmentation rules (see \code{.prefFragPattern()})
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages and objects exportet to current session for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a matrix with fragment sequence, mass, start- and end-position, heading and tailing AA (or NA if terminal fragment)
#' @seealso  \code{\link{makeFragments}}  
#' @examples
#' pepT <- cbind(precAA=c("A","D","D","A","D"),seq=c("AKA","PKA","AKA","PKD","PKD"),
#'   tailAA=c("A","A","D","P","P"))
#' scorePrefFrag(pepT)
#' @export
scorePrefFrag <- function(pepTab, useCol=c("seq","precAA","tailAA"), prefFragPat=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## search for preferential fragmentation sites from 'pepTab' among the 3 colums specified via 'useCol' (for full AA sequence, preceeding AA, tailing AA) 
  ## and return sum of scores (from 3rd column of prefFragPat) for both ends.
  ## Note : proteins must be witten as single lettre code
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="scorePrefFrag")
  if(length(dim(pepTab)) <2) stop(fxNa," 'pepTab' must be matrix or data.frame containing all comumns specified in 'useCol' !")
  errMsg <- " 'useCol' should identify which 3 columns to use from 'pepTab' (either column-names column-number)"
  if(length(useCol) <3) stop(fxNa, errMsg)
  if(debug) {message(fxNa, "xxSco1"); xxSco1 <- list(pepTab=pepTab,useCol=useCol,prefFragPat=prefFragPat)}
  if("index" %in% useCol) {if(!"index" %in% colnames(pepTab)) pepTab <- cbind(index=1:nrow(pepTab),pepTab) }
  if(is.character(useCol)) {
    chCol <- match(useCol,colnames(pepTab))
    if(any(is.na(chCol))) stop(fxNa," Column(s) ",wrMisc::pasteC(useCol[which(is.na(chCol))],quoteC="'")," not found in 'pepTab'\n",errMsg)
    useCol <- chCol
  } else if(max(useCol) >ncol(pepTab)) stop(fxNa,"Some values of 'useCol' seem to high'\n",errMsg )  
  pepTab <- as.matrix(pepTab[,useCol])
  useCol <- 1:length(useCol)
  if(is.null(prefFragPat)) prefFragPat <- .prefFragPattern()
  if(debug) { message(fxNa, "xxSco2"); xxSco2 <- list(pepTab=pepTab,prefFragPat=prefFragPat)}
  ## main
  pepLen <- nchar(pepTab[,useCol[1]])
  out <- cbind(Nt=paste(pepTab[,useCol[2]],substr(pepTab[,useCol[1]],1,1),sep="_"), Ct=paste(substr(pepTab[,useCol[1]], pepLen,pepLen), pepTab[,useCol[3]], sep="_"))
  ref <- paste(prefFragPat[,1],prefFragPat[,2],sep="_")
  out <- apply(out, 2, match, ref)
  out <- rowSums(matrix(as.numeric(prefFragPat[out,3]), ncol=2), na.rm=TRUE) /2
  chNA <- is.na(out)
  if(any(chNA)) out[which(chNA)] <- 0
  if(!is.null(rownames(pepTab))) names(out) <- rownames(pepTab)
  out }
   
