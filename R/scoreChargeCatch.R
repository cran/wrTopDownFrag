#' Scoring of charge catching potential for peptides
#'
#' Make score based on cumulative search for AA with given potential to catch charge (H+, or optionally any charge).
#' Note : at current cumulative scoring large peptides may get priviliged.
#'
#' @param resTab (matrix or data.frame) matrix or data.frame of results for SINGLE protein (here only the column specified with argument 'pepCol' will be used)
#' @param pepCol (character) column name of 'resTab' containing the peptide sequence to be scored
#' @param scale01 (logical) linear rescale output to maximum 1.0
#' @param chargeMode (character) this value may be 'pos' (default) for the positively charged amino-acids K,R and H or, 
#'   if this argument has any other value, than all charged amino-acids (K,R,H, S,T,N,Q, D,E, W and Y) will be considered.
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced 
#' @return numeric vector with score for each peptide of resTab (even if \code{scale01=TRUE} minimum may be >0 if all peptides do contain charge-catching AAs)
#' @seealso \code{\link{fragmentSeq}}
#' @examples
#' resTa <- matrix(c(1:4,"PEPTID","PEPTIK","PEPTRK","AGV"), ncol=2,
#'   dimnames=list(NULL,c("predInd","seq"))) 
#' scoreChargeCatch(resTa)
#' 
#' @export
scoreChargeCatch <- function(resTab, pepCol="seq", scale01=TRUE, chargeMode="pos", silent=FALSE, callFrom=NULL) {
  ## Scoring of charge catching potential for peptides
  ## 
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="scoreChargeCatch")
  chatchCha <- sapply(.chargeCatchingAA(chargeMode=chargeMode)[,1], grep, resTab[,pepCol])
  chLi <- sapply(chatchCha,length) <1
  pepSco <- rep(0,nrow(resTab))
  names(pepSco) <- rownames(resTab)
  if(all(chLi)) {
    if(!silent) message(fxNa,"no charge catching AAs found in peptides")
  } else {  
    if(any(chLi)) chatchCha <- chatchCha[which(!chLi)] 
    for(i in 1:length(chatchCha)) {
      sco <- as.numeric(.chargeCatchingAA(chargeMode=chargeMode)[,2][which(names(chatchCha)[1] ==.chargeCatchingAA(chargeMode=chargeMode)[,1])])
      ## is additive function the right way ? (long peptides will be priviliged)
      pepSco[chatchCha[[i]]] <- pepSco[chatchCha[[i]]] +sco
      }}
  if(scale01) {
    maxSc <- max(pepSco,na.rm=TRUE)
    if(maxSc >0) pepSco <- round(pepSco/maxSc,3) } 
  pepSco } 
      
#' @export
.chargeCatchingAA <- function(chargeMode="pos"){
  ## produce matrix with values for capacity of catching (extra) charges
  chargeCatching <- if(identical(as.character(chargeMode),"pos")) {
    cbind(AA=c("K","R","H"), sco=rep(c(1),c(3)))
  } else {
    cbind(AA=c("D","E", "S","T","N","Q", "K","R","H", "W","Y"), sco=rep(c(1,1,0.7,0.7),c(2,4,3,2)))}  
  chargeCatching }   
