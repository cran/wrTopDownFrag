#' Make decoy mass by full randomization
#'
#' Make full random decoy mass vector (mimick pepTab)   
#' 
#' @param pepTab (matrix) typically table of petides, one column should match 'useCol' for numeric mass values
#' @param randCha (numeric) vector of possible mass alterations for random drawing
#' @param useCol (character)  column from 'pepTab' with mass values to make decoys
#' @param negAvoid (logical) if TRUE try avoiding 0 or negative random mass in result
#' @param sepCol (character) optional column from 'pepTab' 
#' @param inDel (logical) switch to make random mass by insertion/deletion of 1 AA  
#' @param nAlter (integer) number of alterations per peptide  
#' @param setSeed (character) seed for random number generation set.seed()
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages and objects exportet to current session for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a matrix with additional column 'decoyMass', or if 'sepCol'=FALSE as additional lines
#' @seealso \code{\link{randMassByStochastic}}
#' @examples
#' pepTab1 <- cbind(no=11:12, seq=c("YVVDTS","YVVDTSK"), origNa="test.P000",
#'   ty=c("inter","Cter"),mass=c(681.308997360991,809.403960378691))
#' randMassByMut(pepTab1, corMutShift())
#' randMassByMut(pepTab1, corInDelShift(), inDel=TRUE)
#' @export
randMassByMut <- function(pepTab, randCha, useCol="mass", negAvoid=TRUE, sepCol=FALSE, inDel=FALSE, nAlter=1, setSeed=NULL, silent=FALSE, debug=FALSE, callFrom=NULL){
  ## make decoy mass vector by introducing random mutations or in/del varaiants to sequences of pepTab
  ## return matrix similar to input but with additional column 'decoyMass', or if 'sepCol'=FALSE as additional lines
  ## note : increases slightly resultant masses (due to adjustment for init freq of AAs)
  ## 'pepTab'.. (matrix), typically table of petides, one column should match 'useCol' for numeric mass values
  ## 'useCol'.. (character) column from 'pepTab' with mass values to make decoys
  ## 'randCha'.. (numeric)  vector of possible mass alterations for random drawing
  ## 'inDel'.. (logical) to change pure adding of random change (with delta-mass) to add&remove (with insertions/deletions)n will override 'negAvoid' to FALSE
  ## 'nAleter'.. (integer) define numer of changes on given line/value from 'pepTab' (default=1)
  ## 'setSeed'.. (integer) seed for defining random drawings, use via set.seed()
  ## 'sepCol'.. (logical)  if TRUE decoy masses will be added as new column 'decoyMass'
  ## 'negAvoid'.. (logical) if TRUE try avoiding 0 or negative random mass in result
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="randMassByMut")
  if(length(randCha) >0) { randCha <- round(stats::rnorm(1900, sd=30),3)  # mimic precise list of changes
    randCha <- randCha + 3.283 -mean(randCha, na.rm=TRUE)
    ch0 <- randCha ==0 | randCha > 129.1 | randCha < -129.1               # stay within realistic range (max mass change at mutation=129.06)
    if(any(ch0)) randCha <- randCha[which(!ch0)]
    }
  ind <- wrMisc::naOmit(randCha)
  if(length(setSeed)==1) if(is.finite(setSeed)) set.seed(setSeed)
  ind <- sample.int(length(ind), nrow(pepTab), replace=TRUE)
  if(debug) { silent <- FALSE
    message(fxNa,"Change mass by ",wrMisc::pasteC(utils::head(randCha[ind]) ))}
  orient <- if(inDel) rep(c(1,-1),ceiling(nrow(pepTab)/2))[1:nrow(pepTab)] else rep(1,nrow(pepTab))
  decM <- as.numeric(pepTab[,useCol]) +randCha[ind]
  if(nAlter >1) for(i in 2:nAlter) decM <- as.numeric(pepTab[,useCol]) + orient *randCha[ind[(i-1)*nrow(pepTab) +(1:nrow(pepTab))]]
  if(inDel) negAvoid <- FALSE 
  chNeg <- decM <= 0
  if(any(chNeg) && negAvoid) decM[which(chNeg)] <- as.numeric(pepTab[which(chNeg),useCol]) +abs(randCha[ind][which(chNeg)])
  if(!silent) message(fxNa,"Decoy same as input: ",sum(pepTab[,useCol] ==decM),";  any identical mass (input): ",sum(table(pepTab[,useCol]) >1),";  any identical mass to decoy: ", sum(pepTab[,useCol] %in% decM))
  if(sepCol) cbind(pepTab, decoyMass=decM) else {
    tmp <- matrix(nrow=length(decM), ncol=ncol(pepTab), dimnames=list(paste("decoy",1:length(decM), sep="_"), colnames(pepTab)))
    if("no" %in% colnames(pepTab)) tmp[,"no"] <- max(as.numeric(pepTab[,"no"]),na.rm=TRUE) + 1:length(decM)
    if(debug) message(fxNa," dim pepTab",wrMisc::pasteC(dim(pepTab)),"   dim tmp",wrMisc::pasteC(dim(tmp)),"   max n",max(as.numeric(pepTab[,"no"]),na.rm=TRUE))
    chNa <- c("origNa","ty") %in% colnames(pepTab)
    if(any(chNa)) tmp[,wrMisc::naOmit(match(c("origNa","ty"), colnames(tmp)))] <- "decoy"
    tmp[,useCol] <- decM
    rbind(pepTab, tmp)} }
 
