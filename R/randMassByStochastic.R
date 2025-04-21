#' Make Decoy Mass By Full Randomization
#'
#' Make full random decoy mass vector (mimick pepTab)   
#' 
#' @param xChar (numeric vector) characterize main data, must conatain elements 'n','minV','maxV'
#' @param nRepeat (integer) number of time whole randomization process should be repeated
#' @param negAvoid (logical) if TRUE try avoiding 0 or negative random mass in result
#' @param setSeed (integer) seed for random number generation set.seed()
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages and objects exportet to current session for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a matrix with additional column 'decoyMass', or if 'sepCol'=FALSE as additional lines
#' @seealso \code{\link[stats]{Uniform}}, \code{\link{randMassByMut}}, \code{\link[wrMisc]{convToNum}}
#' @examples
#' rand <- randMassByStochastic(c(n=10, minV=2, maxV=7))
#' summary(rand)
#' @export
randMassByStochastic <- function(xChar, nRepeat=5, negAvoid=TRUE, setSeed=NULL, silent=FALSE,debug=TRUE,callFrom=NULL){
  ## make full random decoy mass vector (mimick pepTab)
  ## return input matrix with additional column 'decoyMass', or if 'sepCol'=FALSE as additional lines
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="randMassByStochastic")
  if(!all(c("n","minV","maxV") %in% names(xChar))) stop(fxNa," format of 'xChar' seems incorrect !")
  randRa <- xChar[c("minV","maxV")]
  if(!is.numeric(randRa)) {randRa <- as.numeric(as.character(randRa)); names(randRa) <- c("minV","maxV")}
  if(diff(randRa) <=0) stop(fxNa,"PROBLEM : min & max values of 'xChar' do NOT give DIFFERENCE >0 !! (values inverted ?)")
  if(negAvoid) { if(randRa[2] <0) message(fxNa,"Majority of data negative, impossible to avoid neagtive limits, ignoring argument"); negAvoid <- FALSE}
  if(negAvoid) randRa[1] <- max(randRa[1], 0, na.rm=TRUE)  
  if(length(setSeed)==1) if(is.finite(setSeed)) set.seed(setSeed)
  if(!silent) message(fxNa,"Make ",xChar["n"] *nRepeat," random values, ranging from ",signif(randRa[1],4)," to ",signif(randRa[2],4))
  rand <- stats::runif(as.integer(xChar["n"])*nRepeat, min=randRa[1], max=randRa[2])
  rand }
    
