#' Full combinatorial and cumulative values
#'
#' Use for all preparing all combinations of non-compulsatory, ie variable, mass modifications
#' Variable modifications may or may not be present. Thus, for a given amino-acid with a variable modification two versions of the molecular weight need to be considered.
#' Most (variable) modifications are linked to a type of amino acid, like serine-residues for phosphorlylation.
#' Thus in this case, each instance of the amino acid in question may or may not be modified. 
#' So, for example if there are 2 serines, 0, 1 or 2 phosphorylation modifications may be present. 
#' For this reason the is the argument \code{nMax} to stay within biologically relevant ranges (external knowledge) and reduce complexity significantly.    
#' Some modifications are exclusive to others, argument \code{notSingle} : An (artificially occuring) de-phosphorylation event during fragmentation can only happen if the amino acid was already phosphorylated in the first place.  
#' 
#' 
#' @param nMax (integer or data.frame with 1 line) maximum number of modifications
#' @param modVal (numeric, has to have names !) the change of molecular mass introduced by given modifications (as specified by the name of the value)
#' @param notSingle (character) names of 'modVal' where 1st element of 'notSingle'  cannot happen/appear if 2nd element not present (eg de-phospho/phosphorylation)
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return named (concatenated names of modVal) numeric vector 
#' @seealso \code{\link[wrMisc]{convToNum}}
#' @examples
#' ## to follow easily the results, hypothetical mass-modification values were chosen
#' mo1 <- c(a=10, b=1, c=0.1, d=0.01); nMa1 <- c(1,2,0,3)
#' combinateAllAndSum(nMa1, mo1)
#' ## # like 'b' for phospho & 'd' for de-phospho (which can't happen without phospho event)
#' combinateAllAndSum(nMa1, mo1, notSingle=c("d","b"))    
#' @export
combinateAllAndSum <- function(nMax, modVal, notSingle=NULL, silent=TRUE, callFrom=NULL){
  ## full combinatorial and cumulative values : all values of 'modVal' up to 'nMax' times in all combinations (and summed at end)
  ## use for all combinations of non-compulsatory mass modifications
  ## return named (concatenated names of modVal) numeric vector
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="combinateAllAndSum")
  msg <- "'modVal' should be numeric, 'nMax' integer & of same length"
  if(is.data.frame(nMax)) if(ncol(nMax) ==length(modVal)) {
    nMaNa <- colnames(nMax)
    nMax <- as.integer(nMax)
    names(nMax) <- nMaNa }  
  if(!all(length(nMax) >= length(modVal), is.numeric(nMax), is.numeric(modVal))) stop(msg)
  if(length(nMax) > length(modVal)) {nMax <- nMax[which(nMax >0)]; if(length(nMax) != length(modVal)) stop(msg)}       # adjust if some nMax==1
  if(any(nMax <1)) {modVal <- modVal[which(nMax >0)]; nMax <- nMax[which(nMax >0)] }
  if(length(modVal) <1) return(0)
  if(sum(nMax >0) >1) {
    a3a <- wrMisc::combinatIntTable(nMax, include0=TRUE, asList=FALSE, silent=TRUE, callFrom=fxNa)
    if(length(notSingle)==2) {
      if(is.character(notSingle)) notSingle <- wrMisc::naOmit(match(notSingle, names(modVal)))
      if(!silent) message(fxNa," prohibit more occurances of '",names(modVal)[notSingle[1]],"' than '",names(modVal)[notSingle[2]],"'")
      resCo <- which(a3a[,,notSingle[1]] > a3a[,,notSingle[2]])
      if(length(resCo) >0) a3a[,,notSingle[1]][resCo] <- 0            # set to 0 (ie combination will become redundant, remove later)
    }
    useNa <- apply(a3a, c(1,2), function(x) paste(rep(names(modVal), x), collapse=""))     # for names
    out <- array(rep(modVal, each=prod(dim(a3a)[1:2])), dim=dim(a3a))*a3a                  # still array
    out <- as.numeric(apply(out,c(1,2), sum))
    names(out) <- useNa
    out <- out[wrMisc::firstOfRepeated(names(out))$indUniq]     # problem when 2 modifs with same mass change : filtering for unique value/number -> need to work as unique names !
    ## note: BUT working by names() will allow different type mass changes of same numeric value (if same AA concerned this can't be resolved at this place !)
  } else { out <- modVal[which(nMax >0)]*(0:max(nMax))
    names(out) <- sapply(0:max(nMax), function(x) paste(rep(names(modVal)[which(nMax >0)],x),collapse=""))
  }
  out }
  
#' @export
.parCombinateAllAndSum <- function(uniqCo,massModV,nProc=NULL,firstOfRepeated=NULL,parRegDefault=TRUE,silent=FALSE,callFrom=NULL){
  ## version for multi-processor execution
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".parCombinateAllAndSum")
  chPa <- requireNamespace("BiocParallel", quietly=TRUE)
  if(!chPa) stop(fxNa,": package 'BiocParallel' not installed, please install from Bioconductor")
  maxNProc <- 12
  nCoresAvail <- parallel::detectCores()
  nProcs <- if(is.null(nProc)) round(nCoresAvail*0.8) else min(nProc, nCoresAvail)          # all cores out of 2, n-1 when 4&6, n-2 at 10&12
  nProcs <- as.integer(nProcs)
  isWin <- length(grep("ming.32", R.Version()$platform)) > 0
  out <- NULL
  ## need to cut lines into named list for running bplapply() (instead of apply() over initial lines of uniqCo)
  uniqCoL <- by(uniqCo, 1:nrow(uniqCo), function(x) {names(x) <- colnames(uniqCo); x})    # cut in lines while keeping colnames
  ##
  fxPar <- function(x,massModV,combinateAllAndSum,fxNa,silent) {
    ## functions from outside (like combinateAllAndSum()) need to get 'imported' explicitetly
    combinateAllAndSum(x, massModV, notSingle=c("q","p"), callFrom=fxNa, silent=silent) }
  if(!silent) message(fxNa," ready to configure/launch as ",nProcs," processors (nProcs) ")  
  oldOp <- options()                       # for restoring from backup
  on.exit(options(oldOp))  
  ##based on comment from MMorgan remove quote within options
  options(MulticoreParam=BiocParallel::MulticoreParam(workers=nProcs))
  out <- BiocParallel::bplapply(uniqCoL, fxPar, massModV=massModV, combinateAllAndSum=combinateAllAndSum, fxNa=fxNa, silent=silent)
  BiocParallel::bpstop()  
  out }
    
