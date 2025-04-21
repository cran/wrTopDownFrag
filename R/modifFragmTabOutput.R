#' Change fragment identification output format (for biologists) 
#'
#' Change fragment identification output to format better adopted for biologists 
#'
#' @param datafr (data.frame) initial output from identifyPepFragments()
#' @param addData (matrix or data.frame) suppelemental data
#' @param fuseC (character) columns to exract preceeding and tailing AA to fuse with separator 'sep' to main sequence
#' @param sep (character) separator for concatenation
#' @param modifCol (character) default 'modif'
#' @param replMod (matrix) if names of modifications shoule be renamed : the columns 'old' and 'new' indicata how modifcations should be renamed 
#' @param finCols (character) columns to retain for final output
#' @param supFinCols (character)
#' @param sortTable (character) sort output 1st by name, then by 'beg' or 'end'
#' @param silent (logical) suppress messages
#' @param debug (logical) addtional diagnostic messages
#' @param callFrom (character) allow easier tracking of message produced
#' @return data.frame of reorganized identification results
#' @seealso \code{\link{identifyPepFragments}}
#' @examples
#' protP <- c(protP="PEPTIDE")
#' obsMassX <- cbind(a=c(199.1077,296.1605,397.2082,510.2922,625.3192),
#'   b=c(227.1026,324.1554,425.2031,538.2871,653.3141),
#'   x=c(729.2937,600.2511,503.1984,402.1507,289.0666),
#'   y=c(703.3145,574.2719,477.2191,376.1714,263.0874))
#' rownames(obsMassX) <- c("E","P","T","I","D")      # all 1 & 7 ions not included
#' modTy1 <- list(basMod=c("b","y"), varMod=c("p","o","q"))
#' frag1 <- identifyPepFragments(ex=as.numeric(obsMassX), pe=protP, modTy=modTy1, 
#'   minFragSize=2, chargeCatchFilter=FALSE)
#' (frag1b <- if(length(unlist(frag1$identif)) >0) modifFragmTabOutput(frag1))
#' @export
modifFragmTabOutput <- function(datafr, addData=NULL, fuseC=c("precAA","seq","tailAA"),sep=".",modifCol="mod",replMod=cbind(old="by",new="i"),
  finCols=c("fraNa","origNa","beg","end","seq","ty","mod","modSpec","obsMass","mass","ppmToPred","ambig","runNo","FDR","sco4","sc.prefFrag",
  "sc.chargeCatch","sc.complemFra","sc.sameSite","logInt"), 
  supFinCols=NULL, sortTable="end", silent=FALSE, debug=FALSE, callFrom=NULL){
  ## change output format for biologists 
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="modifFragmTabOutput")
  valid <- FALSE                      # initialize
  if(length(datafr) >0) {
    if(is.list(datafr) && "identif" %in% names(datafr)) datafr <- datafr$identif
    chDim <- dim(datafr)
    if(length(chDim) ==2) {if(chDim[1] >0 && chDim[2] >1) valid <- TRUE }}
  if(!valid) stop(fxNa,"No valid rasult from identifyPepFragments(), please provide data.frame or list with $identif which is NOT empty ...")  
  ## add supl data (ie more cols)
  if(length(addData) >0) { if(length(dim(addData)) <2) addData <- as.matrix(addData)
    if(nrow(datafr) == nrow(addData)) {
      datafr <- cbind(datafr, addData)
    } else {
      if(max(datafr[,"measLi"]) < nrow(addData)) {                           # not sure if this option fully correct !?!
        addData <- as.matrix(addData[as.integer(datafr[,"measLi"]),])        # try to match via index from 'measLi'
        if(is.null(colnames(addData))) colnames(addData) <- paste("suplData",1:ncol(addData),sep="_")
        chColNa <- colnames(addData) %in% colnames(datafr)
         if(any(chColNa)) message(fxNa, "change  colnames(addData)\n")
        if(any(chColNa)) colnames(addData)[which(chColNa)] <- paste(colnames(addData)[which(chColNa)],"sup",sep="_")
        supFinCols <- c(supFinCols, colnames(addData))
        datafr <- cbind(datafr,addData)
        rm(addData)}
  }} 
  chColNa <- match(finCols, colnames(datafr))
  names(chColNa) <- finCols
  if(debug) { message(fxNa, "xxMod0"); xxMod0 <- list(datafr=datafr,addData=addData,fuseC=fuseC,sep=sep,modifCol=modifCol,replMod=replMod,finCols=finCols,supFinCols=supFinCols,chColNa=chColNa)}
  if(any(is.na(chColNa))) { if(all(is.na(chColNa))) stop(fxNa,"NONE of the columns specified in 'finCols' found !!")
    if(!silent) message(fxNa,"TROUBLE AHEAD : ",sum(is.na(chColNa))," column(s) not found :",wrMisc::pasteC(finCols[which(is.na(chColNa))],quoteC="'"))}
  for(i in chColNa) if(is.factor(datafr[,i])) datafr[,i] <- as.character(datafr[,i]) 
  ## sort ascending by query protein & then by start site then by longest 
  if(sortTable=="beg" && nrow(datafr) >1) {
    fragLe <- nchar(datafr[,chColNa[4]]) - as.integer(datafr[,chColNa[3]]) + as.integer(datafr[,chColNa[2]])  # otherwise shortest will come 1st
    datafr <- datafr[order(datafr[,chColNa[1]], as.integer(datafr[,chColNa[2]]), fragLe, decreasing=FALSE),]}
  if(sortTable=="end" && nrow(datafr) >1) {
    fragLe <- as.integer(datafr[,chColNa[4]]) - as.integer(datafr[,chColNa[3]])   # otherwise shortest will come 1st: end-beg ;  include no of peptides ? nchar(datafr[,chColNa[5]]) 
    datafr <- datafr[order(datafr[,chColNa[1]], as.integer(datafr[,chColNa[3]]), fragLe, decreasing=TRUE),]}
  ## fusePrecAAtail
  if(debug) { message(fxNa, "xxMod1"); xxMod1 <- list(datafr=datafr,addData=addData,fuseC=fuseC,sep=sep,modifCol=modifCol,replMod=replMod,finCols=finCols,supFinCols=supFinCols,chColNa=chColNa)}
  if(length(fuseC) <3) if(!silent) message(fxNa,"TROUBLE AHEAD 'fuseC' should point to 3 columns existing in 'datafr', found only ",length(fuseC)) 
  if(length(fuseC) >3) {fuseC <- fuseC[1:3]; if(!silent) message(fxNa," 'fuseC' should point to 3 columns existing in 'datafr', found more ") }
  if(!is.integer(fuseC)) fuseC <- wrMisc::naOmit(match(fuseC,colnames(datafr)))
  if(debug) { message(fxNa, "xxMod2"); xxMod2 <- list(datafr=datafr,addData=addData,fuseC=fuseC,sep=sep,modifCol=modifCol,replMod=replMod,finCols=finCols,supFinCols=supFinCols)}
  if(length(fuseC) >2){
    datafr[,fuseC[2]] <- paste0(datafr[,fuseC[1]],sep,datafr[,fuseC[2]],sep,datafr[,fuseC[3]])
    datafr[,fuseC[2]] <- sub("\\.NA$","",sub("^NA\\.","",datafr[,fuseC[2]]))         # remove tailing .NA etc
    } else {if(!silent) message(fxNa,
      "Can't fuse preceeding AA, etc : 'fuseC' should point to 3 columns existing in 'datafr', found only ",length(fuseC))}
  ## change names of AA-modifications
  if(length(modifCol) >0) if(!is.integer(modifCol)) modifCol <- wrMisc::naOmit(match(modifCol,colnames(datafr)))
  if(length(modifCol) >0 && length(replMod) >0) {
    chReplMod <- lapply(replMod[,1], grep, datafr[,modifCol])           # which lines are subject to replace eg 'by' by 'i'
    for(i in which(sapply(chReplMod, length) >0)) {
      datafr[chReplMod[[i]], modifCol] <- sub(replMod[i,1],replMod[i,2], datafr[chReplMod[[i]],modifCol]) 
      rownames(datafr[chReplMod[[i]],]) <- sub(replMod[i,1],replMod[i,2],rownames(datafr[chReplMod[[i]],])) }}
  if(debug) { message(fxNa, "xxMod3"); xxMod3 <- list(datafr=datafr,addData=addData,fuseC=fuseC,sep=sep,modifCol=modifCol,replMod=replMod,finCols=finCols,supFinCols=supFinCols,chReplMod=chReplMod)}
  ## final selection of cols to report
  if(length(supFinCols) >0) finCols <- wrMisc::naOmit(unique(c(finCols,supFinCols)))
  chColNa <- match(finCols,colnames(datafr))
  if(!any(is.na(chColNa))) {if(!silent) message(fxNa,"Trouble with final selection of columns")}
  chColNa <- wrMisc::naOmit(chColNa)
  if(length(chColNa) >0) datafr <- datafr[,chColNa] else { dataFr <- NULL
    if(!silent) message(fxNa,"No valid input for reorganizing results-table !")}
  datafr }  
   
