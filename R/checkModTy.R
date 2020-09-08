#' Check & complete mixed of variable and fixed modifications
#'
#' Check & complete settings for mixed of variable and fixed modifications.
#' The final format is a list with $basMod, $varMod and $varMo2 
#' 
#' @param modTy (character) list of modification types to be considered
#' @param knownMods (character) optonal custom list of known modifications, default from \code{AAfragSettings(outTy="all")$knownMods} 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return corrected list of mixed of variable and fixed modifications ($basMod, $varMod and $varMo2) 
#' @seealso \code{\link{AAfragSettings}}
#' @examples
#' modTy1 <- list(basMod=c("b","y","h"),varMod=c("p","o","q"))
#' checkModTy(modTy1)
#' @export
checkModTy <- function(modTy,knownMods=NULL,silent=TRUE,callFrom=NULL){
  ## check & complete
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="checkModTy")
  chMod <- which(names(modTy) %in% c("basMod","varMod","varMo2") & !(sapply(modTy, function(x) is.null(x) | identical(x,""))))
  if(length(chMod) <1) stop(" Problem with 'modTy' : either incorrect names or empty !")
  # check for repeated
  chRep <- lapply(modTy, duplicated)
  if(any(unlist(chRep))) for(i in which(sapply(chRep,sum) >0)) {modTy[[i]] <- unique(modTy[[i]])
    if(!silent) message(fxNa," correcting duplicated modification-terms to ",wrMisc::pasteC(chRep[[i]]))}
  ## check for unknown labels
  chMod <- lapply(modTy,function(x) x %in% unlist(knownMods))
  if(any(unlist(chMod)))  for(i in which(sapply(chMod, function(x) sum(!x)) >0)) {
    if(!silent) message(fxNa," removing unknown modification-labels ", wrMisc::pasteC(modTy[[i]][which(!chMod[[i]])], quoteC="'"))  
    modTy[[i]] <- modTy[[i]][which(chMod[[i]])]}  
  if(all(c("varMod","varMo2") %in% names(modTy) ==c(TRUE,FALSE))) {        # if modTy$varMo2 missing -> create new ...
    modTy$varMo2 <- modTy$varMod                                           # modTy$varMo2 for variable modifs really counted (ie wo 'q' since same mass as without 'p')
    chMod <- modTy$varMo2 %in% "q"
    if(!silent) message(fxNa," adding $varMo2 to 'modTy'")
    if(any(chMod) & !"p" %in% modTy$basMod) modTy$varMo2 <- modTy$varMo2[which(!chMod)] }
  modTy }

#' @export
.checkModTy <- function(modTy,knownMods,phoDePho=c("p","q"),modTyGr=c("basMod","varMod"),silent=FALSE,callFrom=NULL) {
  ## checking of 'modTy'
  ## return verified/corrected 'modTy'
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".checkModTy")
  chModFx <- function(mod, possMod=knownMods) {
    ch1 <- which(mod %in% unlist(possMod))
    if(length(ch1) >0) mod[ch1] else "" }
  if(length(unlist(modTy)) >0) {                               # if fragmentation/modification types given, check for known entries
  modTyIni <- modTy
  if(is.list(modTy)) if(any(modTyGr %in% names(modTy))) {      # clean modTy to known modifications only
    modTy <- list(c(),c())
    names(modTy) <- modTyGr
    if(modTyGr[1] %in% names(modTyIni)) {
      modTy[[1]] <- chModFx(modTyIni[[modTyGr[1]]])
      if(identical(phoDePho %in% modTy[[1]], c(FALSE,TRUE))) {
        if(!silent) message(callFrom,"de-phosphorylation without phosphorylation not realistic -> omit")
        modTy[[1]] <- modTy[[1]][which(!modTy[[1]] %in% phoDePho[2])] } }
    if(modTyGr[2] %in% names(modTyIni)) {
      modTy[[2]] <- chModFx(modTyIni[[modTyGr[2]]])
      ## include de-phosho when phospho in $varMod
      if(identical(phoDePho %in% modTy[[2]], c(TRUE,FALSE))) {
        if(!silent) message(callFrom,"add de-phosphorylation to optional modifications")
        modTy[[2]] <- c(modTy[[2]],phoDePho[2])}
      if(identical(phoDePho %in% modTy[[2]], c(FALSE,TRUE))) {
        if(!silent) message(callFrom,"add phosphorylation to optional modifications (since de-phospho found)")
        modTy[[2]] <- c(phoDePho[1],modTy[[2]])}
      modTy$varMo2 <- if(phoDePho[2] %in% modTy[[2]]) modTy[[2]][which(!modTy[[2]] ==phoDePho[2])] else modTy[[2]]  #variant: wo de-phospho for single modif
    }
  } else {modTy <- list(""); names(modTy) <- modTyGr[1];
    if(!silent) message(callFrom,"no fragmentation/modification types recognized, calculate as unmodified 'pep'")}}
  modTy }
        
