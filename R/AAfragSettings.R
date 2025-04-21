#' Settings For AA Fragmentation
#'
#' This function provides basic settings for what types of fragments may accomodate which type of modifications :  $knownMods: information about which modifications may be considered, $specAAMod: specifc AA sites (if applicable),  $specAAMod: specifc AA sites (if applicable).
#' For example, here 'p' codes for gain of mass for HPO3 only at S, T and Y residues.
#' Note: $knownMods$Nterm and $knownMods$Cterm are treated as mutually exclusive
#' @param outTy (character) default "all" or any of the list-elements
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a list ($knownMods, $knspecAAMods, $modChem, $neutralLossOrGain)
#' @seealso  \code{\link{makeFragments}},  \code{\link{fragmentSeq}}, \code{\link[wrProteo]{massDeFormula}}
#' @examples
#' AAfragSettings()
#' @export
AAfragSettings <- function(outTy="all", silent=FALSE, debug=FALSE, callFrom=NULL) {
  knownMods <- list(noMod=c(""), Nterm=c("a","b","c"), Cterm=c("x","y","z"), NCterm=c("d","i"), ny=c(),
    intern=c("a","b","c","f","j","x","y","z"), spcAA=c("p","h","k","o","r","s","m","n","u"), spcNterm=c("g","n"), spcCterm=c("u"))
    specAAMod <- list(p=c("T","S","Y"), h=c("S","T","E","D"), k=c("Q","K","R","N"),o=c("M"), r=c("C"),s="S", m=c("R","K"),n=c("Q"), u=c("R"))
    ## if possible avoid overlapping 'intern' 'NCterm' overlapping, may cause problems in countPotMofifAAs() !!
  modChem <- matrix(c("","+0H","a","-C-2H-2O", "b","-2H-O","c","-O+N+2H", "x","+C+O-2H","y","","z","-N-2H", "d","-2H-O","i","-C-O",
    "f","-2H-O", "g","-C-2H-O-N", "j","-C-O","p","HP3O","h","-2H-O","k","-3H-N", "o","+O",
    "r","+O","s","+2C+2H+O","m","+C+2H","n","-3H-N","u","-C-4H-2N-O","q","-3H-P-4O"), ncol=2, byrow=TRUE, dimnames=list(NULL, c("ty","mod")))
  neutralLossOrGain <- c("b","d","f","h","k","n","z") 
  switch(outTy,
    knownMods=knownMods,
    specAAMod=specAAMod,
    modChem=modChem,
    neutralLossOrGain=neutralLossOrGain,
    all=list(knownMods=knownMods, specAAMod=specAAMod, modChem=modChem, neutralLossOrGain=neutralLossOrGain))} 
    
