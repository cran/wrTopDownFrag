#' Corrective Values For Random Sequences For Mutations
#'
#' The protein sequence and composition of most proteomes cannot be well mimicked by pure random sequences.
#' Thus, random sequences (for comparing the quality of fitting) shhould be corrected accordingly, this vector provides help to do so.
#' This function loads corMutShift- or corInDelShift-values from RData. 
#' The vector contains possible mass alterations for random drawing either by mutating a given aminoacid (\code{}corMutShift),
#' or by making or in/del changes (\code{}corInDelShift). These values are based on simulations in the human proteome (from UniProt).
#' 
#' @param fi (character) file (and path) to RData to read as corMutShift or corInDelShift. The file when opend must contain a numeric vector either called 'corMutShift' or 'corInDelShift'.
#' @return This functions returns a numeric vector  (possible mass alterations for random drawing)
#' @seealso \code{\link{corInDelShift}}, \code{\link{randMassByStochastic}}
#' @examples
#' corMutShift <- corMutShift()
#' str(corMutShift)
#' 
#' corInDelShift <- corInDelShift()
#' str(corInDelShift)
#' @export
corMutShift <- function(fi=NULL) {
  ##
  if(is.null(fi)) fi <- file.path(system.file("extdata", package="wrTopDownFrag"), "corMutShift.RData")
  chFi <- file.exists(fi)
  if(chFi) try(load(fi), silent=TRUE) else stop("File ",fi," does not exist !")
  if("try-error" %in% class(chFi)) stop("Unable to load file '",fi,"', check rights to read or check if file is indeed Rdata")
  if(!"corMutShift" %in% ls()) stop("The file '",fi,"' does not contain a vector named 'corMutShift' !") 
  corMutShift }
 

#' Corrective Values For Random Sequences For In/Dels
#'
#' The protein sequence and composition of most proteomes cannot be well mimicked by pure random sequences.
#' Thus, random sequences (for comparing the quality of fitting) shhould be corrected accordingly, this vector provides help to do so.
#' This function loads corMutShift- or corInDelShift-values from RData. 
#' The vector contains possible mass alterations for random drawing either by mutating a given aminoacid (\code{}corMutShift),
#' or by making or in/del changes (\code{}corInDelShift). These values are based on simulations in the human proteome (from UniProt).
#' 
#' @param fi (character) file (and path) to RData to read as corMutShift or corInDelShift. The file when opend must contain a numeric vector either called 'corMutShift' or 'corInDelShift'.
#' @return This functions returns a numeric vector  (1907 possible mass alterations for random drawing)
#' @seealso  \code{\link{corMutShift}}, \code{\link{randMassByStochastic}}
#' @examples
#' corMutShift <- corMutShift()
#' str(corMutShift)
#' 
#' corInDelShift <- corInDelShift()
#' str(corInDelShift)
#' @export
corInDelShift <- function(fi=NULL) {
  ##
  if(is.null(fi)) fi <- file.path(system.file("extdata", package="wrTopDownFrag"), "corInDelShift.RData")
  chFi <- file.exists(fi)
  if(chFi) try(load(fi), silent=TRUE) else stop("File ",fi," does not exist !")
  if("try-error" %in% class(chFi)) stop("Unable to load file '",fi,"', check rights to read or check if file is indeed Rdata")
  if(!"corInDelShift" %in% ls()) stop("The file '",fi,"' does not contain a vector named 'corInDelShift' !") 
  corInDelShift }
  
