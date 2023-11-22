
#' Convert from basepairs to centimorgan
#'
#' Convert from basepairs to centimorgan
#'
#' @param pos_BP Numeric.  The position in basepairs.
#'
#' @return pos_CM The postion in centiMorgans
#' @keywords internal
#'
convert_BP_to_cM <- function(pos_BP){
  pos_BP/1000000

  }

#' Convert from centiMorgan to basepairs
#'
#' Convert from centiMorgan to basepairs
#'
#' @param pos_CM Numeric.  The position in centiMorgan.
#'
#' @return pos_BP The postion in basepairs
#' @keywords internal
convert_CM_to_BP <- function(pos_CM){

  pos_CM*1000000

  }
