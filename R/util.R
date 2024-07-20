#' Internal Printing Utility
#'
#' @description This internal function provides a flexible way to print messages
#' to the console, optionally including elapsed time information.
#'
#' @param verbose logical. If \code{TRUE}, messages are printed to the console. 
#'   If \code{FALSE}, messages are suppressed.
#' @param msg The message to be printed.
#' @param startTime A timestamp indicating the start time for 
#'   elapsed time calculation. If omitted, no elapsed time is shown.
#'
#' @return No return

#' @keywords internal
printer <- function(verbose, msg, startTime){
    if(verbose){
        if(missing(startTime)){
            message(msg)
        }
        else{
            message(sprintf('[Elapsed time: %.2fs] %s\n',
                        difftime(Sys.time(), startTime, units="secs"), msg))
        }
    }
}