#' Dummy hello function
#'
#' @param who String target
#' @param salutation String kind
#' @return A salutation string
#' @export
#' @examples
#' dmy_hello()
#' dmy_hello("Earth")
#' dmy_hello("Moon", "'Night")
dmy_hello <- function(who = "World", salutation = "Hello") {
  paste(salutation, " ", who, "!", sep = "")
}
