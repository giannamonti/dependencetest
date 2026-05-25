#!/usr/bin/env Rscript

##
# runner script example
#

devtools::load_all(".")

start_time <- Sys.time()

get_timestamp <- function() {
  curr_time <- Sys.time()
  elapsed_secs = curr_time - start_time
  list(
    timestamp = format(curr_time, "%Y%m%d-%H%M%S"),
    elapsed = sprintf("%.3f", elapsed_secs)
  )
}

log <- function(...) {
  msg <- paste("#", ..., sep = " ")
  ts <- get_timestamp()
  text <- paste(ts$timestamp, ts$elapsed, "LOG", msg, sep = "|" )
  message(text)
}

v <- function(...) cat(sprintf(...), "\n", sep = " ", file = stderr())
s <- function(...) do.call(paste, as.list(c(..., sep = ", ")))

scall <- function () {
  c(
    dmy_hello(),
    dmy_hello("Earth"),
    dmy_hello("Moon", "'Night")
  )
}

vcall <- function () {
  dmy_hello(c(
    "Mars",
    "Venus"
  ))
}

random_pause <- function () {
  sec <- as.integer(runif(1, max = 5))
  log(">>> pause:", sec, ", ...")
  Sys.sleep(sec)
  log("<<< pause:", sec, ", done")
}

task <- function() {
  v("scall: %s", s(scall()))
  random_pause()
  v("vcall: %s", s(vcall()))
  0
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  log("> start:", paste(args, collapse = ", "))
  log("? args:", paste(commandArgs(), collapse = " "))
  log("? wdir:", getwd())
  rc <- 0
  print(elapsed <- system.time({
    rc <- task()
  }))
  log("< end:", rc, " -- ", summary(elapsed))
  rc
}

main()
