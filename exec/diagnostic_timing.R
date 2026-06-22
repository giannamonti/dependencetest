# =============================================================================
# Diagnostic script: times Genest (TestIndCopula) and BRS (indep.test)
# separately, on ONE replication of EVERY (n, rho, table) combination, to
# find exactly which cells are pathologically slow before relaunching the
# full simulation.
#
# Run this on a single core, sequentially — it is meant to be fast (80
# combinations x 1 replication x 2 tests = 160 timed calls) and to produce
# a clear diagnostic table, not to be parallelized.
# =============================================================================

library(mvtnorm)
library(MixedIndTests)
library(energy)
library(data.table)

n_v <- c(25, 50, 100, 200)
r_v <- seq(0, 0.9, by = 0.1)

cont.table <- function(n, order, rho = 0.5) {
  mu    <- c(0, 0)
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  x     <- rmvnorm(n, mean = mu, sigma = Sigma)
  
  if (order == 3) {
    x.ord <- as.factor(ifelse(x[, 1] <= -0.6, 1, ifelse(x[, 1] > 0.6, 3, 2)))
    y.ord <- as.factor(ifelse(x[, 2] <= -0.6, 1, ifelse(x[, 2] > 0.6, 3, 2)))
  } else if (order == 5) {
    x.ord <- as.factor(ifelse(x[, 1] <= -1.0, 1,
                              ifelse(x[, 1] >   1.0, 5,
                                     ifelse(x[, 1] <= -0.3, 2,
                                            ifelse(x[, 1] <=  0.0, 3, 4)))))
    y.ord <- as.factor(ifelse(x[, 2] <= -1.0, 1,
                              ifelse(x[, 2] >   1.0, 5,
                                     ifelse(x[, 2] <= -0.3, 2,
                                            ifelse(x[, 2] <=  0.0, 3, 4)))))
  }
  rbind(x.ord, y.ord)
}

results <- data.table(
  n = integer(), rho = numeric(), table = character(),
  n_cells_used = integer(), min_cell_count = integer(),
  genest_time = numeric(), genest_status = character(),
  brs_time = numeric(), brs_status = character()
)

set.seed(999)

total <- length(n_v) * length(r_v) * 2
counter <- 0

for (nk in n_v) {
  for (rj in r_v) {
    for (ord in c(3, 5)) {
      
      counter <- counter + 1
      cat(sprintf("[%d/%d] n=%d rho=%.1f table=%dx%d ... ",
                  counter, total, nk, rj, ord, ord))
      
      xy  <- cont.table(nk, order = ord, rho = rj)
      xnum <- as.numeric(as.factor(xy[1, ]))
      ynum <- as.numeric(as.factor(xy[2, ]))
      
      tab <- table(xnum, ynum)
      min_cell <- min(tab[tab > 0])   # smallest non-zero cell count
      n_cells_used <- sum(tab > 0)
      
      # ---- Genest -------------------------------------------------------
      t0 <- proc.time()["elapsed"]
      genest_status <- "ok"
      genest_res <- tryCatch({
        out <- TestIndCopula(cbind(xnum, ynum), trunc.level = 2, B = 1000,
                             par = FALSE, graph = FALSE)
        out$pvalue$Tn2 / 100
      }, error = function(e) {
        genest_status <<- paste("ERROR:", conditionMessage(e))
        NA_real_
      })
      genest_time <- proc.time()["elapsed"] - t0
      
      # ---- BRS ------------------------------------------------------------
      t0 <- proc.time()["elapsed"]
      brs_status <- "ok"
      brs_res <- tryCatch({
        indep.test(xnum, ynum, method = "mvI", R = 199)$p.value
      }, error = function(e) {
        brs_status <<- paste("ERROR:", conditionMessage(e))
        NA_real_
      })
      brs_time <- proc.time()["elapsed"] - t0
      
      cat(sprintf("Genest=%.2fs BRS=%.2fs min_cell=%d\n",
                  genest_time, brs_time, min_cell))
      
      results <- rbind(results, data.table(
        n = nk, rho = rj, table = paste0(ord, "x", ord),
        n_cells_used = n_cells_used, min_cell_count = min_cell,
        genest_time = genest_time, genest_status = genest_status,
        brs_time = brs_time, brs_status = brs_status
      ))
    }
  }
}

fwrite(results, "diagnostic_timing_results.csv")

cat("\n\n=== SUMMARY: slowest 10 combinations by BRS time ===\n")
print(results[order(-brs_time)][1:10])

cat("\n=== SUMMARY: slowest 10 combinations by Genest time ===\n")
print(results[order(-genest_time)][1:10])

cat("\n=== Any errors? ===\n")
print(results[genest_status != "ok" | brs_status != "ok"])

cat("\n=== Correlation between min_cell_count and timing ===\n")
cat("cor(min_cell_count, brs_time)    =", cor(results$min_cell_count, results$brs_time), "\n")
cat("cor(min_cell_count, genest_time) =", cor(results$min_cell_count, results$genest_time), "\n")

cat("\nFull results saved to diagnostic_timing_results.csv\n")