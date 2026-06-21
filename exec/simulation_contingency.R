# =============================================================================
# Simulation study: contingency tables (3x3 and 5x5), nominal level / power
#
# Tests compared:
#   Bn, Pn    – dependence package (basis = "dummy")
#   ChiSqr    – Pearson's chi-squared test
#   Hoef      – Hoeffding's D (Hmisc)
#   Genest    – Cramér-von Mises copula-based test, Tn2 (MixedIndTests)
#   BRS       – Bakirov, Rizzo & Szekely (2006) test (energy::indep.test, mvI)
#
# Design:
#   n      = 25, 50, 100, 200
#   rho    = seq(0, 0.9, 0.1)
#   nsim   = 10,000 replications per (n, rho) cell
#   tables = 3x3 and 5x5 ordinal discretizations of a bivariate normal
#
# -----------------------------------------------------------------------------
# IMPORTANT — fine-grained parallelization (read before changing nsim/chunk):
#
# A previous version of this script assigned one entire (n, rho) cell — all
# nsim = 1e4 replications — to a single worker. This caused severe load
# imbalance: an isolated timing test showed that BRS (energy::indep.test,
# R = 199) costs ~2.8 seconds per replication on the heaviest cell (5x5,
# n = 200, rho = 0.9), i.e. ~7.8 hours for that cell ALONE on one core. Since
# parLapply assigns whole list elements to workers, that single cell pinned
# one core for days while lighter cells finished quickly and left other
# cores idle — the job did not hang, it was simply bottlenecked on the
# single slowest worker, and no amount of additional cores would have
# helped as long as one task = one full cell.
#
# Fix: each task sent to a worker is now a (cell, chunk) pair — a cell
# split into chunks of `chunk_size` replications. All chunks across all
# cells are flattened into one long task list and distributed via
# parLapply, so the heaviest cell's workload is spread across many workers
# instead of pinning a single one. Per-chunk results are re-aggregated by
# cell after the parallel run.
# =============================================================================

library(Hmisc)
library(mvtnorm)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(dependence)
library(MixedIndTests)
library(energy)
library(parallel)

# -----------------------------------------------------------------------------
# Global constants
# -----------------------------------------------------------------------------
nsim       <- 1e4
chunk_size <- 200     # replications per task; tune so a chunk takes a few
# minutes even on the heaviest cell (5x5, rho=0.9):
# ~2.8s/rep for BRS there => ~9.3 min per 200-rep chunk
n    <- c(25, 50, 100, 200)
r    <- seq(0, 0.9, by = 0.1)

# -----------------------------------------------------------------------------
# Directory setup
# -----------------------------------------------------------------------------
dir_results <- "results_contingency"
dir_plots   <- "plots_contingency"
dir.create(dir_results, showWarnings = FALSE)
dir.create(dir_plots,   showWarnings = FALSE)

res_path  <- function(...) file.path(dir_results, ...)
plot_path <- function(...) file.path(dir_plots,   ...)

# -----------------------------------------------------------------------------
# Data-generating function
# -----------------------------------------------------------------------------
# Generates an ordinal contingency table of size order x order from a
# bivariate normal with correlation rho, by discretizing each margin.
# Cutpoints:
#   order == 3: (-Inf, -0.6], (-0.6, 0.6], (0.6, +Inf)
#   order == 5: (-Inf, -1], (-1, -0.3], (-0.3, 0], (0, 1], (1, +Inf)
cont.table <- function(n, order, rho = 0.5) {
  mu    <- c(0, 0)
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  x     <- rmvnorm(n, mean = mu, sigma = Sigma)
  
  if (order == 3) {
    x.ord <- as.factor(ifelse(x[, 1] <= -0.6, 1,
                              ifelse(x[, 1] >   0.6, 3, 2)))
    y.ord <- as.factor(ifelse(x[, 2] <= -0.6, 1,
                              ifelse(x[, 2] >   0.6, 3, 2)))
  } else if (order == 5) {
    x.ord <- as.factor(ifelse(x[, 1] <= -1.0, 1,
                              ifelse(x[, 1] >   1.0, 5,
                                     ifelse(x[, 1] <= -0.3, 2,
                                            ifelse(x[, 1] <=  0.0, 3,
                                                   4)))))
    y.ord <- as.factor(ifelse(x[, 2] <= -1.0, 1,
                              ifelse(x[, 2] >   1.0, 5,
                                     ifelse(x[, 2] <= -0.3, 2,
                                            ifelse(x[, 2] <=  0.0, 3,
                                                   4)))))
  }
  rbind(x.ord, y.ord)
}

# -----------------------------------------------------------------------------
# Pre-generate all simulation data IN MEMORY, split into chunks of
# `chunk_size` replications, and flatten into one long task list.
# Each task = one (n, rho, chunk) combination, carrying its own pre-generated
# data for both the 3x3 and 5x5 tables.
# -----------------------------------------------------------------------------
cat("Generating simulation data in memory...\n")
set.seed(42)

grid        <- expand.grid(k = seq_along(n), j = seq_along(r))
n_chunks    <- ceiling(nsim / chunk_size)
task_list   <- vector("list", nrow(grid) * n_chunks)
task_idx    <- 1L

for (cell_idx in seq_len(nrow(grid))) {
  k  <- grid$k[cell_idx]
  j  <- grid$j[cell_idx]
  nk <- n[k]
  rj <- r[j]
  
  for (ch in seq_len(n_chunks)) {
    this_chunk_size <- min(chunk_size, nsim - (ch - 1L) * chunk_size)
    
    x3 <- matrix(NA_real_, nrow = this_chunk_size, ncol = nk)
    y3 <- matrix(NA_real_, nrow = this_chunk_size, ncol = nk)
    x5 <- matrix(NA_real_, nrow = this_chunk_size, ncol = nk)
    y5 <- matrix(NA_real_, nrow = this_chunk_size, ncol = nk)
    
    for (i in seq_len(this_chunk_size)) {
      xy3 <- cont.table(nk, order = 3, rho = rj)
      x3[i, ] <- as.numeric(as.factor(xy3[1, ]))
      y3[i, ] <- as.numeric(as.factor(xy3[2, ]))
      
      xy5 <- cont.table(nk, order = 5, rho = rj)
      x5[i, ] <- as.numeric(as.factor(xy5[1, ]))
      y5[i, ] <- as.numeric(as.factor(xy5[2, ]))
    }
    
    task_list[[task_idx]] <- list(
      cell_id = cell_idx, chunk_id = ch,
      n = nk, rho = rj,
      x3 = x3, y3 = y3, x5 = x5, y5 = y5
    )
    task_idx <- task_idx + 1L
  }
}

cat("Data generation complete. ", length(task_list), "tasks across",
    nrow(grid), "cells (chunk size =", chunk_size, ").\n")

# -----------------------------------------------------------------------------
# Parallel cluster setup
# -----------------------------------------------------------------------------
n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
cat("Using", n_cores, "cores\n")

cl <- makeCluster(n_cores, rscript_args = "--no-init-file")

# Make sure workers see the same library paths as the master (renv-safe)
lib_path <- .libPaths()
clusterExport(cl, "lib_path")
clusterEvalQ(cl, .libPaths(lib_path))

clusterEvalQ(cl, {
  library(Hmisc)
  library(dependence)
  library(MixedIndTests)
  library(energy)
})

# -----------------------------------------------------------------------------
# Worker function: processes one (n, rho, chunk) task for both 3x3 and 5x5
# tables, returning per-replication p-values (not yet aggregated to a
# rejection rate, since a cell's full nsim is now split across many tasks).
# Entire body wrapped in tryCatch; individual fragile test calls wrapped
# too so a single failure yields NA rather than crashing the worker.
# -----------------------------------------------------------------------------
run_chunk <- function(task) {
  tryCatch({
    
    nk         <- task$n
    rj         <- task$rho
    chunk_n    <- nrow(task$x3)
    
    Btest.pval3     <- numeric(chunk_n)
    Ptest.pval3     <- numeric(chunk_n)
    chisqtest.pval3 <- numeric(chunk_n)
    hoeffd.pval3    <- numeric(chunk_n)
    genest.pval3    <- numeric(chunk_n)
    brs.pval3       <- numeric(chunk_n)
    
    Btest.pval5     <- numeric(chunk_n)
    Ptest.pval5     <- numeric(chunk_n)
    chisqtest.pval5 <- numeric(chunk_n)
    hoeffd.pval5    <- numeric(chunk_n)
    genest.pval5    <- numeric(chunk_n)
    brs.pval5       <- numeric(chunk_n)
    
    for (i in seq_len(chunk_n)) {
      
      # ---- 3 x 3 table ------------------------------------------------------
      x3i <- task$x3[i, ]
      y3i <- task$y3[i, ]
      
      res3 <- indeptest(as.factor(x3i), as.factor(y3i), basis = "dummy")
      Btest.pval3[i]     <- res3$B_pvalue
      Ptest.pval3[i]     <- res3$P_pvalue
      chisqtest.pval3[i] <- chisq.test(as.factor(x3i), as.factor(y3i))$p.value
      hoeffd.pval3[i]    <- hoeffd(x3i, y3i)$P[1, 2]
      
      genest.pval3[i] <- tryCatch({
        out3 <- TestIndCopula(cbind(x3i, y3i),
                              trunc.level = 2, B = 1000,
                              par = FALSE, graph = FALSE)
        out3$pvalue$Tn2 / 100
      }, error = function(e) NA_real_)
      
      brs.pval3[i] <- tryCatch({
        indep.test(x3i, y3i, method = "mvI", R = 199)$p.value
      }, error = function(e) NA_real_)
      
      # ---- 5 x 5 table ------------------------------------------------------
      x5i <- task$x5[i, ]
      y5i <- task$y5[i, ]
      
      res5 <- indeptest(as.factor(x5i), as.factor(y5i), basis = "dummy")
      Btest.pval5[i]     <- res5$B_pvalue
      Ptest.pval5[i]     <- res5$P_pvalue
      chisqtest.pval5[i] <- chisq.test(as.factor(x5i), as.factor(y5i))$p.value
      hoeffd.pval5[i]    <- hoeffd(x5i, y5i)$P[1, 2]
      
      genest.pval5[i] <- tryCatch({
        out5 <- TestIndCopula(cbind(x5i, y5i),
                              trunc.level = 2, B = 1000,
                              par = FALSE, graph = FALSE)
        out5$pvalue$Tn2 / 100
      }, error = function(e) NA_real_)
      
      brs.pval5[i] <- tryCatch({
        indep.test(x5i, y5i, method = "mvI", R = 199)$p.value
      }, error = function(e) NA_real_)
    }
    
    list(
      cell_id = task$cell_id, n = nk, rho = rj,
      pvals3 = data.frame(Bn = Btest.pval3, Pn = Ptest.pval3,
                          ChiSqr = chisqtest.pval3, Hoef = hoeffd.pval3,
                          Genest = genest.pval3, BRS = brs.pval3),
      pvals5 = data.frame(Bn = Btest.pval5, Pn = Ptest.pval5,
                          ChiSqr = chisqtest.pval5, Hoef = hoeffd.pval5,
                          Genest = genest.pval5, BRS = brs.pval5)
    )
    
  }, error = function(e) {
    message("ERROR cell_id=", task$cell_id, " chunk_id=", task$chunk_id,
            " n=", task$n, " rho=", task$rho, ": ", conditionMessage(e))
    NULL
  })
}

clusterExport(cl, "run_chunk")

set.seed(20230101)
clusterSetRNGStream(cl, 20230101)

# -----------------------------------------------------------------------------
# Run parallel computation over the flattened (cell, chunk) task list
# -----------------------------------------------------------------------------
cat("Starting parallel computation on", n_cores, "cores,",
    length(task_list), "tasks...\n")
t_start <- proc.time()

chunk_results <- parLapply(cl, task_list, run_chunk)

t_elapsed <- proc.time() - t_start
cat(sprintf("Done. Wall time: %.1f min\n", t_elapsed["elapsed"] / 60))

stopCluster(cl)

# -----------------------------------------------------------------------------
# Drop failed chunks, then re-aggregate p-values by cell (n, rho) before
# computing rejection rates, since a cell's nsim replications are now
# spread across multiple chunk results.
# -----------------------------------------------------------------------------
n_null <- sum(sapply(chunk_results, is.null))
if (n_null > 0L) warning(n_null, " chunks returned NULL and were dropped.")
chunk_results <- Filter(Negate(is.null), chunk_results)

cells <- split(chunk_results, sapply(chunk_results, function(r) r$cell_id))

results_list <- lapply(cells, function(chunks) {
  nk <- chunks[[1]]$n
  rj <- chunks[[1]]$rho
  
  pvals3 <- rbindlist(lapply(chunks, function(c) c$pvals3))
  pvals5 <- rbindlist(lapply(chunks, function(c) c$pvals5))
  
  test_cols <- c("Bn", "Pn", "ChiSqr", "Hoef", "Genest", "BRS")
  
  rbind(
    data.table(
      rejection_rate = sapply(test_cols, function(tc) mean(pvals3[[tc]] < 0.05, na.rm = TRUE)),
      n = nk, rho = rj, table = "3 x 3", test = test_cols
    ),
    data.table(
      rejection_rate = sapply(test_cols, function(tc) mean(pvals5[[tc]] < 0.05, na.rm = TRUE)),
      n = nk, rho = rj, table = "5 x 5", test = test_cols
    )
  )
})

# -----------------------------------------------------------------------------
# Collect and save results
# -----------------------------------------------------------------------------
dt <- rbindlist(results_list)
fwrite(dt, res_path("rejection_rates_nominal.csv"))
cat("Results saved to", res_path("rejection_rates_nominal.csv"), "\n")

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
dt$n <- factor(paste0("n = ", dt$n),
               levels = paste0("n = ", c(25, 50, 100, 200)))

p_cont <- dt |>
  ggplot(aes(x = rho, y = rejection_rate,
             color = test, linetype = test, shape = test)) +
  facet_grid(vars(n), vars(table)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, lty = 2) +
  ylab("rejection rate") +
  xlab(expression(rho)) +
  theme_bw()

ggsave(plot_path("power_nominal_plot.pdf"), plot = p_cont, width = 7, height = 7)

cat("Plot saved to", plot_path("power_nominal_plot.pdf"), "\n")

# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------
# Note: rho = 0 is the null case (level / size), not power; it is kept in
# the table below as a reference row rather than split out separately —
# when reporting, the rho = 0 row should be read/labelled as empirical size.
test_cols <- c("Bn", "Pn", "ChiSqr", "Hoef", "Genest", "BRS")

# Mean rejection rate and SD by rho, test and table size
# (rho = 0 row = empirical size; rho > 0 rows = power)
mean_power <- dt |>
  group_by(table, rho, test) |>
  summarise(mean_power = mean(rejection_rate),
            sd_power   = sd(rejection_rate),
            .groups    = "drop") |>
  arrange(table, rho, desc(mean_power))

cat("\n--- Mean rejection rate by test, rho and table size ---\n")
cat("    (rho = 0 is empirical size, not power)\n")
print(mean_power)
fwrite(mean_power, res_path("mean_power.csv"))

# Mean rank across (n, rho) cells, separately for each table size
# (lower = better). Computed on rho > 0 only, since ranking "best under H0"
# is not a meaningful notion of power.
dt_power <- dt |> filter(rho > 0)

mean_ranks_list <- list()
for (tbl in unique(dt_power$table)) {
  wdt <- dt_power |>
    filter(table == tbl) |>
    select(-table) |>
    pivot_wider(names_from = test, values_from = rejection_rate)
  
  rdt <- t(apply(
    as.matrix(wdt[, test_cols]), 1,
    data.table::frankv, order = -1L, ties.method = "min"
  ))
  colnames(rdt) <- paste0("rank_", test_cols)
  mean_ranks_list[[tbl]] <- sort(colMeans(rdt))
}

cat("\n--- Mean rank by test, 3 x 3 table, rho > 0 only (lower = better) ---\n")
print(mean_ranks_list[["3 x 3"]])
cat("\n--- Mean rank by test, 5 x 5 table, rho > 0 only (lower = better) ---\n")
print(mean_ranks_list[["5 x 5"]])

fwrite(as.data.frame(t(mean_ranks_list[["3 x 3"]])), res_path("mean_ranks_3x3.csv"))
fwrite(as.data.frame(t(mean_ranks_list[["5 x 5"]])), res_path("mean_ranks_5x5.csv"))

# -----------------------------------------------------------------------------
# Quick check: n = 100 cross-tab
# -----------------------------------------------------------------------------
dt2 <- fread(res_path("rejection_rates_nominal.csv"))

print(
  dt2 |>
    filter(n == 100) |>
    pivot_wider(names_from = test, values_from = rejection_rate) |>
    arrange(table, rho)
)