# =============================================================================
# Supplementary simulation: replicating the design of Genest, Neslehova,
# Remillard & Murphy (2019, Biometrika), Supplementary Material S4.4.
#
# Purpose: this is NOT a replacement for the contingency-table simulation
# in simulation_contingency.R. It is an additional, faithful replication of
# the design used by Genest et al. to assess the behaviour of Sn relative to
# chi-squared, G2 and the Bakirov-Rizzo-Szekely (BRS) test, with Bn and Pn
# added as new competitors. Reported regardless of outcome.
#
# Design (matching Genest et al. 2019, Supplementary Material):
#   N    = 1000 simulation runs
#   B    = 1000 multiplier bootstrap replicates (for Genest's Sn/Tn2)
#   n    = 100, 250
#   tau  = 0 (H0), 0.1, 0.2
#   copula family: Clayton, Gumbel
#   margins (paired, same margin on both variables):
#     P2  : Poisson(mean = 2)
#     Bin : Binomial(size = 4, p = 0.5)
#     Geo : Geometric(p = 0.5)
#     NB  : Negative Binomial(size = 5, p = 5/7)
#
# Tests compared:
#   Bn, Pn   – dependence package (basis = "dummy", on the resulting table)
#   Genest   – Sn / Tn2, MixedIndTests::TestIndCopula
#   BRS      – Bakirov et al. (2006), energy::indep.test(method = "mvI")
#   ChiSqr   – Pearson's chi-squared
#   G2       – likelihood-ratio test of independence
#
# NOTE: Zelterman's (1987) test was considered but deliberately omitted.
# Its bivariate independence-testing form requires the analytic mean and
# variance of D^2 under H0 for the asymptotic normal approximation, which
# could not be verified against a citable source. Rather than include an
# unverified ad hoc formula, it was dropped; G2 and chi-squared already
# cover the likelihood-based competitors used by Genest et al. (2019).
#
# -----------------------------------------------------------------------------
# IMPORTANT — fine-grained parallelization (read before changing N/chunk):
#
# A previous version assigned one entire (n, family, tau, margin) cell — all
# N = 1000 replications — to a single worker via parLapply. This is the same
# load-imbalance issue diagnosed in simulation_contingency.R: BRS
# (energy::indep.test, R = 199) and Genest (TestIndCopula, B = 1000) are both
# bootstrap-based and their cost can vary substantially across cells
# (e.g. sparse tables from Geo/NB margins at high tau). With one task = one
# full cell, the heaviest cell pins a single core for the entire run while
# lighter cells finish quickly and leave other cores idle — the job does not
# hang, it is bottlenecked on the single slowest worker.
#
# Fix: each task is now a (cell, chunk) pair — a cell split into chunks of
# `chunk_size` replications. All chunks across all cells are flattened into
# one task list distributed via parLapply, so the heaviest cell's workload
# is spread across many workers. Per-chunk p-values are re-aggregated by
# cell after the parallel run, before computing rejection rates.
# =============================================================================

library(copula)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dependence)
library(MixedIndTests)
library(energy)
library(parallel)

# -----------------------------------------------------------------------------
# Global constants
# -----------------------------------------------------------------------------
N          <- 1000     # simulation runs (matches Genest et al.: N = 1000)
chunk_size <- 50        # replications per task; tune so a chunk takes a few
# minutes even on the heaviest cell
n_v   <- c(100, 250)    # sample sizes, matching Genest et al. Tables S1-S2
tau_v <- c(0, 0.1, 0.2) # 0 = independence (level), else power
families <- c("clayton", "gumbel")
margins  <- c("P2", "Bin", "Geo", "NB")

# -----------------------------------------------------------------------------
# Directory setup
# -----------------------------------------------------------------------------
dir_results <- "results_genest_replication"
dir_plots   <- "plots_genest_replication"
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots,   showWarnings = FALSE, recursive = TRUE)

res_path  <- function(...) file.path(dir_results, ...)
plot_path <- function(...) file.path(dir_plots,   ...)

# -----------------------------------------------------------------------------
# Margin sampling functions (paired: same margin for X and Y)
# -----------------------------------------------------------------------------
margin_quantile <- function(u, margin) {
  switch(margin,
         P2  = qpois(u, lambda = 2),
         Bin = qbinom(u, size = 4, prob = 0.5),
         Geo = qgeom(u, prob = 0.5),
         NB  = qnbinom(u, size = 5, prob = 5/7),
         stop("Unknown margin: ", margin)
  )
}

# -----------------------------------------------------------------------------
# Data-generating function: draws n pairs (X, Y) with the requested copula
# (Clayton or Gumbel, parametrized by Kendall's tau) and the requested
# discrete margin applied to both X and Y via the quantile transform.
# tau = 0 returns independent uniforms (the H0 case).
# -----------------------------------------------------------------------------
gen_pair <- function(n, family, tau, margin) {
  if (tau <= 0) {
    u <- cbind(runif(n), runif(n))
  } else {
    cop <- switch(family,
                  clayton = claytonCopula(iTau(claytonCopula(), tau)),
                  gumbel  = gumbelCopula(iTau(gumbelCopula(),  tau)),
                  stop("Unknown family: ", family)
    )
    u <- rCopula(n, cop)
  }
  x <- margin_quantile(u[, 1], margin)
  y <- margin_quantile(u[, 2], margin)
  cbind(x, y)
}

# -----------------------------------------------------------------------------
# Pre-generate all simulation data IN MEMORY, split into chunks of
# `chunk_size` replications, and flatten into one long task list.
# Each task = one (n, family, tau, margin, chunk) combination.
# Note: tau = 0 is family-agnostic (pure independence), so it is generated
# only once per (n, margin) and labelled family = "H0" to avoid duplicating
# identical cells under both "clayton" and "gumbel".
# -----------------------------------------------------------------------------
cat("Generating simulation data in memory...\n")
set.seed(20240601)

cells <- list()
for (nk in n_v) {
  for (mg in margins) {
    cells[[length(cells) + 1]] <- list(n = nk, family = "H0", tau = 0, margin = mg)
    for (fam in families) {
      for (tv in tau_v[tau_v > 0]) {
        cells[[length(cells) + 1]] <- list(n = nk, family = fam, tau = tv, margin = mg)
      }
    }
  }
}

n_chunks  <- ceiling(N / chunk_size)
task_list <- vector("list", length(cells) * n_chunks)
task_idx  <- 1L

for (cell_idx in seq_along(cells)) {
  cl_ <- cells[[cell_idx]]
  nk  <- cl_$n
  fam <- cl_$family
  tv  <- cl_$tau
  mg  <- cl_$margin
  fam_gen <- if (fam == "H0") "clayton" else fam   # tau=0 makes family irrelevant
  
  for (ch in seq_len(n_chunks)) {
    this_chunk_size <- min(chunk_size, N - (ch - 1L) * chunk_size)
    
    x_mat <- matrix(NA_real_, nrow = this_chunk_size, ncol = nk)
    y_mat <- matrix(NA_real_, nrow = this_chunk_size, ncol = nk)
    
    for (i in seq_len(this_chunk_size)) {
      pair <- gen_pair(nk, fam_gen, tv, mg)
      x_mat[i, ] <- pair[, 1]
      y_mat[i, ] <- pair[, 2]
    }
    
    task_list[[task_idx]] <- list(
      cell_id = cell_idx, chunk_id = ch,
      n = nk, family = fam, tau = tv, margin = mg,
      x = x_mat, y = y_mat
    )
    task_idx <- task_idx + 1L
  }
}

cat("Data generation complete. ", length(task_list), "tasks across",
    length(cells), "cells (chunk size =", chunk_size, ").\n")

# -----------------------------------------------------------------------------
# Parallel cluster setup
# -----------------------------------------------------------------------------
n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
cat("Using", n_cores, "cores\n")

cl <- makeCluster(n_cores, rscript_args = "--no-init-file")

lib_path <- .libPaths()
clusterExport(cl, "lib_path")
clusterEvalQ(cl, .libPaths(lib_path))

clusterEvalQ(cl, {
  library(dependence)
  library(MixedIndTests)
  library(energy)
})

# -----------------------------------------------------------------------------
# G2 (likelihood-ratio) test implementation
# (not bundled in any of the packages already loaded, so implemented here)
# -----------------------------------------------------------------------------
g2_test <- function(x, y) {
  tab <- table(x, y)
  n_tot <- sum(tab)
  row_p <- rowSums(tab) / n_tot
  col_p <- colSums(tab) / n_tot
  expected <- outer(row_p, col_p) * n_tot
  obs <- as.vector(tab)
  exp_ <- as.vector(expected)
  keep <- obs > 0
  g2 <- 2 * sum(obs[keep] * log(obs[keep] / exp_[keep]))
  df <- (nrow(tab) - 1) * (ncol(tab) - 1)
  p_value <- pchisq(g2, df = df, lower.tail = FALSE)
  list(statistic = g2, p.value = p_value)
}

clusterExport(cl, "g2_test")

# -----------------------------------------------------------------------------
# Worker function: processes one (n, family, tau, margin, chunk) task,
# returning per-replication p-values (not yet aggregated to a rejection
# rate, since a cell's full N is now split across many tasks).
# Entire body wrapped in tryCatch; individual fragile test calls wrapped
# too so a single failure yields NA rather than crashing the worker.
# -----------------------------------------------------------------------------
run_chunk <- function(task) {
  tryCatch({
    
    nk  <- task$n
    fam <- task$family
    tv  <- task$tau
    mg  <- task$margin
    chunk_n <- nrow(task$x)
    
    pval_Bn     <- numeric(chunk_n)
    pval_Pn     <- numeric(chunk_n)
    pval_Genest <- numeric(chunk_n)
    pval_BRS    <- numeric(chunk_n)
    pval_ChiSqr <- numeric(chunk_n)
    pval_G2     <- numeric(chunk_n)
    
    for (i in seq_len(chunk_n)) {
      
      xi <- task$x[i, ]
      yi <- task$y[i, ]
      
      bn_pn <- tryCatch({
        indeptest(as.factor(xi), as.factor(yi), basis = "dummy")
      }, error = function(e) NULL)
      
      pval_Bn[i] <- if (!is.null(bn_pn)) bn_pn$B_pvalue else NA_real_
      pval_Pn[i] <- if (!is.null(bn_pn)) bn_pn$P_pvalue else NA_real_
      
      pval_Genest[i] <- tryCatch({
        out <- TestIndCopula(cbind(xi, yi), trunc.level = 2, B = 1000,
                             par = FALSE, graph = FALSE)
        out$pvalue$Tn2 / 100
      }, error = function(e) NA_real_)
      
      pval_BRS[i] <- tryCatch({
        indep.test(xi, yi, method = "mvI", R = 199)$p.value
      }, error = function(e) NA_real_)
      
      pval_ChiSqr[i] <- tryCatch({
        suppressWarnings(chisq.test(xi, yi)$p.value)
      }, error = function(e) NA_real_)
      
      pval_G2[i] <- tryCatch({
        g2_test(xi, yi)$p.value
      }, error = function(e) NA_real_)
    }
    
    list(
      cell_id = task$cell_id, n = nk, family = fam, tau = tv, margin = mg,
      pvals = data.frame(Bn = pval_Bn, Pn = pval_Pn, Genest = pval_Genest,
                         BRS = pval_BRS, ChiSqr = pval_ChiSqr, G2 = pval_G2)
    )
    
  }, error = function(e) {
    message("ERROR cell_id=", task$cell_id, " chunk_id=", task$chunk_id,
            " n=", task$n, " family=", task$family,
            " tau=", task$tau, " margin=", task$margin,
            ": ", conditionMessage(e))
    NULL
  })
}

clusterExport(cl, "run_chunk")

set.seed(20240602)
clusterSetRNGStream(cl, 20240602)

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
# Drop failed chunks, then re-aggregate p-values by cell before computing
# rejection rates, since a cell's N replications are now spread across
# multiple chunk results.
# -----------------------------------------------------------------------------
n_null <- sum(sapply(chunk_results, is.null))
if (n_null > 0L) warning(n_null, " chunks returned NULL and were dropped.")
chunk_results <- Filter(Negate(is.null), chunk_results)

cells_split <- split(chunk_results, sapply(chunk_results, function(r) r$cell_id))

test_cols <- c("Bn", "Pn", "Genest", "BRS", "ChiSqr", "G2")

results_list <- lapply(cells_split, function(chunks) {
  nk  <- chunks[[1]]$n
  fam <- chunks[[1]]$family
  tv  <- chunks[[1]]$tau
  mg  <- chunks[[1]]$margin
  
  pvals <- rbindlist(lapply(chunks, function(c) c$pvals))
  
  data.table(
    rejection_rate = sapply(test_cols, function(tc) mean(pvals[[tc]] < 0.05, na.rm = TRUE)),
    n = nk, family = fam, tau = tv, margin = mg, test = test_cols
  )
})

# -----------------------------------------------------------------------------
# Collect and save results
# -----------------------------------------------------------------------------
dt <- rbindlist(results_list)
fwrite(dt, res_path("genest_replication_results.csv"))
cat("Results saved to", res_path("genest_replication_results.csv"), "\n")

# -----------------------------------------------------------------------------
# Plotting: rejection rate vs tau, faceted by margin and n, one plot per family
# (H0 cells, tau = 0, are shown on every facet as the common baseline)
# -----------------------------------------------------------------------------
dt_h0  <- dt |> filter(family == "H0")
dt_dep <- dt |> filter(family != "H0")

for (fam in families) {
  dt_fam <- bind_rows(
    dt_h0  |> mutate(family = fam),
    dt_dep |> filter(family == fam)
  )
  dt_fam$n <- factor(paste0("n = ", dt_fam$n), levels = paste0("n = ", n_v))
  
  p <- dt_fam |>
    ggplot(aes(x = tau, y = rejection_rate,
               color = test, linetype = test, shape = test)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.05, lty = 2) +
    facet_grid(vars(n), vars(margin)) +
    labs(title = paste("Copula family:", fam),
         x = expression(tau), y = "Rejection rate") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(plot_path(paste0("genest_replication_", fam, ".pdf")),
         plot = p, width = 9, height = 5)
}

cat("Plots saved to", dir_plots, "\n")

# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------
# Note: family == "H0" (tau == 0) is the null case (level / size), not
# power; it is kept in the table below as a reference row — when reporting,
# the H0 rows should be read/labelled as empirical size.
mean_power <- dt |>
  group_by(n, family, tau, margin, test) |>
  summarise(mean_power = mean(rejection_rate),
            sd_power   = sd(rejection_rate),
            .groups    = "drop") |>
  arrange(n, family, tau, margin, desc(mean_power))

cat("\n--- Mean rejection rate by test, n, family, tau and margin ---\n")
cat("    (family = H0 is empirical size, not power)\n")
print(mean_power)
fwrite(mean_power, res_path("mean_power.csv"))

# Mean rank across (n, margin, tau) cells, separately by copula family
# (lower = better). Computed on dependence cells only (family != "H0"),
# since ranking "best under H0" is not a meaningful notion of power.
dt_power <- dt |> filter(family != "H0")

mean_ranks_list <- list()
for (fam in unique(dt_power$family)) {
  wdt <- dt_power |>
    filter(family == fam) |>
    select(-family) |>
    pivot_wider(names_from = test, values_from = rejection_rate)
  
  rdt <- t(apply(
    as.matrix(wdt[, test_cols]), 1,
    data.table::frankv, order = -1L, ties.method = "min"
  ))
  colnames(rdt) <- paste0("rank_", test_cols)
  mean_ranks_list[[fam]] <- sort(colMeans(rdt))
}

cat("\n--- Mean rank by test, Clayton copula, tau > 0 only (lower = better) ---\n")
print(mean_ranks_list[["clayton"]])
cat("\n--- Mean rank by test, Gumbel copula, tau > 0 only (lower = better) ---\n")
print(mean_ranks_list[["gumbel"]])

fwrite(as.data.frame(t(mean_ranks_list[["clayton"]])), res_path("mean_ranks_clayton.csv"))
fwrite(as.data.frame(t(mean_ranks_list[["gumbel"]])),  res_path("mean_ranks_gumbel.csv"))