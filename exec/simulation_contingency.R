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
# Parallelization mirrors simulation_power_parallel.R:
#   - all data pre-generated in memory and embedded in the parameter list
#   - workers started with --no-init-file to bypass renv/.Rprofile issues
#   - library paths exported explicitly to workers
#   - entire worker body wrapped in tryCatch; TestIndCopula / indep.test
#     calls wrapped individually so a single failure yields NA, not a crash
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
nsim <- 1e4
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
#   order == 4: (-Inf, -0.8], (-0.8, 0], (0, 0.8], (0.8, +Inf)
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
  } else if (order == 4) {
    x.ord <- as.factor(ifelse(x[, 1] <= -0.8, 1,
                              ifelse(x[, 1] >   0.8, 4,
                                     ifelse(x[, 1] >   0.0, 3, 2))))
    y.ord <- as.factor(ifelse(x[, 2] <= -0.8, 1,
                              ifelse(x[, 2] >   0.8, 4,
                                     ifelse(x[, 2] >   0.0, 3, 2))))
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
# Pre-generate all simulation data IN MEMORY and embed into param_list
# Each cell (n_k, rho_j) gets its own pre-generated set of nsim tables for
# both the 3x3 and 5x5 discretizations, stored as numeric matrices.
# -----------------------------------------------------------------------------
cat("Generating simulation data in memory...\n")
set.seed(42)

grid       <- expand.grid(k = seq_along(n), j = seq_along(r))
param_list <- vector("list", nrow(grid))

for (idx in seq_len(nrow(grid))) {
  k  <- grid$k[idx]
  j  <- grid$j[idx]
  nk <- n[k]
  rj <- r[j]
  
  # Pre-generate nsim tables for 3x3 and 5x5, store as numeric (not factor)
  x3 <- matrix(NA_real_, nrow = nsim, ncol = nk)
  y3 <- matrix(NA_real_, nrow = nsim, ncol = nk)
  x5 <- matrix(NA_real_, nrow = nsim, ncol = nk)
  y5 <- matrix(NA_real_, nrow = nsim, ncol = nk)
  
  for (i in seq_len(nsim)) {
    xy3 <- cont.table(nk, order = 3, rho = rj)
    x3[i, ] <- as.numeric(as.factor(xy3[1, ]))
    y3[i, ] <- as.numeric(as.factor(xy3[2, ]))
    
    xy5 <- cont.table(nk, order = 5, rho = rj)
    x5[i, ] <- as.numeric(as.factor(xy5[1, ]))
    y5[i, ] <- as.numeric(as.factor(xy5[2, ]))
  }
  
  param_list[[idx]] <- list(
    k = k, j = j, n = nk, rho = rj,
    x3 = x3, y3 = y3, x5 = x5, y5 = y5
  )
}

cat("Data generation complete.\n")

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
# Worker function: processes one (n, rho) cell for both 3x3 and 5x5 tables
# Entire body wrapped in tryCatch; individual test calls wrapped too so a
# single failed replication yields NA rather than crashing the worker.
# -----------------------------------------------------------------------------
run_one <- function(params) {
  tryCatch({
    
    nk         <- params$n
    rj         <- params$rho
    nsim_local <- nrow(params$x3)
    
    Btest.pval3     <- numeric(nsim_local)
    Ptest.pval3     <- numeric(nsim_local)
    chisqtest.pval3 <- numeric(nsim_local)
    hoeffd.pval3    <- numeric(nsim_local)
    genest.pval3    <- numeric(nsim_local)
    brs.pval3       <- numeric(nsim_local)
    
    Btest.pval5     <- numeric(nsim_local)
    Ptest.pval5     <- numeric(nsim_local)
    chisqtest.pval5 <- numeric(nsim_local)
    hoeffd.pval5    <- numeric(nsim_local)
    genest.pval5    <- numeric(nsim_local)
    brs.pval5       <- numeric(nsim_local)
    
    for (i in seq_len(nsim_local)) {
      
      # ---- 3 x 3 table ------------------------------------------------------
      x3i <- params$x3[i, ]
      y3i <- params$y3[i, ]
      
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
        indep.test(x3i, y3i, method = "mvI")$p.value
      }, error = function(e) NA_real_)
      
      # ---- 5 x 5 table ------------------------------------------------------
      x5i <- params$x5[i, ]
      y5i <- params$y5[i, ]
      
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
        indep.test(x5i, y5i, method = "mvI")$p.value
      }, error = function(e) NA_real_)
    }
    
    rbind(
      data.frame(rejection_rate = mean(Btest.pval3     < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "3 x 3", test = "Bn"),
      data.frame(rejection_rate = mean(Ptest.pval3     < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "3 x 3", test = "Pn"),
      data.frame(rejection_rate = mean(chisqtest.pval3 < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "3 x 3", test = "ChiSqr"),
      data.frame(rejection_rate = mean(hoeffd.pval3    < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "3 x 3", test = "Hoef"),
      data.frame(rejection_rate = mean(genest.pval3    < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "3 x 3", test = "Genest"),
      data.frame(rejection_rate = mean(brs.pval3       < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "3 x 3", test = "BRS"),
      data.frame(rejection_rate = mean(Btest.pval5     < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "5 x 5", test = "Bn"),
      data.frame(rejection_rate = mean(Ptest.pval5     < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "5 x 5", test = "Pn"),
      data.frame(rejection_rate = mean(chisqtest.pval5 < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "5 x 5", test = "ChiSqr"),
      data.frame(rejection_rate = mean(hoeffd.pval5    < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "5 x 5", test = "Hoef"),
      data.frame(rejection_rate = mean(genest.pval5    < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "5 x 5", test = "Genest"),
      data.frame(rejection_rate = mean(brs.pval5       < 0.05, na.rm = TRUE), n = nk, rho = rj, table = "5 x 5", test = "BRS")
    )
    
  }, error = function(e) {
    message("ERROR n=", params$n, " rho=", params$rho, ": ", conditionMessage(e))
    NULL
  })
}

clusterExport(cl, "run_one")

set.seed(20230101)
clusterSetRNGStream(cl, 20230101)

# -----------------------------------------------------------------------------
# Run parallel computation
# -----------------------------------------------------------------------------
cat("Starting parallel computation on", n_cores, "cores...\n")
t_start <- proc.time()

results_list <- parLapply(cl, param_list, run_one)

t_elapsed <- proc.time() - t_start
cat(sprintf("Done. Wall time: %.1f min\n", t_elapsed["elapsed"] / 60))

stopCluster(cl)

# -----------------------------------------------------------------------------
# Collect and save results
# -----------------------------------------------------------------------------
n_null <- sum(sapply(results_list, is.null))
if (n_null > 0L) warning(n_null, " cells returned NULL and were dropped.")
results_list <- Filter(Negate(is.null), results_list)

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
# Quick check: n = 100 cross-tab
# -----------------------------------------------------------------------------
dt2 <- fread(res_path("rejection_rates_nominal.csv"))

print(
  dt2 |>
    filter(n == 100) |>
    pivot_wider(names_from = test, values_from = rejection_rate) |>
    arrange(table, rho)
)