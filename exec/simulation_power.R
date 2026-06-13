# =============================================================================
# Simulation study: power comparison of independence tests
#
# Tests compared:
#   HHG   – Heller-Heller-Gorfine
#   dCov  – Distance covariance (energy)
#   MIC   – Maximal Information Coefficient (minerva)
#   Hoef  – Hoeffding's D (Hmisc)
#   Bn    – B-statistic, polynomial basis (dependence)
#   Pn    – P-statistic, polynomial basis (dependence)
#   Genest– Cramér-von Mises copula-based test, Tn2 (MixedIndTests)
#
# Design:
#   n     = 100 observations per sample
#   nsim  = 10,000 replications per (type, noise) cell
#   types = 10 functional relationships
#   noise = seq(0.1, 2, 0.1)  →  200 cells total
#
# NOTE on Genest timing:
#   TestIndCopula uses B = 1000 Gaussian multipliers internally.
#   With nsim = 10,000 this is very slow; consider nsim = 100–500 for
#   a quick pilot run (change the `nsim` constant below).
#
# Changes vs previous version:
#   - Data passed in-memory via param_list (no CSV I/O on workers)
#   - tryCatch around entire run_one body + around TestIndCopula calls
#   - na.rm = TRUE in Genest power calculation
# =============================================================================

library(data.table)
library(parallel)
library(HHG)
library(energy)
library(minerva)
library(Hmisc)
library(dependence)
library(MixedIndTests)
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Global constants
# -----------------------------------------------------------------------------
n     <- 100
nsim  <- 1e4          # set to e.g. 200 for a pilot run

types      <- 1:10
noises     <- seq(0.1, 2, 0.1)
type_names <- c("Linear", "Quadratic", "Cubic", "Sine",
                "Fourth root", "Circle", "Two curves", "X", "Diamond", "XsdY")

# -----------------------------------------------------------------------------
# Directory helpers
# -----------------------------------------------------------------------------
dir_results <- "results_power"
dir_plots   <- "plots_power"
dir.create(dir_results, showWarnings = FALSE)
dir.create(dir_plots,   showWarnings = FALSE)

res_path  <- function(...) file.path(dir_results, ...)
plot_path <- function(...) file.path(dir_plots,   ...)

# -----------------------------------------------------------------------------
# Data-generating functions
# -----------------------------------------------------------------------------
datagen.noise <- function(n, typ, noise = 0.1) {
  switch(typ,
         linear(n, noise),
         quadratic(n, noise),
         cubic(n, noise),
         sine(n, noise),
         x14(n, noise),
         circle(n, noise),
         twocurves(n, noise),
         Xfun(n, noise),
         Diamond(n, noise),
         XsdY(n, noise)
  )
}

linear <- function(n, noise) {
  x <- runif(n); y <- x + noise * rnorm(n)
  rbind(x, y)
}
quadratic <- function(n, noise) {
  x <- runif(n); y <- 4 * (x - 0.5)^2 + noise * rnorm(n)
  rbind(x, y)
}
cubic <- function(n, noise) {
  x <- runif(n)
  y <- 128*(x - 1/3)^3 - 48*(x - 1/3)^2 - 12*(x - 1/3) + noise * rnorm(n, sd = 3)
  rbind(x, y)
}
sine <- function(n, noise) {
  x <- runif(n); y <- sin(4*pi*x) + noise * rnorm(n, sd = 2)
  rbind(x, y)
}
x14 <- function(n, noise) {
  x <- runif(n); y <- x^(1/4) + noise * rnorm(n, sd = 0.5)
  rbind(x, y)
}
circle <- function(n, noise) {
  x <- runif(n); r <- rbinom(n, 1, 0.5)
  y <- (2*r - 1) * sqrt(1 - (2*x - 1)^2) + noise * rnorm(n, sd = 0.5)
  rbind(x, y)
}
twocurves <- function(n, noise) {
  x <- runif(n); r <- rbinom(n, 1, 0.5)
  y <- 2*r*x + (1 - r)*sqrt(x)/2 + noise * rnorm(n)
  rbind(x, y)
}
Xfun <- function(n, noise) {
  x <- runif(n); r <- rbinom(n, 1, 0.5)
  y <- r*x + (1 - r)*(1 - x) + noise * rnorm(n, sd = 1/5)
  rbind(x, y)
}
Diamond <- function(n, noise) {
  x  <- runif(n)
  r1 <- runif(n, 0.5 - x, 0.5 + x)
  r2 <- runif(n, x - 0.5, 1.5 - x)
  y  <- r1*(x < 0.5) + r2*(x >= 0.5) + noise * rnorm(n, sd = 1/10)
  rbind(x, y)
}
XsdY <- function(n, noise) {
  x <- sqrt(rchisq(n, 1)); y <- rnorm(n)*x + rnorm(n)*noise
  rbind(x, y)
}

# -----------------------------------------------------------------------------
# Pre-generate all simulation data IN MEMORY and embed into param_list
# Workers receive data directly via the param payload — no disk I/O needed.
# -----------------------------------------------------------------------------
cat("Generating simulation data in memory...\n")
set.seed(202210)

grid       <- expand.grid(ty = types, no = noises)
param_list <- vector("list", nrow(grid))

for (k in seq_len(nrow(grid))) {
  ty <- grid$ty[k]
  no <- grid$no[k]
  sims        <- matrix(nrow = n * nsim, ncol = 3,
                        dimnames = list(NULL, c("x", "y", "yperm")))
  sims[, 1:2] <- t(datagen.noise(n * nsim, ty, no))
  sims[, 3]   <- sims[sample.int(nrow(sims)), 2]   # global permutation → H0
  param_list[[k]] <- list(ty = ty, no = no, sims = sims)
}

cat("Data generation complete.\n")

# -----------------------------------------------------------------------------
# Parallel cluster setup
# -----------------------------------------------------------------------------
n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
cat("Using", n_cores, "cores\n")
cl <- makeCluster(n_cores, rscript_args = "--no-init-file")

lib_path <- .libPaths()   # prende i path dal master
clusterExport(cl, "lib_path")
clusterEvalQ(cl, .libPaths(lib_path))

clusterEvalQ(cl, {
  library(HHG)
  library(energy)
  library(minerva)
  library(Hmisc)
  library(dependence)
  library(MixedIndTests)
})

clusterExport(cl, c("n", "nsim", "type_names"))

# -----------------------------------------------------------------------------
# Worker function: processes one (type, noise) cell
# Data arrives in params$sims — no file reading required.
# The entire body is wrapped in tryCatch so a worker never crashes silently.
# -----------------------------------------------------------------------------
run_one <- function(params) {
  tryCatch({
    
    ty   <- params$ty
    no   <- params$no
    sims <- params$sims   # matrix (n*nsim) x 3: x, y, yperm
    
    # Polynomial orders for the Bn / Pn test
    p_ord <- max(1L, floor(100^0.3) - 1L)   # = 3 for n = 100
    q_ord <- p_ord
    
    # ------- Storage vectors ------------------------------------------------
    hhg_dep  <- numeric(nsim); hhg_ind  <- numeric(nsim); hhg_tim  <- numeric(nsim)
    dcov_dep <- numeric(nsim); dcov_ind <- numeric(nsim); dcov_tim <- numeric(nsim)
    mic_dep  <- numeric(nsim); mic_ind  <- numeric(nsim); mic_tim  <- numeric(nsim)
    hoef_dep <- numeric(nsim); hoef_ind <- numeric(nsim); hoef_tim <- numeric(nsim)
    ourB_dep <- numeric(nsim); ourB_ind <- numeric(nsim)
    ourP_dep <- numeric(nsim); ourP_ind <- numeric(nsim); our_tim  <- numeric(nsim)
    genest_dep <- numeric(nsim); genest_ind <- numeric(nsim)
    genest_tim <- numeric(nsim)
    
    # ------- Main simulation loop -------------------------------------------
    for (i in seq_len(nsim)) {
      
      X <- sims[((i - 1L)*n + 1L):(i*n), ]   # rows: obs; cols: x, y, yperm
      
      # ---- HHG --------------------------------------------------------------
      t0 <- proc.time()["elapsed"]
      Dx  <- as.matrix(dist(X[, 1], diag = TRUE, upper = TRUE))
      Dy  <- as.matrix(dist(X[, 2], diag = TRUE, upper = TRUE))
      Dy0 <- as.matrix(dist(X[, 3], diag = TRUE, upper = TRUE))
      hhg_dep[i] <- hhg.test(Dx, Dy,  nr.perm = 0)[[1]]
      hhg_ind[i] <- hhg.test(Dx, Dy0, nr.perm = 0)[[1]]
      hhg_tim[i] <- proc.time()["elapsed"] - t0
      
      # ---- dCov -------------------------------------------------------------
      t0 <- proc.time()["elapsed"]
      dcov_dep[i] <- dcov.test(X[, 1], X[, 2])$statistic
      dcov_ind[i] <- dcov.test(X[, 1], X[, 3])$statistic
      dcov_tim[i] <- proc.time()["elapsed"] - t0
      
      # ---- MIC --------------------------------------------------------------
      t0 <- proc.time()["elapsed"]
      mic_dep[i] <- mine_stat(X[, 1], X[, 2])
      mic_ind[i] <- mine_stat(X[, 1], X[, 3])
      mic_tim[i] <- proc.time()["elapsed"] - t0
      
      # ---- Hoeffding's D ----------------------------------------------------
      t0 <- proc.time()["elapsed"]
      hoef_dep[i] <- hoeffd(X[, 1], X[, 2])$D[2]
      hoef_ind[i] <- hoeffd(X[, 1], X[, 3])$D[2]
      hoef_tim[i] <- proc.time()["elapsed"] - t0
      
      # ---- Bn / Pn (dependence package) ------------------------------------
      t0 <- proc.time()["elapsed"]
      res_dep  <- indeptest(X[, 1], X[, 2], basis = "poly", p = p_ord, q = q_ord)
      res_ind  <- indeptest(X[, 1], X[, 3], basis = "poly", p = p_ord, q = q_ord)
      ourB_dep[i] <- res_dep$B_stat;  ourB_ind[i] <- res_ind$B_stat
      ourP_dep[i] <- res_dep$P_stat;  ourP_ind[i] <- res_ind$P_stat
      our_tim[i]  <- proc.time()["elapsed"] - t0
      
      # ---- Genest (MixedIndTests) ------------------------------------------
      t0 <- proc.time()["elapsed"]
      genest_dep[i] <- tryCatch({
        res <- TestIndCopula(cbind(X[, 1], X[, 2]),
                             trunc.level = 2, B = 1000,
                             par = FALSE, graph = FALSE)
        res$pvalue$Tn2 / 100
      }, error = function(e) NA_real_)
      
      genest_ind[i] <- tryCatch({
        res <- TestIndCopula(cbind(X[, 1], X[, 3]),
                             trunc.level = 2, B = 1000,
                             par = FALSE, graph = FALSE)
        res$pvalue$Tn2 / 100
      }, error = function(e) NA_real_)
      genest_tim[i] <- proc.time()["elapsed"] - t0
      
    }  # end simulation loop
    
    # ------- Power estimates ------------------------------------------------
    data.frame(
      power = c(
        mean(hhg_dep  > quantile(hhg_ind,  0.95)),
        mean(dcov_dep > quantile(dcov_ind, 0.95)),
        mean(mic_dep  > quantile(mic_ind,  0.95)),
        mean(hoef_dep > quantile(hoef_ind, 0.95)),
        mean(ourB_dep > quantile(ourB_ind, 0.95)),
        mean(ourP_dep > quantile(ourP_ind, 0.95)),
        mean(genest_dep < 0.05, na.rm = TRUE)
      ),
      time = c(
        mean(hhg_tim),
        mean(dcov_tim),
        mean(mic_tim),
        mean(hoef_tim),
        mean(our_tim),
        mean(our_tim),   # Bn and Pn share the same timing call
        mean(genest_tim)
      ),
      noise = no,
      type  = type_names[ty],
      test  = c("HHG", "dCov", "MIC", "Hoef", "Bn", "Pn", "Genest")
    )
    
  }, error = function(e) {
    message("ERROR type=", params$ty, " noise=", params$no, ": ", conditionMessage(e))
    NULL
  })
}

# -----------------------------------------------------------------------------
# Export worker function to cluster
# -----------------------------------------------------------------------------
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

# Drop any NULL entries from cells that errored
n_null <- sum(sapply(results_list, is.null))
if (n_null > 0L) warning(n_null, " cells returned NULL and were dropped.")
results_list <- Filter(Negate(is.null), results_list)

tests <- rbindlist(results_list)
tests$type <- factor(tests$type, levels = type_names)
fwrite(tests, res_path("test_power.csv"))
cat("Results saved to", res_path("test_power.csv"), "\n")

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
dt <- fread(res_path("test_power.csv"))
dt$type <- factor(dt$type, levels = type_names)

# All tests
p_all <- dt |>
  ggplot(aes(x = noise, y = power,
             color = test, linetype = test, shape = test)) +
  geom_line() +
  geom_point() +
  facet_wrap(~type, ncol = 2) +
  labs(color = NULL, linetype = NULL, shape = NULL,
       y = "Power", x = "Noise") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(plot_path("power_all_tests.pdf"), plot = p_all, width = 7, height = 9)

# Bn vs competitors (Pn omitted — near-identical power to Bn)
p_no_pn <- dt |>
  filter(test != "Pn") |>
  ggplot(aes(x = noise, y = power,
             color = test, linetype = test, shape = test)) +
  geom_line() +
  geom_point() +
  facet_wrap(~type, ncol = 2) +
  labs(color = NULL, linetype = NULL, shape = NULL,
       y = "Power", x = "Noise") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(plot_path("power_Bn_vs_competitors.pdf"), plot = p_no_pn, width = 7, height = 9)

cat("Plots saved to", dir_plots, "\n")

# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------

# Mean power and SD across all (type, noise) cells
mean_power <- dt |>
  group_by(test) |>
  summarise(mean_power = mean(power),
            sd_power   = sd(power),
            .groups    = "drop") |>
  arrange(desc(mean_power))

cat("\n--- Mean power by test ---\n")
print(mean_power)
fwrite(mean_power, res_path("mean_power.csv"))

# Mean rank across cells (lower = better)
wdt       <- dt |> select(-time) |> pivot_wider(names_from = test, values_from = power)
test_cols <- c("HHG", "dCov", "MIC", "Hoef", "Bn", "Pn", "Genest")
rdt       <- t(apply(
  as.matrix(wdt[, test_cols]), 1,
  data.table::frankv, order = -1L, ties.method = "min"
))
colnames(rdt) <- paste0("rank_", test_cols)
mean_ranks    <- sort(colMeans(rdt))

cat("\n--- Mean rank by test (lower = better) ---\n")
print(mean_ranks)
fwrite(as.data.frame(t(mean_ranks)), res_path("mean_ranks.csv"))

# Mean computation time per replication
mean_time <- dt |>
  group_by(test) |>
  summarise(mean_time_sec = mean(time), .groups = "drop") |>
  arrange(mean_time_sec)

cat("\n--- Mean time per replication (seconds) ---\n")
print(mean_time)
fwrite(mean_time, res_path("mean_time.csv"))