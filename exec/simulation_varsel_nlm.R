library(MASS)
library(mvtnorm)
library(energy)
library(dependence)
library(data.table)
library(parallel)

source("gendataRegMod.R")

# -----------------------------------------------------------------------------
# Directory setup
# -----------------------------------------------------------------------------
dir.create("results",        showWarnings = FALSE)
dir.create("results/varsel", showWarnings = FALSE)
dir.create("plots",          showWarnings = FALSE)
dir.create("plots/varsel",   showWarnings = FALSE)

res_path  <- function(...) file.path("results", "varsel", ...)
plot_path <- function(...) file.path("plots",   "varsel", ...)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
nSim    <- 1000L
nsample <- 100L
p_vals  <- c(1000L, 2500L, 5000L)
rho     <- 0.5

# Rule for p (and q) — consistent with Equation (eq:p_rule) in the paper
p_order <- function(n) max(1L, floor(n^0.3) - 1L)

# Cases 1-4: active set = first 5 covariates; Cases 5-8: first 4
active_cases1_4 <- 1:5
active_cases5_8 <- 1:4

# -----------------------------------------------------------------------------
# Data-generating functions (one per case)
# Cases 1-4: from Liu et al. (2022), Example 1
# Cases 5-8: from Zhu et al. (2011) Ex.3 and Liu et al. (2022) Ex.2
# -----------------------------------------------------------------------------
gendata <- function(case, n, p, rho) {
  sig        <- rho^abs(outer(1:p, 1:p, "-"))
  diag(sig)  <- 1
  
  switch(as.character(case),
         
         # --- Linear / Poisson models (Cases 1-4) ---------------------------------
         "1" = {
           beta <- c(rep(1, 5), rep(0, p - 5))
           X    <- mvrnorm(n, rep(0, p), sig)
           Y    <- as.vector(X %*% beta) + rnorm(n)
           list(X = X, Y = Y)
         },
         "2" = {
           beta <- c(rep(1, 5), rep(0, p - 5))
           U    <- rmvt(n, sigma = diag(p), df = 1)
           X    <- U %*% chol(sig)
           Y    <- as.vector(X %*% beta) + rnorm(n)
           list(X = X, Y = Y)
         },
         "3" = {
           beta <- c(rep(2, 5), rep(0, p - 5))
           X    <- mvrnorm(n, rep(0, p), sig)
           Y    <- as.vector(exp(X %*% beta)) + rnorm(n)
           list(X = X, Y = Y)
         },
         "4" = {
           beta <- c(rep(2, 5), rep(0, p - 5))
           X    <- mvrnorm(n, rep(0, p), sig)
           Y    <- as.numeric(rpois(n, as.vector(exp(X %*% beta))))
           list(X = X, Y = Y)
         },
         
         # --- Nonlinear models (Cases 5-8) ----------------------------------------
         # Case 5: transformation model (Zhu et al. 2011, Ex. 3.a)
         "5" = {
           u    <- runif(4)
           beta <- c(2 - u[1], 2 - u[2], 2 - u[3], 2 - u[4], rep(0, p - 4))
           X    <- mvrnorm(n, rep(0, p), sig)
           Y    <- as.numeric(exp(as.vector(X %*% beta) + rnorm(n)))
           list(X = X, Y = Y)
         },
         # Case 6: multiple-index model (Zhu et al. 2011, Ex. 3.b)
         "6" = {
           u     <- runif(4)
           beta1 <- c(2 - u[1], 2 - u[2], rep(0, p - 2))
           beta2 <- c(0, 0, 2 + u[3], 2 + u[4], rep(0, p - 4))
           X     <- mvrnorm(n, rep(0, p), sig)
           Y     <- as.vector(X %*% beta1) + as.vector(exp(X %*% beta2)) + rnorm(n)
           list(X = X, Y = Y)
         },
         # Case 7: heteroskedastic model (Zhu et al. 2011, Ex. 3.c)
         "7" = {
           u     <- runif(4)
           beta1 <- c(2 - u[1], 2 - u[2], rep(0, p - 2))
           beta2 <- c(0, 0, 2 + u[3], 2 + u[4], rep(0, p - 4))
           X     <- mvrnorm(n, rep(0, p), sig)
           Y     <- as.vector(X %*% beta1) + as.vector(exp(X %*% beta2 + rnorm(n)))
           list(X = X, Y = Y)
         },
         # Case 8: challenging nonlinear structure (Liu et al. 2022, Ex. 2.d)
         "8" = {
           X <- mvrnorm(n, rep(0, p), sig)
           Y <- 1 - 5*(X[,2] + X[,3])^(-3) *
             exp(1 + 10*sin(pi*X[,1]/2) + 5*X[,4]) + rnorm(n)
           list(X = X, Y = Y)
         }
  )
}

# -----------------------------------------------------------------------------
# Cluster setup
# -----------------------------------------------------------------------------
n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
cat("Using", n_cores, "cores\n")
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(MASS); library(mvtnorm); library(energy); library(dependence)
})
clusterExport(cl, c("nSim", "nsample", "rho", "p_order",
                    "active_cases1_4", "active_cases5_8",
                    "gendata", "res_path"))

# -----------------------------------------------------------------------------
# Worker function: one call = one (case, p) combination
# -----------------------------------------------------------------------------
run_one <- function(params) {
  
  case <- params$case
  pval <- params$pval
  
  active <- if (case <= 4L) active_cases1_4 else active_cases5_8
  pord   <- p_order(nsample)
  # q=1 for additive models (Cases 1-4), q=p for nonlinear (Cases 5-8)
  qord   <- if (case <= 4L) 1L else pord
  
  pMat  <- numeric(nSim)
  dMat  <- numeric(nSim)
  BnMat <- numeric(nSim)
  
  set.seed(123 + case * 100L + pval)
  
  for (i in seq_len(nSim)) {
    
    dat <- gendata(case, nsample, pval, rho)
    X0  <- dat$X
    Y0  <- dat$Y
    
    # Pearson correlation screening
    pCorrs       <- abs(cor(X0, Y0, method = "pearson"))
    pCorrOrder   <- order(pCorrs, decreasing = TRUE)
    pMat[i]      <- max(which(pCorrOrder %in% active))
    
    # Distance correlation screening
    dCorrs       <- abs(apply(X0, 2, dcor, y = Y0))
    dCorrOrder   <- order(dCorrs, decreasing = TRUE)
    dMat[i]      <- max(which(dCorrOrder %in% active))
    
    # Bn screening
    BnCorrs      <- apply(X0, 2, function(Xk)
      indeptest(Xk, Y0, basis = "poly",
                p = pord, q = qord)$B_stat)
    BnCorrOrder  <- order(BnCorrs, decreasing = TRUE)
    BnMat[i]     <- max(which(BnCorrOrder %in% active))
  }
  
  out <- data.frame(pMat = pMat, dMat = dMat, BnMat = BnMat)
  fname <- res_path(paste0("case_", case, "_p", pval, ".csv"))
  data.table::fwrite(out, fname)
  cat("Saved:", fname, "\n")
  
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Run in parallel across all (case, p) combinations
# -----------------------------------------------------------------------------
param_list <- Map(
  function(case, pval) list(case = case, pval = pval),
  rep(1:8, each = length(p_vals)),
  rep(p_vals, times = 8L)
)

clusterSetRNGStream(cl, 20230101)

cat("Starting feature screening simulation...\n")
t_start <- proc.time()

parLapply(cl, param_list, run_one)

t_elapsed <- proc.time() - t_start
cat(sprintf("Done. Wall time: %.1f min\n", t_elapsed["elapsed"] / 60))

stopCluster(cl)

# -----------------------------------------------------------------------------
# Summary function
# -----------------------------------------------------------------------------
out_table_screen <- function(mat, d1 = 20, d2 = 40) {
  R <- t(apply(mat, 2, quantile,
               probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
  colnames(R) <- c("5%", "25%", "50%", "75%", "95%")
  
  S1 <- colMeans(mat <= d1)
  S2 <- colMeans(mat <= d2)
  
  list(R = R, S1 = S1, S2 = S2)
}

# -----------------------------------------------------------------------------
# Collect and summarise results
# -----------------------------------------------------------------------------
cases  <- 1:8
p_labs <- c("1000", "2500", "5000")

results <- list()

for (case in cases) {
  results[[case]] <- list()
  for (pval in p_vals) {
    mat <- as.matrix(data.table::fread(
      res_path(paste0("case_", case, "_p", pval, ".csv"))
    ))
    results[[case]][[as.character(pval)]] <- out_table_screen(mat)
  }
}

# Save summary as CSV for each case
for (case in cases) {
  R_all  <- do.call(rbind, lapply(p_vals, function(pv)
    results[[case]][[as.character(pv)]]$R))
  S1_all <- do.call(rbind, lapply(p_vals, function(pv)
    t(results[[case]][[as.character(pv)]]$S1)))
  S2_all <- do.call(rbind, lapply(p_vals, function(pv)
    t(results[[case]][[as.character(pv)]]$S2)))
  
  rownames(R_all) <- rownames(S1_all) <- rownames(S2_all) <- p_labs
  
  data.table::fwrite(as.data.frame(R_all),
                     res_path(paste0("summary_R_case", case, ".csv")),
                     row.names = TRUE)
  data.table::fwrite(as.data.frame(S1_all),
                     res_path(paste0("summary_S1_case", case, ".csv")),
                     row.names = TRUE)
  data.table::fwrite(as.data.frame(S2_all),
                     res_path(paste0("summary_S2_case", case, ".csv")),
                     row.names = TRUE)
}

cat("Summary files written to results/varsel/\n")