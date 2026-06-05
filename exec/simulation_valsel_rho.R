library(MASS)
library(mvtnorm)
library(energy)
library(dependence)
library(data.table)
library(parallel)

# -----------------------------------------------------------------------------
# Fixed setup
# -----------------------------------------------------------------------------
nSim    <- 1000L
nsample <- 100L
p_vals  <- c(1000L, 2500L, 5000L)
rho_vals <- c(0.2, 0.5, 0.8)

p_order <- function(n) max(1L, floor(n^0.3) - 1L)

active_cases1_4 <- 1:5
active_cases5_8 <- 1:4

# -----------------------------------------------------------------------------
# Data-generating functions
# -----------------------------------------------------------------------------
gendata <- function(case, n, p, rho) {
  sig       <- rho^abs(outer(1:p, 1:p, "-"))
  diag(sig) <- 1
  
  switch(as.character(case),
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
         "5" = {
           u    <- runif(4)
           beta <- c(2 - u[1], 2 - u[2], 2 - u[3], 2 - u[4], rep(0, p - 4))
           X    <- mvrnorm(n, rep(0, p), sig)
           Y    <- as.vector(exp(as.vector(X %*% beta) + rnorm(n)))
           list(X = X, Y = Y)
         },
         "6" = {
           u     <- runif(4)
           beta1 <- c(2 - u[1], 2 - u[2], rep(0, p - 2))
           beta2 <- c(0, 0, 2 + u[3], 2 + u[4], rep(0, p - 4))
           X     <- mvrnorm(n, rep(0, p), sig)
           Y     <- as.vector(X %*% beta1) + as.vector(exp(X %*% beta2)) + rnorm(n)
           list(X = X, Y = Y)
         },
         "7" = {
           u     <- runif(4)
           beta1 <- c(2 - u[1], 2 - u[2], rep(0, p - 2))
           beta2 <- c(0, 0, 2 + u[3], 2 + u[4], rep(0, p - 4))
           X     <- mvrnorm(n, rep(0, p), sig)
           Y     <- as.vector(X %*% beta1) + as.vector(exp(X %*% beta2 + rnorm(n)))
           list(X = X, Y = Y)
         },
         "8" = {
           X <- mvrnorm(n, rep(0, p), sig)
           Y <- 1 - 5*(X[,2] + X[,3])^(-3) *
             exp(1 + 10*sin(pi*X[,1]/2) + 5*X[,4]) + rnorm(n)
           list(X = X, Y = Y)
         }
  )
}

# -----------------------------------------------------------------------------
# Worker function
# -----------------------------------------------------------------------------
run_one <- function(params) {
  case     <- params$case
  pval     <- params$pval
  rho      <- params$rho
  res_path <- params$res_path
  
  active <- if (case <= 4L) active_cases1_4 else active_cases5_8
  pord   <- p_order(nsample)
  qord   <- if (case <= 4L) 1L else pord
  
  pMat  <- numeric(nSim)
  dMat  <- numeric(nSim)
  BnMat <- numeric(nSim)
  
  set.seed(123 + case * 100L + pval)
  
  for (i in seq_len(nSim)) {
    dat <- gendata(case, nsample, pval, rho)
    X0  <- dat$X
    Y0  <- dat$Y
    
    pCorrs      <- abs(cor(X0, Y0, method = "pearson"))
    pCorrOrder  <- order(pCorrs, decreasing = TRUE)
    pMat[i]     <- max(which(pCorrOrder %in% active))
    
    dCorrs      <- abs(apply(X0, 2, dcor, y = Y0))
    dCorrOrder  <- order(dCorrs, decreasing = TRUE)
    dMat[i]     <- max(which(dCorrOrder %in% active))
    
    BnCorrs     <- apply(X0, 2, function(Xk)
      indeptest(Xk, Y0, basis = "poly", p = pord, q = qord)$B_stat)
    BnCorrOrder <- order(BnCorrs, decreasing = TRUE)
    BnMat[i]    <- max(which(BnCorrOrder %in% active))
  }
  
  out   <- data.frame(pMat = pMat, dMat = dMat, BnMat = BnMat)
  fname <- file.path(res_path, paste0("case_", case, "_p", pval, ".csv"))
  data.table::fwrite(out, fname)
  cat("Saved:", fname, "\n")
  
  invisible(NULL)
}

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
# Cluster setup — opened once, reused across rho values
# -----------------------------------------------------------------------------
n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
cat("Using", n_cores, "cores\n")
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(MASS); library(mvtnorm); library(energy); library(dependence)
})
clusterExport(cl, c("nSim", "nsample", "p_order",
                    "active_cases1_4", "active_cases5_8", "gendata"))

clusterSetRNGStream(cl, 20230101)

# -----------------------------------------------------------------------------
# Loop over rho values
# -----------------------------------------------------------------------------
t_total <- proc.time()

for (rho in rho_vals) {
  
  rho_tag  <- paste0("rho", gsub("\\.", "", sprintf("%.1f", rho)))
  res_base <- file.path("results", "varsel", rho_tag)
  p_labs   <- as.character(p_vals) 
  dir.create(res_base, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("plots", "varsel", rho_tag),
             recursive = TRUE, showWarnings = FALSE)
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("rho =", rho, "| saving to:", res_base, "\n")
  cat(strrep("=", 60), "\n", sep = "")
  
  # Build param list — res_path and rho passed explicitly to each worker
  param_list <- Map(
    function(case, pval) list(case    = case,
                              pval    = pval,
                              rho     = rho,
                              res_path = res_base),
    rep(1:8, each = length(p_vals)),
    rep(p_vals, times = 8L)
  )
  
  t_start <- proc.time()
  parLapply(cl, param_list, run_one)
  t_elapsed <- proc.time() - t_start
  cat(sprintf("rho = %.1f done. Wall time: %.1f min\n",
              rho, t_elapsed["elapsed"] / 60))
  
  # ---- Summarise results ---------------------------------------------------
  cases  <- 1:8
  p_labs <- as.character(p_vals)
  results <- list()
  
  for (case in cases) {
    results[[case]] <- list()
    for (pval in p_vals) {
      mat <- as.matrix(data.table::fread(
        file.path(res_base, paste0("case_", case, "_p", pval, ".csv"))
      ))
      results[[case]][[as.character(pval)]] <- out_table_screen(mat)
    }
  }
  
  for (case in cases) {
    R_all  <- do.call(rbind, lapply(p_vals, function(pv)
      results[[case]][[as.character(pv)]]$R))
    S1_all <- do.call(rbind, lapply(p_vals, function(pv)
      results[[case]][[as.character(pv)]]$S1))
    S2_all <- do.call(rbind, lapply(p_vals, function(pv)
      results[[case]][[as.character(pv)]]$S2))
    
    # R_all: nrow = length(p_vals) * 3 methods — rownames auto da rbind
    # S1/S2: nrow = length(p_vals)
    rownames(S1_all) <- rownames(S2_all) <- p_labs
    
    data.table::fwrite(as.data.frame(R_all),
                       file.path(res_base, paste0("summary_R_case", case, ".csv")),
                       row.names = TRUE)
    data.table::fwrite(as.data.frame(S1_all),
                       file.path(res_base, paste0("summary_S1_case", case, ".csv")),
                       row.names = TRUE)
    data.table::fwrite(as.data.frame(S2_all),
                       file.path(res_base, paste0("summary_S2_case", case, ".csv")),
                       row.names = TRUE)
  }
  cat("Summary CSVs written to", res_base, "\n")
}

stopCluster(cl)

t_total_elapsed <- proc.time() - t_total
cat(sprintf("\nAll rho values done. Total wall time: %.1f min\n",
            t_total_elapsed["elapsed"] / 60))