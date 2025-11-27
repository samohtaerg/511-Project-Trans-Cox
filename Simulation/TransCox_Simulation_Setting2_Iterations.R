################################################################################
# TransCox Simulation Setting 2 
################################################################################

# Load required packages
library(survival)

# Create output directory
output_dir <- "simulation_results_setting2"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat(" Created directory:", output_dir, "\n\n")
} else {
  cat(" Using existing directory:", output_dir, "\n\n")
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

GenSimData <- function(nprim = 200, naux = 500, setting = 1) {
  X1_p <- runif(nprim, 0, 1)
  X2_p <- rbinom(nprim, size = 1, prob = 0.5)
  X1_a <- runif(naux, 0, 1)
  X2_a <- rbinom(naux, size = 1, prob = 0.5)

  b1_p <- -0.5
  b2_p <- 0.5
  XB_p <- cbind(X1_p, X2_p) %*% c(b1_p, b2_p)
  scale_p <- exp(-XB_p/2)
  shape_p <- 2

  if(setting %in% c(1, 3)) {
    b1_a <- -0.5
    b2_a <- 0.5
  } else {
    b1_a <- -0.5
    b2_a <- 0.2
  }

  if(setting %in% c(1, 2)) {
    XB_a <- cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
    scale_a <- exp(-XB_a/2)
    shape_a <- 2
  } else {
    XB_a <- cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
    scale_a <- exp(-XB_a/2) * 3/2
    shape_a <- 2
  }

  T_p <- rweibull(nprim, shape = shape_p, scale = scale_p)
  C_p <- runif(nprim, 0, 4.55)
  survTime_p <- ifelse(T_p < C_p, T_p, C_p)
  event_p <- ifelse(T_p < C_p, 2, 1)

  T_a <- rweibull(naux, shape = shape_a, scale = scale_a)
  C_a <- runif(naux, 0, 4.55)
  survTime_a <- ifelse(T_a < C_a, T_a, C_a)
  event_a <- ifelse(T_a < C_a, 2, 1)

  primData <- data.frame(X1 = X1_p, X2 = X2_p, time = survTime_p, status = event_p)
  auxData <- data.frame(X1 = X1_a, X2 = X2_a, time = survTime_a, status = event_a)

  return(list(primData = primData, auxData = auxData, true_beta = c(b1_p, b2_p)))
}

GetAuxSurv <- function(auxData, weights = NULL, cov = c("X1", "X2")) {
  functext <- paste0("res.cox <- coxph(Surv(time, status == 2) ~ ",
                     paste(cov, collapse = "+"),
                     ", data = auxData, weights = weights)")
  res.cox <- NULL
  eval(parse(text = functext))
  bhest <- basehaz(res.cox, centered = FALSE)
  estR <- res.cox$coefficients
  q <- data.frame(cumHazards = bhest$hazard, breakPoints = bhest$time)
  return(list(estR = estR, q = q))
}

deltaQ <- function(primData, q) {
  primData <- primData[order(primData$time), ]
  event_times <- primData$time[primData$status == 2]
  dQ_vec <- rep(NA, length(event_times))
  cumQ_vec <- rep(NA, length(event_times))

  for(i in 1:length(event_times)) {
    t_i <- event_times[i]
    if(t_i <= min(q$breakPoints)) {
      cumQ_i <- 0
    } else if(t_i >= max(q$breakPoints)) {
      cumQ_i <- max(q$cumHazards)
    } else {
      idx <- which.min(abs(q$breakPoints - t_i))
      if(q$breakPoints[idx] > t_i && idx > 1) idx <- idx - 1
      cumQ_i <- q$cumHazards[idx]
    }
    cumQ_vec[i] <- cumQ_i
    if(i == 1) {
      dQ_vec[i] <- cumQ_i
    } else {
      dQ_vec[i] <- cumQ_i - cumQ_vec[i-1]
    }
  }
  return(data.frame(t = event_times, dQ = dQ_vec, cumQ_upd = cumQ_vec))
}

GetPrimaryParam <- function(primData, q, estR) {
  primData <- primData[order(primData$time), ]
  dQ <- deltaQ(primData, q)
  fullcum <- rep(0, nrow(primData))
  idx0 <- match(dQ$t, primData$time)
  fullcum[idx0] <- dQ$cumQ_upd
  for(i in 1:length(fullcum)) {
    if(fullcum[i] == 0 & i > 1) fullcum[i] <- fullcum[i-1]
  }
  primData$fullCumQ <- fullcum

  Ximat <- matrix(0, nrow(primData), nrow(primData))
  for(i in 1:nrow(primData)) {
    tmpidx <- rep(0, nrow(primData))
    for(j in 1:i) {
      if(primData$status[j] == 2) tmpidx[j] <- 1
    }
    Ximat[i,] <- tmpidx
  }
  Xinn <- Ximat[, which(primData$status == 2)]
  return(list(primData = primData, Xinn = Xinn, dQ = dQ, estR = estR))
}

TransCox_R <- function(CovData, cumH, hazards, status, estR, Xinn,
                       lambda1, lambda2, maxit = 1000) {
  n_params <- length(estR)
  n_hazards <- length(hazards)
  n_total <- n_params + n_hazards

  CovData <- as.matrix(CovData)
  Xinn <- as.matrix(Xinn)
  event_idx <- which(status == 2)

  objective <- function(params) {
    eta <- params[1:n_params]
    xi <- params[(n_params + 1):n_total]
    beta_adj <- estR + eta
    haz_adj <- hazards + xi

    if(any(haz_adj <= 0)) return(1e10)

    XB <- as.vector(CovData %*% beta_adj)
    cum_haz_adj <- cumH + as.vector(Xinn %*% xi)

    ll_event <- sum(XB[event_idx])
    ll_baseline <- sum(log(haz_adj))
    ll_risk <- -sum(cum_haz_adj * exp(XB))
    log_lik <- ll_event + ll_baseline + ll_risk
    penalty <- lambda1 * sum(abs(eta)) + lambda2 * sum(abs(xi))

    return(-log_lik + penalty)
  }

  init_params <- rep(0, n_total)
  result <- optim(par = init_params, fn = objective, method = "L-BFGS-B",
                  control = list(maxit = maxit))

  eta <- result$par[1:n_params]
  xi <- result$par[(n_params + 1):n_total]

  return(list(eta = eta, xi = xi, convergence = result$convergence))
}

runTransCox_one <- function(Pout, l1 = 1, l2 = 1, maxit = 1000, cov = c('X1', 'X2')) {
  CovData <- Pout$primData[, cov]
  status <- Pout$primData[, "status"]
  cumH <- Pout$primData$fullCumQ
  hazards <- Pout$dQ$dQ

  result <- TransCox_R(CovData = as.matrix(CovData), cumH = cumH,
                       hazards = hazards, status = status,
                       estR = Pout$estR, Xinn = Pout$Xinn,
                       lambda1 = l1, lambda2 = l2, maxit = maxit)

  return(list(eta = result$eta, xi = result$xi,
              new_beta = Pout$estR + result$eta,
              new_IntH = Pout$dQ$dQ + result$xi,
              time = Pout$primData[status == 2, "time"],
              convergence = result$convergence))
}

# =============================================================================
# RUN 100 ITERATIONS
# =============================================================================

cat("Running 100 iterations of Setting 1 simulation...\n")
cat("Parameters: nprim=200, naux=500, lambda1=0.5, lambda2=0.5\n\n")

N_ITERATIONS <- 100
SETTING <- 2
N_PRIMARY <- 200
N_AUXILIARY <- 500
LAMBDA1 <- 0.5
LAMBDA2 <- 0.5

# Storage for results
results_list <- list()
results_summary <- data.frame(
  iteration = integer(),
  true_beta1 = numeric(),
  true_beta2 = numeric(),
  aux_beta1 = numeric(),
  aux_beta2 = numeric(),
  cox_beta1 = numeric(),
  cox_beta2 = numeric(),
  pooled_beta1 = numeric(),
  pooled_beta2 = numeric(),
  transcox_beta1 = numeric(),
  transcox_beta2 = numeric(),
  aux_mse = numeric(),
  cox_mse = numeric(),
  pooled_mse = numeric(),
  transcox_mse = numeric(),
  n_events_prim = integer(),
  n_events_aux = integer(),
  convergence = integer()
)

# Run iterations
set.seed(2025)  # For reproducibility
pb <- txtProgressBar(min = 0, max = N_ITERATIONS, style = 3)

for(iter in 1:N_ITERATIONS) {
  setTxtProgressBar(pb, iter)

  # Generate data
  sim_data <- GenSimData(nprim = N_PRIMARY, naux = N_AUXILIARY, setting = SETTING)
  pData <- sim_data$primData
  aData <- sim_data$auxData
  true_beta <- sim_data$true_beta

  # Fit auxiliary model
  Cout <- GetAuxSurv(aData, cov = c("X1", "X2"))

  # Prepare primary data
  Pout <- GetPrimaryParam(pData, q = Cout$q, estR = Cout$estR)

  # Run TransCox
  Tres <- runTransCox_one(Pout, l1 = LAMBDA1, l2 = LAMBDA2,
                          maxit = 1000, cov = c("X1", "X2"))

  # Standard Cox
  standard_cox <- coxph(Surv(time, status == 2) ~ X1 + X2, data = pData)
  
  # Pooled Cox (combine primary + auxiliary)
  pooledData <- rbind(pData, aData)
  pooled_cox <- coxph(Surv(time, status == 2) ~ X1 + X2, data = pooledData)

  # Calculate MSE
  aux_mse <- sum((Cout$estR - true_beta)^2)
  cox_mse <- sum((coef(standard_cox) - true_beta)^2)
  pooled_mse <- sum((coef(pooled_cox) - true_beta)^2)
  transcox_mse <- sum((Tres$new_beta - true_beta)^2)

  # Store results
  results_summary <- rbind(results_summary, data.frame(
    iteration = iter,
    true_beta1 = true_beta[1],
    true_beta2 = true_beta[2],
    aux_beta1 = Cout$estR[1],
    aux_beta2 = Cout$estR[2],
    cox_beta1 = coef(standard_cox)[1],
    cox_beta2 = coef(standard_cox)[2],
    pooled_beta1 = coef(pooled_cox)[1],
    pooled_beta2 = coef(pooled_cox)[2],
    transcox_beta1 = Tres$new_beta[1],
    transcox_beta2 = Tres$new_beta[2],
    aux_mse = aux_mse,
    cox_mse = cox_mse,
    pooled_mse = pooled_mse,
    transcox_mse = transcox_mse,
    n_events_prim = sum(pData$status == 2),
    n_events_aux = sum(aData$status == 2),
    convergence = Tres$convergence
  ))

  # Store detailed results
  results_list[[iter]] <- list(
    data = sim_data,
    aux_fit = Cout,
    transcox = Tres,
    standard_cox = standard_cox
  )
}

close(pb)
cat("\n All iterations completed\n\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save summary table
write.csv(results_summary,
          file = file.path(output_dir, "results_summary.csv"),
          row.names = FALSE)
cat(" Saved results_summary.csv\n")

# Save detailed results
saveRDS(results_list, file = file.path(output_dir, "results_detailed.rds"))
cat(" Saved results_detailed.rds\n\n")