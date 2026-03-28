# --- 0. Data Generation ----

### load in genetic data----------
load("cleaned_data.Rdata")
library(MASS)
### scale age to 0 - 2
X_SNPs_cleaned[,3] <- 2*(X_SNPs_cleaned[,3] - min(X_SNPs_cleaned[,3]))/(max(X_SNPs_cleaned[,3] - min(X_SNPs_cleaned[,3])))
X_SNPs_cleaned[,2] <- X_SNPs_cleaned[,2] - 1

##### Data generation---------
## data generation for low/high dimensional mimic data
n_sub <- 250
p_dim <- 1500
p <- p_dim + 2
M = 20 # Number of replication
p0 = 20

sim_params <- list(
  n = n_sub,      # Sample size
  p = p,      # Total number of predictors
  p0 = p0,      # Number of true predictors
  rho = 0,    # Correlation for AR-1 covariance matrix
  sigma = NULL,    # Standard deviation of the error term
  k = M,       # Number of simulation replicates
  seed = NULL,  # Seed for reproducibility
  cores = 40    # Number of CPU cores for parallel processing
)

### define the coefficients used to generate Y, original version
beta_mimic_high_dim <- rep(0, p) # beta, including SNPs, age and sex
beta_mimic_high_dim[c(1,2)] = 1 # sex and age, coefficients set as 1
### use SNPs not nearby to each other.
#high_dim_SNPs_index <- (1:10  - 1) * 100 + 3
#beta_mimic_high_dim[high_dim_SNPs_index] = 1 # SNP3, SNP103, SNP203, SNP303,..., SNP903: coefficients set as 1

### New version, first 20
true_beta <- rep(0, p)
true_beta[1:20] <- 1
beta_mimic_high_dim <- true_beta

## Construct Cov matrix (Sigma) for error, AR-1
Sig_mat <- matrix(,nrow = n_sub, ncol = n_sub)
for (i in 1:n_sub) {
  for (j in 1:n_sub) {
    Sig_mat[i,j] <- sim_params$rho^( abs(i - j) )
  }
}

# randomly select subjects and SNPs
set.seed(1234)
mimic_high_dim_selected_SNPs <- sample(4:NCOL(X_SNPs_cleaned), p_dim)
mimic_high_dim_selected_subjects <- sample(1:NROW(X_SNPs_cleaned), n_sub)
X_SNPs_high_dim <- X_SNPs_cleaned[mimic_high_dim_selected_subjects, c(2,3,mimic_high_dim_selected_SNPs)]

#TESTING part
#replace SNPs with random 0,1,2
#X_SNPs_high_dim[,3:(p)] <- sample(c(0,1,2),(p-2)*n_sub,replace = TRUE)
#X_SNPs_high_dim[,3:(p)] <- rnorm((p-2)*n_sub)

## we have an unchanged X and random error terms
## data generation function for parallel computation
Y_mimic_high_dim_mat <- matrix(,nrow = M,ncol = n_sub)
for(i in 1:M){
  err <- mvrnorm(n = 1,mu = rep(0, n_sub),Sigma  = Sig_mat)
  Y_mimic_high_dim_mat[i,] <-  X_SNPs_high_dim %*% beta_mimic_high_dim + 0.6*err
}


# Title: Main Simulation Execution Script
# Description: This script configures and runs the KIQR simulation study.
#              It sources the core package functions.

# --- 1. SETUP ----
# Load all libraries and source all necessary functions
cat("Loading required packages and functions...\n")
library(MASS)
library(mvtnorm)
library(progress)
library(parallel)
library(quantreg)
library(rqPen)
library(dplyr)

# Source the refactored package functions
source("R/loss_penalty_functions.R")
source("R/tuning_criteria.R")
source("R/solvers.R")
source("R/kiqr_api.R")
source("R/utils.R")

# HELPER FUNCTION (only for the baseline rq.pen model)
run_pen_qr_replicate <- function(i, X_input, Y_input, tau, penalty_factor = NULL, pre_selected = NULL) {
  pre_selected <- pre_selected[[i]]
  if(!is.null(pre_selected)){
    if (!is.null(penalty_factor)) {
      #if pre_selected and prior, combine them
      penalty_index <- which(penalty_factor == 0)
      pre_selected_prior <- union(pre_selected, penalty_index)
      penalty_factor_new <- rep(1,length(pre_selected_prior))
      penalty_factor_new[match(penalty_index,pre_selected_prior)] <- 0
      #re-write to penalty_factor
      penalty_factor <- penalty_factor_new
      #browser()
      X_input = X_input[,pre_selected_prior]
    }else{
      X_input = X_input[,pre_selected]
      penalty_factor <- rep(1, NCOL(X_input))
    }
  }else{
    if (is.null(penalty_factor)) {
      penalty_factor <- rep(1, NCOL(X_input))
    }
  }

  fit <- rqPen::rq.pen(
    x = X_input,
    y = as.vector(Y_input[i, ]),
    tau = tau,
    penalty.factor = penalty_factor,
    penalty = "SCAD"
  )
  best_lambda <- select_lambda_hbic(fit)
  as.vector(stats::coef(fit, lambda = best_lambda))
}

run_gwas_replicate <- function(i, X_input, Y_input, tau, gwas_threshold = 5e-8) {
  y_rep        <- Y_input[i, ]
  x_covariates <- X_input[, 1:2]
  beta_gwas    <- rep(0, sim_params$p)

  # Estimate and store base covariate effects (gender and age)
  fit_cov <- quantreg::rq(y_rep ~ x_covariates, tau = tau)
  beta_gwas[1:2] <- coef(fit_cov)[-1] # Get coefficients, excluding the intercept

  # Screen all SNPs via marginal regression
  snp_indices <- 3:sim_params$p
  for (j in snp_indices) {
    X_marginal <- cbind(x_covariates, X_data[, j])
    fit_snp <- quantreg::rq(y_rep ~ X_marginal, tau = tau)
    summary_snp <- summary(fit_snp, se = "nid")
    p_val <- summary_snp$coefficients[4, 4]

    # If significant, store its marginal coefficient; otherwise, it remains zero
    if (p_val < gwas_threshold) {
      beta_gwas[j] <- summary_snp$coefficients[4, 1]
    }
  }
  return(beta_gwas)
}


# --- 2. CONFIGURATION ----

run_config <- list(
  taus = c(0.5, 0.8),
  lambda_grid = seq(0.01, 0.2, length.out = 20),
  zeta_grid = seq(0, 1, length.out = 11),
  high_dim = FALSE
)


# --- 3. DATA GENERATION ----
if(run_config$high_dim){
  # Apply SIS to retain 2000 ------------------------------------------------
  library(SIS)
  SIS_list <- lapply(X = 1:M,function(x){
    SIS_i <- SIS::SIS(x = X_SNPs_high_dim, y = Y_mimic_high_dim_mat[i, ], nsis = 2000)
    SIS_index_i <-union(SIS_i$ix0, c(1,2)) # add 1,2 if not selected.
    SIS_index_i
  })
}else{
  SIS_list <- lapply(X = 1:M,function(x){
  SIS_index_i <- 1:p})
}
X_data <- X_SNPs_high_dim
Y_data <- Y_mimic_high_dim_mat

cat("SKip this step...\n")



# --- 4. DEFINE PRIOR SETS ----
cat("Defining prior information sets...\n")
priors <- list(
  S1 = rep(1, sim_params$p), S2 = rep(1, sim_params$p), S3 = rep(1, sim_params$p),
  S4 = rep(1, sim_params$p)
)
priors$S1[1:20] <- 0                                     # 20 correct variables
priors$S2[c(1:10, (sim_params$p - 1):sim_params$p)] <- 0 # 10 correct, 2 wrong
priors$S3[16:35] <- 0                                    # 5 correct, 15 wrong, now S3, used to be S6
priors$S4[(sim_params$p - 19):sim_params$p] <- 0         # 20 wrong, last 20
#priors$S3[c(1:10, (sim_params$p - 4):sim_params$p)] <- 0 # 10 correct, 5 wrong




# --- 5. RUN SIMULATIONS ----
cat("Starting simulations...\n")
all_results <- list()

for (tau_val in run_config$taus) {
  cat(sprintf("\n--- Running for tau = %.1f ---\n", tau_val))
  tau_id <- paste0("tau_", tau_val)
  all_results[[tau_id]] <- list()

  # 0. GWAS

  cat("Running GWAS QR...\n")
  all_results[[tau_id]]$GWAS_qr <- do.call(
    rbind,
    mclapply(1:sim_params$k, FUN = run_gwas_replicate,
             X_input = X_data, Y_input = Y_data, tau = tau_val,
             mc.cores = sim_params$cores)
  )


  # 1. Traditional Penalized QR
  cat("Running Traditional Penalized QR...\n")
  all_results[[tau_id]]$trad_qr <- do.call(
    rbind,
    mclapply(1:sim_params$k, FUN = run_pen_qr_replicate,
             X_input = X_data, Y_input = Y_data, tau = tau_val,
             pre_selected = SIS_list,
             mc.cores = sim_params$cores)
  )

  # 2. Loop through each prior set
  for (prior_name in names(priors)) {
    cat(sprintf("Running models for prior set: %s\n", prior_name))

    # 2a. Prior-Informed Penalized QR
    all_results[[tau_id]]$prior_qr[[prior_name]] <- do.call(
      rbind,
      mclapply(1:sim_params$k, FUN = run_pen_qr_replicate,
               X_input = X_data, Y_input = Y_data, tau = tau_val,
               penalty_factor = priors[[prior_name]],
               pre_selected = SIS_list,
               mc.cores = sim_params$cores)
    )

    # 2b. Prior-Informed KIQR (using the new API function)
    prior_beta_for_kiqr <- all_results[[tau_id]]$prior_qr[[prior_name]]
    pb_kiqr <- progress_bar$new(total = sim_params$k, format = "Replicate [:bar] :percent")

    kiqr_raw_results <- lapply(1:sim_params$k, function(i) {
      pb_kiqr$tick()
      #if SIS
      penalty_factor = priors[[prior_name]]
      pre_selected = SIS_list[[i]]

      penalty_index <- which(penalty_factor == 0)
      pre_selected_prior <- union(pre_selected, penalty_index)
      penalty_factor_new <- rep(1,length(pre_selected_prior))
      penalty_factor_new[match(penalty_index,pre_selected_prior)] <- 0
      #re-write to penalty_factor
      penalty_factor <- penalty_factor_new
      X_input = X_data[,pre_selected_prior]

      # Call the new, clean API function for the 2D search.
      KIQR_2D_search(
        X = X_input,
        Y = Y_data[i, ],
        beta_p = prior_beta_for_kiqr[i, ],
        tau = tau_val,
        lambda_grid = run_config$lambda_grid,
        zeta_grid = run_config$zeta_grid,
        half_life = 200,
        max_iter = 5000,
        cores = sim_params$cores
      )
    })

    # Unpack results
    all_results[[tau_id]]$kiqr[[prior_name]] <- list(
      betas = do.call(rbind, lapply(kiqr_raw_results, `[[`, "beta")),
      lambdas = sapply(kiqr_raw_results, `[[`, "lambda"),
      zetas = sapply(kiqr_raw_results, `[[`, "zeta")
    )
  }
}

# --- 6. SAVE RESULTS ----
cat("\nSaving results...\n")
save(
  all_results,
  sim_params,
  true_beta,
  file = "mimic_simulation_results_lowdim_test_n250_rho0_M200_error06_gwas.Rdata"
)
cat("Results saved.\n")


