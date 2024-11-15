# -------------------------------- TEAM CONTRIBUTIONS ----------------------------------#
#
# - Xu Guo: S2743905
#           50% - Led the development of the linear mixed model structure 
#                  and optimized the fitting process using optimization techniques.
#                        
# - Xie Ying: S2510813
#           35% - Implemented the random effects design matrix, contributed to 
#                  matrix manipulations, and ensured the accuracy of model computations.
#
# - Sagar Udasi: S2606876
#           15% - Contributed to the overall testing of the code.

# -------------------------- PROBLEM STATEMENT AND APPROACH --------------------------- #
# This code implements a Linear Mixed Model (LMM) for estimating fixed and random 
# effect coefficients. It constructs design matrices, applies optimization 
# to minimize the negative log-likelihood, and estimates random effect variances 
# and fixed effect coefficients. Key functions like LMMsetup, generate_Z_matrix,
# and compute_beta handle model setup, random effect matrix generation, 
# and fixed effect estimation, supporting multiple random effect groups.
# The output includes estimates for both fixed effects and random effect variances.

# ------------------------------------------ CODE ------------------------------------- #
## 'lmm' function is the main function for fitting a Linear Mixed Model (LMM)
lmm <- function(form, dat, ref = list()) {
  
  # Check if 'form' is a valid formula (model specification)
  if (!inherits(form, "formula")) {
    stop("'form' must be a valid formula.")
  }
  
  # Check if 'dat' is a data frame (input data)
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data frame.")
  }
  
  # Check if 'ref' is a list of valid variable names (for random effects)
  if (!is.list(ref)) {
    stop("'ref' must be a list.")
  }
  
  # Check if all variables in each group of 'ref' exist in the data columns
  for (group in ref) {
    if (!all(group %in% colnames(dat))) {
      stop(paste("All variables in a group of 'ref' must exist in the data columns. Invalid group: ", paste(group, collapse = ":")))
    }
  }
  
  # If no random effects specified (ref is empty), fit a simple linear model
  if (length(ref) == 0) {
    lm_model <- lm(form, data = dat) # Fit linear model
    beta_out <- coef(lm_model) # Extract coefficients from the linear model
    names(beta_out) <- names(coef(lm_model)) # Assign names to the coefficients
    return(list(beta = beta_out, theta = NULL)) # Return beta estimates, no random effects
  }

  # Initialize model matrices and parameters using function 'LMMsetup'
  setup <- LMMsetup(form, dat, ref)

  # Define the initial values of the theta parameter vector, assuming initial values are all 0s
  # The first element is log(sigma) and the remaining elements are the log standard deviations of the random effect variances
  initial_theta <- rep(0, length(ref) + 1)
  

  # Use 'optim' to minimize 'LMMprof', obtaining optimal theta values
  # The 'par' argument sets initial values of theta, and 'LMMprof' is the objective function
  # 'setup' and 'group_sizes' provide data structure details required by 'LMMprof'
  theta_opt <- optim(par = initial_theta, fn = LMMprof, setup = setup)
  
  # Calculate beta estimates using the optimized theta values from 'theta_opt'
  beta <- compute_beta(theta_opt$par, setup)
  
  # Extract the estimated beta coefficients as a vector for clarity
  beta_out <- as.vector(beta$beta_hat)
  # Assign meaningful names to beta coefficients: "Residual" for the residual term,
  # and concatenated random effect names (e.g., "Worker:Machine") for each random effect
  names(beta_out) <- c("Residual", sapply(ref, function(x) paste(x, collapse = ":")))

  # This transformation converts log-variance terms back to variance scale for interpretation
  theta_out <- exp(theta_opt$par)
  # Assign names to theta parameters in the same way as beta coefficients
  names(theta_out) <- c("Intercept", sapply(ref, function(x) paste(x, collapse = ":")))
  # Return a list containing both 'beta' and 'theta' estimates
  
  return(list(beta = beta_out, theta = theta_out))
}

## 'LMMsetup' function prepares the necessary matrices and decompositions 
## to facilitate the estimation of a Linear Mixed Model.
LMMsetup <- function(form, dat, ref) {
  
  # Check if 'form' is a valid formula
  if (!inherits(form, "formula")) {
    stop("'form' must be a valid formula.")
  }
  
  # Check if 'dat' is a data frame
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data frame.")
  }
  # Check if 'ref' is a list and each item is a valid vector of variable names
  if (!is.list(ref)) {
    stop("'ref' must be a list.")
  }
  # Check if all variables in each group of 'ref' exist in the data columns
  for (group in ref) {
    if (!all(group %in% colnames(dat))) {
      stop(paste("All variables in a group of 'ref' must exist in the data columns. Invalid group: ", paste(group, collapse = ":")))
    }
  }
  
  # Creates the design matrix X for fixed effects from 'form' and 'dat'
  X <- model.matrix(form, dat)
  # Extracts the response vector y from the model formula and data
  y <- model.response(model.frame(form, dat))
  # Calls function 'generate_Z_matrix' to create the random effects design matrix Z 
  # based on the list of grouping variables in 'ref'
  Z <- generate_Z_matrix(ref, dat)
  # QR decomposition of Z
  z_decom <- qr(Z)
  # Extract R from the QR decomposition
  R <- qr.R(z_decom)
  # Number of observations
  n <- nrow(Z)
  # Number of random effect parameters
  p <- ncol(Z)
  # Multiplies Q^T (from the QR decomposition of Z) with the matrix X.
  qtx <- qr.qty(z_decom,X)
  # Multiplies Q^T (from the QR decomposition of Z) with the response vector y.
  qty <- qr.qty(z_decom,y)
  
  # 'sapply' applies a function to each element in 'ref' and returns a vector.
  # This function calculates the number of unique combinations 
  # of the variables in 'vars' within the dataset 'dat' for each item in 'ref'.
  group_sizes <- sapply(ref, function(vars) {
    nrow(unique(dat[vars]))
  })
  
  # Return structured list with matrices and parameters for model fitting
  return(list(X = X, y = y, Z = Z, z_decom = z_decom, R = R, n = n, p = p, qtx = qtx, qty =  qty, group_sizes = group_sizes))
}

## 'LMMprof' function is to compute the negative log-likelihood for a Linear Mixed Model.
LMMprof <- function(theta, setup) {
  
  # Check if 'theta' is a numeric vector
  if (!is.numeric(theta) || length(theta) != length(setup$group_sizes) + 1) {
    stop("'theta' must be a numeric vector of length equal to the number of random effects plus 1.")
  }
 
  # Extract matrices and dimensions from setup
  X <- setup$X # Fixed effects design matrix (X)
  y <- setup$y # Response vector (y)
  Z <- setup$Z # Random effects design matrix (Z)
  n <- setup$n # Number of observations (n)
  p <- setup$p # Number of random effect parameters (p)
  z_decom <- setup$z_decom # QR decomposition of (Z)
  R <- setup$R # Upper triangular matrix from QR decomposition
  group_sizes <- setup$group_sizes # Sizes of the random effect groups
  
  # Generate the psi matrix (diagonal matrix of random effect variances)
  psi <- generate_psi(theta,group_sizes)

  # Calculate log determinant term |V|(determinant of covariance matrix) using Cholesky factorization
  # 'theta_sigma' is the variance-covariance matrix of random effects and residuals
  # Residual variance + random effect variances
  theta_sigma <- R %*% psi %*% t(R) +  diag(rep(exp(theta[1])^2,p)) 
  # Cholesky decomposition of theta_sigma
  L <- chol(theta_sigma) 
  # Log determinant |V| 
  log_def_V <- 2 * sum(log(diag(L))) + 2*(n-p)*theta[1] 
  
  # Estimate beta and other values required for likelihood
  res <- compute_beta(theta,setup,L)
  # Fixed effect estimates (beta)
  beta_hat <- res$beta_hat
  # Matrix for mixing the random effects and fixed effects
  mix_L <- res$mix_L
    
  # Calculate the weighted residual sum
  # Residuals after fixed effect adjustment
  residual <- y - X %*% beta_hat
  # QR decomposition of residuals
  qtresidual <- qr.qty(z_decom,residual)
  # Weighted residuals using mix_L
  mix_L_qtresidual <- mix_L %*% qtresidual
  # Sum of squared weighted residuals
  weighted_residual_sum <- t(mix_L_qtresidual) %*% mix_L_qtresidual
  
  # Compute the negative log-likelihood (LMM objective function)
  neg_log_likelihood <- 0.5 * (log_def_V  + weighted_residual_sum)
  
  # Return the negative log-likelihood value
  neg_log_likelihood
}

## 'generate_Z_matrix' function creates a design matrix (Z) for random effects based on variable groups in 'ref'
generate_Z_matrix <- function(ref, dat) {
  
  # Check if 'ref' is a non-empty list
  if (length(ref) == 0) {
    stop("'ref' must contain at least one variable group.")
  }
  
  # Ensure all variables in 'ref' exist in the data columns
  for (group in ref) {
    if (!all(group %in% colnames(dat))) {
      stop(paste("All variables in a group of 'ref' must exist in the data columns. Invalid group: ", paste(group, collapse = ":")))
    }
  }
  
  # Generate model matrices for each variable group in 'ref'
  Z_blocks <- lapply(ref, function(vars) {
    # Create formula string for the interaction of variables in 'vars' without an intercept term ('- 1')
    # 'paste(vars, collapse = ":")' combines the variable names into a single string, separated by ':' for interaction terms.
    # The '~' symbol indicates the start of a model formula in R. The '- 1' part removes the intercept term.
    formula_str <- paste("~", paste(vars, collapse = ":"), "- 1")
    # This generates a matrix where each column corresponds to the levels of the factors specified in 'vars' (with no intercept).
    # 'contrasts.arg = lapply(dat[vars], contrasts, contrasts = FALSE)' ensures that no extra contrasts are applied.
    # It applies default factor contrasts to categorical variables in the 'vars' group (if any).
    model.matrix(as.formula(formula_str), dat, contrasts.arg = lapply(dat[vars], contrasts, contrasts = FALSE))
  })
  
  # Combine all model matrices (blocks) into one complete Z matrix
  Z <- do.call(cbind, Z_blocks)
  
  return(Z)
}

## 'compute_beta' function is designed to compute the estimated fixed effect coefficients (beta_hat) for a Linear Mixed Model.
compute_beta<- function(theta,setup,L = NULL){
  
  # Check if 'theta' is a numeric vector
  if (!is.numeric(theta)) {
    stop("'theta' must be a numeric vector.")
  }
  
  # Extract relevant matrices and parameters from 'setup'
  X <- setup$X # Fixed effects design matrix (X)
  y <- setup$y # Response vector (y)
  n <- setup$n # Number of observations (n)
  p <- setup$p # Number of random effect parameters (p)
  z_decom <- setup$z_decom # QR decomposition of (Z)
  group_sizes <- setup$group_sizes # Sizes of the random effect groups
  # Multiplies Q^T (from the QR decomposition of Z) with the matrix X.
  qtx <- qr.qty(z_decom,X)
  # Multiplies Q^T (from the QR decomposition of Z) with the response vector y.
  qty <- qr.qty(z_decom,y)
  
  # If 'L' is not provided, compute it
  if(is.null(L)){ 
  # Upper triangular matrix from QR decomposition
  R <- setup$R 
  # Compute psi (random effect variances)
  psi <-  generate_psi(theta,group_sizes)
  # Construct theta's covariance matrix
  theta_sigma <- R %*% psi %*% t(R) +  diag(rep(exp(theta[1])^2,p))
  # Cholesky decomposition of the covariance matrix
  L <- chol(theta_sigma)
  }
  
  ## Create complete mix matrix
  # Solve L * X = I for X, where I is the identity matrix
  L_inv<-backsolve(L,diag(ncol(L)))
  # Compute the inverse of theta's covariance matrix
  theta_sigma_i <- (L_inv) %*% t(L_inv)
  # Create a diagonal matrix for sigma (residual variance)
  # Residual variance (inverse of standard deviation squared)
  sigma <- diag(rep(exp(theta[1])^(-2),n-p)) 
  # Initialize a matrix for the mixed effects model
  mix_matrix <- matrix(0, nrow = n, ncol = n)
  # Fill the top-left block of the matrix with theta's covariance inverse
  mix_matrix[1:p, 1:p] <-  theta_sigma_i
  # Fill the bottom-right block with sigma (residual variance)
  mix_matrix[(p + 1):n, (p + 1):n] <- sigma
  
  
  ## Compute beta
  # Compute the Cholesky decomposition of the mixed effects matrix
  mix_L <- chol(mix_matrix)
  # Compute transformed versions of qtx and qty using the Cholesky factor
  mix_L_qtx <- mix_L %*% qtx # Multiply by qtx (fixed effect components)
  mix_L_qty <- mix_L %*% qty # Multiply by qty (response components)
  # Perform Cholesky factorization on the matrix of X'X
  L_XTWX <- chol(t(mix_L_qtx) %*% mix_L_qtx)
  # Compute the inverse of X'X using backsubstitution
  # Solve L_XTWX * X = I for X
  L_XTWX_inv <- backsolve(L_XTWX,diag(ncol(L_XTWX))) 
  # Estimate the beta coefficients
  beta_hat<-  L_XTWX_inv %*% t(L_XTWX_inv) %*% t(mix_L_qtx) %*% mix_L_qty
  # Return the estimated beta and Cholesky factor
  # The purpose of return 'mix_L' is to avoid computing repeatedly
  return(list(beta_hat = beta_hat, mix_L = mix_L)) # 
}

## ‘generate_psi’ function is to generate the covariance matrix (psi) 
## for the random effects in a Linear Mixed Model. 
generate_psi <- function(theta, group_sizes) {
  
  # Check if 'theta' is a numeric vector
  if (!is.numeric(theta)) {
    stop("'theta' must be a numeric vector.")
  }
  
  # Check if 'group_sizes' is a vector
  if (!is.vector(group_sizes)) {
    stop("'group_sizes' must be a vector.")
  }
  
  # Extract random effect variances (excluding the first element)
  random_effect_vars <- exp(theta[-1])^2
  
  # Check if the length of 'random_effect_vars' matches 'group_sizes'
  if (length(random_effect_vars) != length(group_sizes)) {
    stop("Length of 'group_sizes' must match the length of random effect variances in 'theta'.")
  }
  
  # Extract random effect variances (excluding the first element of theta)
  random_effect_vars <- exp(theta[-1])^2
  
  
  # Create a diagonal matrix where each random effect variance is repeated according to group sizes
  # mapply() applies a function across pairs of elements from random_effect_vars and group_sizes 
  # and repeats each variance according to the group size.
  # unlist() flattens the resulting list of repeated variances into a single vector.
  diag_values <- unlist(mapply(function(var, size) rep(var, size), random_effect_vars, group_sizes))
  
  # Construct a diagonal matrix with the repeated variances
  psi <- diag(diag_values)
  
  # Return the generated psi matrix
  return(psi)
}

# Test
library(nlme)
library(lme4)
lmm(score ~ Machine,Machines,list(c("Worker"),c("Worker","Machine")))
lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),data=Machines,
     REML=FALSE)
