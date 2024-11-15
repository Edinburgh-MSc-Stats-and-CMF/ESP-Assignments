# ######################### INDEPENDENT SUBMISSION (V3) ########################
# - Xu Guo: S2743905
#           45% - Led the development of the linear mixed model structure 
#                  and optimized the fitting process using optimization techniques.
#                        
# - Xie Ying: S2510813
#           35% - Implemented the random effects design matrix, contributed to 
#                  matrix manipulations, and ensured the accuracy of model computations.
#
# - Sagar Udasi: S2606876
#           20% - Contributed to the overall testing of the code and reconstructed the code.
# ##############################################################################


# ==============================================================================
# PROBLEM STATEMENT:
# The goal is to implement a Linear Mixed Model (LMM) for regression with both  
# fixed and random effects. The model should validate inputs, generate necessary 
# design matrices, and optimize variance parameters using a custom approach, with 
# an option for comparison to standard LMM implementations (e.g., from 'lme4').


# SOLUTION APPROACH:
# The random effects design matrix (Z) is generated based on specified reference 
# groups, which indicate the hierarchical structure of the data. To estimate the 
# variance components, a log-likelihood function is maximized using an 
# optimization procedure. The optimization estimates both the fixed effects 
# (beta) and the variance components (theta), which represent the residual 
# variance and the variance due to random effects. The model is fitted using QR 
# decomposition for efficient computation, and the results are interpreted in 
# terms of the estimated coefficients for fixed effects and the variance 
# components for random effects.
# ==============================================================================


# ==============================================================================
# ============================== HELPER FUNCTIONS ==============================
# ==============================================================================


# Function to validate user inputs
validate_inputs <- function(form, dat, ref) {
  "
  This function checks if the provided arguments are valid for the LMM implementation.
 
  Args:
      form (formula): The formula for the linear mixed model.
      dat (data.frame): The data frame containing the model variables.
      ref (list): A list of grouping variables for random effects.
 
  Returns:
      NULL: If all inputs are valid.
      stop: Raises an error message if any input is invalid.
  "
  # Check if 'form' is a valid formula (model specification)
  if (!inherits(form, "formula")) stop("'form' must be a valid formula.")
  # Check if 'dat' is a data frame (input data)
  if (!is.data.frame(dat)) stop("'dat' must be a data frame.")
  # Check if 'ref' is a list of valid variable names (for random effects)
  if (!is.list(ref)) stop("'ref' must be a list.")
  
  # Check if all grouping variables exist in the data frame
  for (group in ref) {
    if (!all(group %in% colnames(dat))) {
      stop(paste("Invalid group variables:", paste(group, collapse = ":")))
    }
  }
}


# Function to generate design matrix for random effects
generate_Z_matrix <- function(ref, dat) {
  "
  This function creates a design matrix (Z) for the random effects based on the 
  grouping variables.
 
  Args:
      ref (list): A list of grouping variables for random effects.
      dat (data.frame): The data frame containing the model variables.
 
  Returns:
      matrix: The design matrix (Z) for random effects.
 
  Raises:
      stop: If 'ref' is empty (no random effects specified).
  "
  # Create model matrices for each group of random effects
  Z_blocks <- lapply(ref, function(vars) {
    # Create a formula for the current group, joining the variable names with ':' 
    # to represent interactions between the variables in that group 
    # (e.g., "variable1:variable2"). We use '- 1' to exclude the intercept term, 
    # as we want to capture the specific effects of these variables.
    formula_str <- paste("~", paste(vars, collapse = ":"), "- 1")
    
    # Generate the model matrix for the group using the formula created above.
    # The 'model.matrix' function converts the formula into a design matrix, 
    # where the columns correspond to the interactions between the group variables.
    # The 'contrasts.arg' argument ensures that factors are treated as categorical, 
    # without applying default contrasts.
    model.matrix(as.formula(formula_str), dat, 
                 contrasts.arg = lapply(dat[vars], contrasts, contrasts = FALSE))
  })
  
  # Combine the model matrices for all random effects groups
  do.call(cbind, Z_blocks)
}


# ==============================================================================
# =========================== LMM IMPLEMENTATION ===============================
# ==============================================================================


# Function to set up the components of the Linear Mixed Model (LMM)
LMMsetup <- function(form, dat, ref) {
  "
  This function prepares the necessary components for the LMM model fitting process.
 
  Args:
      form (formula): The formula for the linear mixed model.
      dat (data.frame): The data frame containing the model variables.
      ref (list): A list of grouping variables for random effects.
 
  Returns:
      list: A list containing various model components (X, y, Z, etc.).
  "
  
  
  validate_inputs(form, dat, ref)
  
  # Creates the design matrix X for fixed effects from 'form' and 'dat'
  X <- model.matrix(form, dat)
  # Extracts the response vector y
  y <- model.response(model.frame(form, dat))
  
  # Calls function 'generate_Z_matrix' to create the random effects design matrix Z 
  # based on the list of grouping variables in 'ref'
  Z <- generate_Z_matrix(ref, dat)
  
  # Perform QR decomposition of Z for efficiency
  z_decom <- qr(Z)
  R <- qr.R(z_decom)
  
  # Pre-compute QR transformations for X and y
  qtx <- qr.qty(z_decom, X)
  qty <- qr.qty(z_decom, y)
  
  # Calculate model dimensions
  n <- nrow(Z)
  p <- ncol(Z)
  # Calculate group size
  # 'sapply' applies a function to each element in 'ref' and returns a vector.
  # This function calculates the number of unique combinations 
  # of the variables in 'vars' within the dataset 'dat' for each item in 'ref'.
  group_sizes <- sapply(ref, function(vars) nrow(unique(dat[vars])))
  
  # Return a list containing all the model components
  list(X = X, y = y, Z = Z, z_decom = z_decom, R = R,
       n = n, p = p, qtx = qtx, qty = qty, group_sizes = group_sizes)
}


# Function to calculate the negative log-likelihood of the LMM model
LMMprof <- function(theta, setup) {
  "
  This function calculates the negative log-likelihood of the linear mixed model 
  given the model parameters.
 
  Args:
      theta (numeric vector): A vector of model parameters (log-variance components).
      setup (list): A list containing the model components (X, y, Z, etc.).
 
  Returns:
      numeric: The negative log-likelihood of the model.
  "
  # Generate variance-covariance matrix for random effects
  psi <- diag(unlist(mapply(function(var, size) rep(exp(var)^2, size),
                            theta[-1], setup$group_sizes)))
  
  # Calculate log determinant term |V|(determinant of covariance matrix) using Cholesky factorization
  # 'theta_sigma' is the variance-covariance matrix of random effects and residuals
  # Residual variance + random effect variances
  theta_sigma <- setup$R %*% psi %*% t(setup$R) + 
    diag(rep(exp(theta[1])^2, setup$p))
  # Cholesky decomposition of theta_sigma
  L <- chol(theta_sigma)
  
  # Calculate the log-determinant of the variance-covariance matrix
  log_det_V <- 2 * sum(log(diag(L))) + 2 * (setup$n - setup$p) * theta[1]
  
  # Compute the inverse of theta_sigma and set up the mixed effects matrix
  # Solve L * X = I for X, where I is the identity matrix
  L_inv <- backsolve(L, diag(ncol(L))) 
  # Compute the inverse of theta's covariance matrix
  theta_sigma_i <- L_inv %*% t(L_inv) 
  # Initialize a matrix for the mixed effects model
  mix_matrix <- matrix(0, nrow = setup$n, ncol = setup$n)
  # Fill the top-left block of the matrix with theta's covariance inverse
  mix_matrix[1:setup$p, 1:setup$p] <- theta_sigma_i
  # Fill the bottom-right block with sigma (residual variance)
  mix_matrix[(setup$p + 1):setup$n, (setup$p + 1):setup$n] <- 
    diag(rep(exp(theta[1])^(-2), setup$n - setup$p))
  
  # Perform Cholesky decomposition on the mixed effects matrix
  mix_L <- chol(mix_matrix)
  mix_L_qtx <- mix_L %*% setup$qtx # Multiply by qtx (fixed effect components)
  mix_L_qty <- mix_L %*% setup$qty # Multiply by qty (response components)
  
  # Perform Cholesky factorization on the matrix of X'X
  L_XTWX <- chol(t(mix_L_qtx) %*% mix_L_qtx)
  # Compute the inverse of X'X using backsubstitution
  # Solve L_XTWX * X = I for X
  L_XTWX_inv <- backsolve(L_XTWX, diag(ncol(L_XTWX)))
  # Estimate the beta coefficients
  beta_hat <- L_XTWX_inv %*% t(L_XTWX_inv) %*% t(mix_L_qtx) %*% mix_L_qty
  
  # Calculate residuals and their weighted sum for the log-likelihood
  # Residuals after fixed effect adjustment
  residual <- setup$y - setup$X %*% beta_hat
  # QR decomposition of residuals
  qtresidual <- qr.qty(setup$z_decom, residual)
  # Weighted residuals using mix_L
  mix_L_qtresidual <- mix_L %*% qtresidual
  # Sum of squared weighted residuals
  weighted_residual_sum <- sum(mix_L_qtresidual^2)
  
  # Store beta_hat as an attribute of the residual sum for later retrieval
  attr(weighted_residual_sum, "beta_hat") <- beta_hat
  
  # Return the negative log-likelihood
  0.5 * (log_det_V + weighted_residual_sum)
}


# Function to fit the Linear Mixed Model (LMM) using optimized matrix operations
lmm <- function(form, dat, ref = list()) {
  "
  Fits a linear mixed-effects model (LMM) to the specified data.
 
  Args:
    form: A formula object specifying the fixed effects model.
    dat: A data frame containing the variables in the model.
    ref: A list of grouping variables for random effects.
 
  Returns:
    A list containing:
      - beta: A vector of estimated fixed effects coefficients.
      - theta: A vector of estimated variance components for the random effects 
               and residual error.
 
  The function first validates the input arguments. If no random effects are 
  specified (i.e., `ref` is empty), it falls back to a simple linear regression 
  using `lm`. For LMMs, it sets up the necessary model components, including 
  fixed effects, random effects, and variance-covariance matrices. The model is 
  fitted by optimizing the log-likelihood function with respect to the variance 
  components (theta). The optimization is performed using the `optim` function.
  
  Once the optimal variance components are found, the fixed effects coefficients 
  (beta) are estimated using the generalized least squares approach. The function 
  returns a list containing the estimated fixed effects and variance components.
  "
  validate_inputs(form, dat, ref)
  
  # If no random effects specified (ref is empty), fit a simple linear model
  if (length(ref) == 0) {
    lm_model <- lm(form, data = dat) 
    beta_out <- coef(lm_model) 
    names(beta_out) <- names(coef(lm_model)) 
    return(list(beta = beta_out, theta = NULL)) 
  }
  
  # Set up the components for the LMM (fixed effects, random effects, etc.)
  setup <- LMMsetup(form, dat, ref)
  
  # Define the initial values of the theta parameter vector, assuming initial values are all 0s
  # The first element is log(sigma) and the remaining elements are the log standard deviations of the random effect variances
  initial_theta <- rep(0, length(ref) + 1)
  opt_result <- optim(par = initial_theta, fn = LMMprof, setup = setup)
  
  # Extract the optimized beta estimates and variance parameters
  beta_hat <- attr(LMMprof(opt_result$par, setup), "beta_hat")
  theta_out <- exp(opt_result$par)
  
  # Name the beta estimates and variance components for clarity
  names(beta_hat) <- colnames(setup$X)
  names(theta_out) <- c("Residual", 
                        sapply(ref, function(x) paste(x, collapse = ":")))
  
  # Return the results as a list of beta estimates and variance components
  list(beta = as.vector(beta_hat), theta = theta_out)
}
