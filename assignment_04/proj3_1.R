
# Main function for fitting a Linear Mixed Model (LMM)
lmm <- function(form, dat, ref = list()) {
  # 使用 LMMsetup 初始化矩阵和其他参数
  # Calculate the number of unique combinations for each random effect group based on variables

  # Initialize model matrices and parameters using LMMsetup
  setup <- LMMsetup(form, dat, ref)
  # 定义初始 theta 参数，其中第一个元素为 log(sigma)，其余为随机效应的 log 标准差
  # Initialize theta parameters: log(sigma) and log SDs for random effects
  initial_theta <- rep(0, length(ref) + 1)  # 假设初始为 0
  
  # 使用 optim 优化 LMMprof 函数，以获得 theta 的最优值
  # Use `optim` to minimize `LMMprof`, obtaining optimal theta values
  theta_opt <- optim(par = initial_theta, fn = LMMprof, setup = setup)
  # Calculate beta estimates using optimized theta values
  beta <- compute_beta(theta_opt$par, setup)
  # Return estimated beta coefficients and transformed theta parameters
  beta_out <- as.vector(beta$beta_hat)
  names(beta_out) <- c("Intercept",sapply(ref, function(x) paste(x, collapse = ":")))

  theta_out <- exp(theta_opt$par)
  names(theta_out) <- c("Intercept",sapply(ref, function(x) paste(x, collapse = ":")))
  
  return(list(beta = beta_out, theta = theta_out))
}

# 辅助函数：LMMsetup
LMMsetup <- function(form, dat, ref) {
  # 解析模型公式，生成固定效应设计矩阵 X 和响应变量 y
  # From model formula, creating fixed effects matrix X and response variable y
  X <- model.matrix(form, dat)
  y <- model.response(model.frame(form, dat))
  # Create random effects matrix Z based on the ref list.
  Z <- generate_Z_matrix(ref, dat)
  # QR decomposition of Z
  z_decom <- qr(Z)
  # Extract R from the QR decomposition
  R <- qr.R(z_decom)
  # Number of observations
  n <- nrow(Z)
  # Number of random effect parameters
  p <- ncol(Z)
  qtx <- qr.qty(z_decom,X)
  qty <- qr.qty(z_decom,y)
  group_sizes <- sapply(ref, function(vars) {
    # 按照随机效应组的变量进行分组，统计唯一组合的数量
    nrow(unique(dat[vars]))
  })
  # 返回包含 X, y, Z 及其他需要的结构的列表
  # Return structured list with matrices and parameters for model fitting
  return(list(X = X, y = y, Z = Z, z_decom = z_decom, R = R, n = n, p = p,qtx = qtx,qty =  qty,group_sizes = group_sizes))
}

# Object Function：LMMprof
LMMprof <- function(theta, setup) {
  # 从 setup 中提取矩阵信息
  # Extract matrices and dimensions from setup
  X <- setup$X
  y <- setup$y
  Z <- setup$Z
  n <- setup$n
  p <- setup$p
  z_decom <- setup$z_decom
  R <- setup$R
  group_sizes <- setup$group_sizes
  ####psi
  psi <-  generate_psi(theta,group_sizes)

  
  # Calculate log determinant |V| term using Cholesky factorization
  theta_sigma <- R %*% psi %*% t(R) +  diag(rep(exp(theta[1])^2,p))
  L <- chol(theta_sigma)
  log_def_V <- 2 * sum(log(diag(L))) + 2*(n-p)*theta[1]
  
  #计算beta估计
  # Estimate beta and other values required for likelihood
  res <- compute_beta(theta,setup,L)
  beta_hat <- res$beta_hat
  mix_L <- res$mix_L
    
  #计算weighted residual sum
  residual <- y - X %*% beta_hat
  qtresidual <- qr.qty(z_decom,residual)
  mix_L_qtresidual <- mix_L %*% qtresidual
  weighted_residual_sum <- t(mix_L_qtresidual) %*% mix_L_qtresidual
  
  ##计算负log
  neg_log_likelihood <- 0.5 * (log_def_V  + weighted_residual_sum)
  neg_log_likelihood
}

# 生成 Z 矩阵的辅助函数
generate_Z_matrix <- function(ref, dat) {
  # 为 ref 中的每个变量组合生成模型矩阵块
  Z_blocks <- lapply(ref, function(vars) {
    # 创建公式字符串，确保没有截距项
    formula_str <- paste("~", paste(vars, collapse = ":"), "- 1")
    # 生成模型矩阵，确保不会添加多余的列
    model.matrix(as.formula(formula_str), dat, contrasts.arg = lapply(dat[vars], contrasts, contrasts = FALSE))
  })
  
  # 合并各模型矩阵块，形成完整的 Z 矩阵
  Z <- do.call(cbind, Z_blocks)
  return(Z)
}

compute_beta<- function(theta,setup,L = NULL){
  
  X <- setup$X
  y <- setup$y
  n <- setup$n
  p <- setup$p
  z_decom <- setup$z_decom
  qtx <- setup$qtx
  qty <- setup$qty
  group_sizes <- setup$group_sizes
  
  if(is.null(L)){ 
  R <- setup$R
  psi <-  generate_psi(theta,group_sizes)
  theta_sigma <- R %*% psi %*% t(R) +  diag(rep(exp(theta[1])^2,p))
  L <- chol(theta_sigma)
  }
  
  #创建完整的混合矩阵
  theta_sigma_i <- chol2inv(L)
  sigma <- diag(rep(exp(theta[1])^(-2),n-p)) 
  #创建空矩阵
  mix_matrix <- matrix(0, nrow = n, ncol = n)
  # 将 theta_sigma 放到左上角
  mix_matrix[1:p, 1:p] <-  theta_sigma_i
  # 将 sigma 放到右下角
  mix_matrix[(p + 1):n, (p + 1):n] <- sigma
  
  
  #计算beta
  
  mix_L <- chol(mix_matrix)
  mix_L_qtx <- mix_L %*% qtx
  mix_L_qty <- mix_L %*% qty
  beta_hat<-  chol2inv(chol(t(mix_L_qtx) %*% mix_L_qtx)) %*% t(mix_L_qtx) %*% mix_L_qty
  
  return(list(beta_hat = beta_hat, mix_L = mix_L))
}

generate_psi <- function(theta, group_sizes) {
  # 提取每个随机效应组的方差，排除第一个元素
  random_effect_vars <- exp(theta[-1])^2
  
  # 检查 group_sizes 和 random_effect_vars 是否匹配
  if (length(random_effect_vars) != length(group_sizes)) {
    stop("Length of 'group_sizes' must match the length of random effect variances in 'theta'.")
  }
  
  # 创建对角矩阵，每个方差重复指定次数
  diag_values <- unlist(mapply(function(var, size) rep(var, size), random_effect_vars, group_sizes))
  psi <- diag(diag_values)
  return(psi)
}

library(nlme)
library(lme4)
lmm(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),data=Machines,
     REML=FALSE)
