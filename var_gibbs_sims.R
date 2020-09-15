########################################################################
# This example assumes q = 1, C = Identity, D = 0, m = 0, a = 0, X = 1_n
########################################################################


# FUNCTION DEFS
get_L_fix <- function(A)
{
  r <- ncol(A)
  L <- sum(A^2) + r * r
  return(L)
}

get_eps_fix <- function(Y, Z, S, Tconst, logarithm = TRUE){
  r <- ncol(Y)
  n <- nrow(Y)
  e_M <- eigen(n * S)
  out <- sum(log(e_M$values))
  h <- norm(scale(Y, scale = F), "2") + 
    sqrt(Tconst) * norm(scale(Z, scale = F), "2")
  out <- out - 2 * r * log(h)
  out <- out * (n - r - 2) / 2
  if(!logarithm) out <- exp(out)
  return(out)
}

get_lambda <- function(A, n){
  r <- ncol(A)
  out <- (r + sum(A^2)) / (n - 2 * r - 3)
  return(out)
}

get_L <- function(lambda, S, n){
  return(lambda * n * sum(diag(S)))
}

get_epsilon <- function(S, Tconst, n, logarithm = T){
  r <- ncol(S)
  e_M <- eigen(n * S)
  out <- sum(log(e_M$values)) - sum(log(e_M$values + Tconst))
  out <- out * (n - 2 - r) / 2
  if(!logarithm) out <- exp(out)
  return(out)
}

# SIMULATION
FIG1 <- FALSE # Produces Fig 2 if FALSE
PDF <- FALSE  # Plot in R instead of producing pdf if FALSE
out_dir <- "~/GDrive/Research/var_gibbs/sims/" # end in "/", used if PDF = TRUE

set.seed(3)
r <- 10
n_max <- ifelse(FIG1, 200, 4000)
n_seq <- floor(seq(2 * r + 2 + 10, n_max, length.out = 100))


# Data generating parameters
A <- matrix(runif(r^2, -.5, .5), r, r)
A <- A + t(A) + diag(1, r)
A <- A / (norm(A, "2") + 0.1)
sigma <- 1
Sigma <- diag(sigma^2, r)
R <- chol(Sigma)


# Generate VAR(1)
W <- matrix(rnorm(n_max * r), n_max, r)
for(ii in 2:n_max){
  W[ii, ] <- W[ii - 1, ] %*% A + t(crossprod(R, W[ii, ]))
}

# Allocate
lambda <- rep(0, length(n_seq))
L <- rep(0, length(n_seq))
L_fix <- L
log_eps <- rep(0, length(n_seq))
log_eps_fix <- log_eps
lambda_pop <- rep(0, length(n_seq))
L_pop <- rep(0, length(n_seq))
log_eps_pop <- rep(0, length(n_seq))

# Compute quantities to plot for range of n
for(ii in 1:length(n_seq)){
  Y <- W[2:n_seq[ii], ]
  Z <- W[1:(n_seq[ii] - 1), ]
  fit <- lm(Y ~ Z)
  A_hat <- coef(fit)[-1, , drop = F]
  n <- nobs(fit)
  S <- crossprod(residuals(fit)) / n

  # lambda
  lambda[ii] <- get_lambda(A = A_hat, n = n)
  lambda_pop[ii] <- get_lambda(A = A, n = n)
  
  # L
  L[ii] <- get_L(lambda = lambda[ii], S = S, n = n)
  L_pop[ii] <-  get_L(lambda = lambda_pop[ii], S = Sigma, n = n)
  
  # epsilon
  Tval <- 2 * L[ii] / (1 - lambda[ii]) + 1e-6
  if(Tval <= 0){
    log_eps[ii] <- NA
  } else{
    log_eps[ii] <- get_epsilon(S = S, Tconst = Tval, n = n, logarithm = T)
  }
  
  Tval_pop <- 2 * L_pop[ii] / (1 - lambda_pop[ii]) + 1e-6
  if(Tval_pop <= 0){
    log_eps_pop[ii] <- NA
  } else{
    log_eps_pop[ii] <- get_epsilon(S = Sigma, Tconst = Tval_pop, n = n, logarithm = T)
  }
  
  # Fixed data quantities
  L_fix[ii] <- get_L_fix(A_hat)
  Tval_fix <- 2 * L_fix[ii] + 1e-6
  log_eps_fix[ii] <- get_eps_fix(Y = Y, Z = Z, S = S, Tconst = Tval_fix)
}

# PRODUCE PLOTS
par(cex.axis = 1.3, cex.lab = 1.3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))

if(FIG1){
  if(PDF) pdf(paste(out_dir, "fig_1.pdf", sep = ""), width = 12.5, height = 5)
  par(mfrow = c(1, 2))
  plot(n_seq, lambda, xlab = "n", ylab = expression(lambda[n]))
  abline(v = n_seq[min(which(lambda < 1))], lwd = 2)
  lines(n_seq, lambda_pop, col = "red", lwd = 2, lty = 2)
  plot(n_seq, L,
       xlab = "n",
       ylab = expression(L[n]),
       ylim = c(min(c(L_pop, L, L_fix)), max(c(L, L_fix))))
  abline(a = sigma, b = 0, lwd = 2)
  lines(n_seq, L_pop, col = "red", lwd = 2, lty = 2)
  lines(n_seq, L_fix, col = "darkgreen", lwd = 2, lty = 3)
  if(PDF) dev.off()
  par(mfrow = c(1, 1))
}

if(!FIG1){
  if(PDF) pdf(paste(out_dir, "fig_2.pdf", sep = ""), width = 12.5, height = 5)
  true_log_eps <- -r^2 * (r + sum(A^2))
  plot(n_seq, log_eps,
       xlab = "n",
       ylab = expression(paste("log ", epsilon[n])),
       ylim = c(true_log_eps, true_log_eps / 1.5))
  lines(n_seq,log_eps_pop, col = "red", lwd = 2, lty = 2)
  lines(n_seq,log_eps_fix, col = "darkgreen", lwd = 2, lty = 3)
  if(PDF) dev.off()
}
