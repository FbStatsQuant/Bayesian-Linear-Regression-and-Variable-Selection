# Horseshoe linear regression on synthetic data (Gibbs sampler, base R)
# Model and updates follow ../Notes/horseshoe_linear_regression.tex
# Run from the Scripts folder; metrics and plots are written to results/.

set.seed(1)
dir.create("results", showWarnings = FALSE)

# ---- Simulate data ----------------------------------------------------------
n <- 100
p <- 50
beta_true <- rep(0, p)
beta_true[c(1, 5, 10, 20, 35)] <- c(3, -2, 1.5, 1, -0.5)

Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
X <- scale(X) * sqrt(n / (n - 1))       # centered, X_j'X_j = n
y <- 2 + X %*% beta_true + rnorm(n)
y <- drop(y)

# ---- Gibbs sampler ----------------------------------------------------------
rinvgamma <- function(k, a, b) b / rgamma(k, shape = a, rate = 1)

horseshoe <- function(y, X, n_iter = 20000, n_burn = 2000,
                      a_sigma = 2, b_sigma = 1, s_tau = 1) {
  n <- nrow(X); p <- ncol(X)
  XtX <- crossprod(X)
  Xty <- drop(crossprod(X, y))

  beta <- rep(0, p); sigma2 <- var(y); tau2 <- 1
  lambda2 <- rep(1, p); nu <- rep(1, p); xi <- 1

  keep <- n_iter - n_burn
  B <- matrix(NA, keep, p)
  S2 <- T2 <- A0 <- numeric(keep)

  for (it in 1:n_iter) {
    alpha <- rnorm(1, mean(y - X %*% beta), sqrt(sigma2 / n))

    A <- XtX
    diag(A) <- diag(A) + 1 / (tau2 * lambda2)
    U <- chol(A)
    m <- backsolve(U, forwardsolve(t(U), Xty - alpha * colSums(X)))
    beta <- m + sqrt(sigma2) * backsolve(U, rnorm(p))

    r <- y - alpha - X %*% beta
    sigma2 <- rinvgamma(1, a_sigma + (n + p) / 2,
                        b_sigma + (sum(r^2) + sum(beta^2 / (tau2 * lambda2))) / 2)
    lambda2 <- rinvgamma(p, 1, 1 / nu + beta^2 / (2 * sigma2 * tau2))
    tau2 <- rinvgamma(1, (p + 1) / 2, 1 / xi + sum(beta^2 / lambda2) / (2 * sigma2))
    nu <- rinvgamma(p, 1, 1 + 1 / lambda2)
    xi <- rinvgamma(1, 1, 1 / s_tau^2 + 1 / tau2)

    if (it > n_burn) {
      k <- it - n_burn
      B[k, ] <- beta; S2[k] <- sigma2; T2[k] <- tau2; A0[k] <- alpha
    }
  }
  list(beta = B, sigma2 = S2, tau2 = T2, alpha = A0)
}

fit <- horseshoe(y, X, s_tau = 5 / 45 / sqrt(n))   # tau scale for ~5 signals

# ---- Metrics ----------------------------------------------------------------
beta_hat <- colMeans(fit$beta)
ci <- apply(fit$beta, 2, quantile, c(0.025, 0.975))
selected <- ci[1, ] > 0 | ci[2, ] < 0
ols <- coef(lm(y ~ X))[-1]

sink("results/metrics.txt", split = TRUE)
cat("MSE horseshoe:", mean((beta_hat - beta_true)^2), "\n")
cat("MSE OLS:      ", mean((ols - beta_true)^2), "\n")
cat("95% CI coverage:", mean(ci[1, ] <= beta_true & beta_true <= ci[2, ]), "\n")
cat("Selected (CI excludes 0):", which(selected), "\n")
cat("True signals:            ", which(beta_true != 0), "\n")
cat("Intercept:", mean(fit$alpha), "  sigma:", mean(sqrt(fit$sigma2)), "\n")

print(round(data.frame(true = beta_true, horseshoe = beta_hat, ols = ols,
                       lower = ci[1, ], upper = ci[2, ])[beta_true != 0, ], 3))
sink()

# ---- Plots ------------------------------------------------------------------
# Each plot is saved as a PNG in results/ and also shown on screen in R.
save_plot <- function(file, draw) {
  png(file.path("results", file), width = 1400, height = 900, res = 150)
  draw()
  dev.off()
  if (interactive()) draw()
}

# 1. Estimates vs truth
save_plot("estimates.png", function() {
  plot(1:p, beta_hat, pch = 19, ylim = range(ci, ols),
       xlab = "j", ylab = "beta", main = "Horseshoe estimates and 95% CIs")
  segments(1:p, ci[1, ], 1:p, ci[2, ])
  points(1:p, beta_true, pch = 4, col = "red", cex = 1.5)
  points(1:p, ols, col = "grey")
  abline(h = 0, lty = 2)
  legend("topright", c("horseshoe", "truth", "OLS"),
         pch = c(19, 4, 1), col = c("black", "red", "grey"))
})

# 2. Trace plots
save_plot("traces.png", function() {
  par(mfrow = c(2, 2))
  plot(fit$beta[, 1], type = "l", main = "beta 1 (true 3)", ylab = "")
  plot(fit$beta[, 35], type = "l", main = "beta 35 (true -0.5)", ylab = "")
  plot(log(sqrt(fit$tau2)), type = "l", main = "log tau", ylab = "")
  plot(sqrt(fit$sigma2), type = "l", main = "sigma", ylab = "")
  par(mfrow = c(1, 1))
})

# 3. Posterior densities: a signal vs a null coefficient
save_plot("densities.png", function() {
  plot(density(fit$beta[, 20]), xlim = c(-0.5, 1.5),
       main = "Posterior: beta 20 (signal) vs beta 2 (null)")
  lines(density(fit$beta[, 2]), col = "blue")
  abline(v = c(0, 1), col = c("blue", "black"), lty = 2)
})
