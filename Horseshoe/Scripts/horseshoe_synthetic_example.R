# Horseshoe linear regression on synthetic data: a worked example.
#
# Implements the main sampler of ../Notes/horseshoe_linear_regression.tex:
#   y_i = alpha + x_i' beta + eps_i,      eps_i ~ N(0, sigma^2)
#   beta_j | sigma^2, tau^2, lambda_j^2 ~ N(0, sigma^2 tau^2 lambda_j^2)
#   lambda_j ~ HC(0, 1),  tau ~ HC(0, s_tau),  sigma^2 ~ IG(a_sigma, b_sigma)
#   p(alpha) propto 1  (unpenalized intercept)
# using the Makalic-Schmidt inverse-gamma augmentation and the ordered Gibbs
# sweep of the notes (intercept, slopes, sigma^2, lambda^2, tau^2, nu, xi).
#
# This is an illustrative example, not a validated production sampler.
# Base R only. Run from any directory:
#   Rscript horseshoe_synthetic_example.R
# Outputs (draws, summaries, plot, session info) go to ./output/<run id>/
# next to this script.

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

seed       <- 20261020L
n          <- 100L    # observations
p          <- 50L     # penalized slopes (intercept not counted)
rho        <- 0.5     # AR(1) correlation between neighbouring raw predictors
alpha_true <- 2
sigma_true <- 1

# True slopes on the standardized-predictor scale: a few strong, a few weak,
# the rest exactly zero.
beta_true <- numeric(p)
beta_true[c(1, 5, 10, 20, 35)] <- c(3, -2, 1.5, 1, -0.5)

# Prior hyperparameters. IG(2, 1) for sigma^2 has mean 1 and infinite variance
# (an illustrative choice, see the notes' noise-prior section). s_tau follows
# the sparsity calibration s_tau = p0 / (p - p0) / sqrt(n); the guess p0 = 5
# happens to equal the true number of nonzero slopes here. The reference
# choice s_tau = 1 is a valid alternative.
a_sigma  <- 2
b_sigma  <- 1
p0_guess <- 5
s_tau    <- p0_guess / (p - p0_guess) / sqrt(n)

# MCMC settings
n_chains <- 4L
n_iter   <- 22000L   # tau mixes slowly; shorter runs gave R-hat(log tau) > 1.01
n_burn   <- 2000L

# Decision settings (Rules 1-3 of the notes)
cred_level <- 0.95
delta      <- 0.25   # practical-relevance threshold, response units per SD
q_cutoff   <- 0.5    # equal false-positive / false-negative costs
r_cutoff   <- 0.5    # conventional shrinkage-score cutoff

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Z ~ IG(a, b) with density b^a / Gamma(a) z^(-a-1) exp(-b/z); vectorized in b.
rinvgamma <- function(k, a, b) b / rgamma(k, shape = a, rate = 1)

# Rank-normalized split R-hat (Vehtari et al. 2021) for an iterations x chains
# matrix. Returns the maximum of the bulk and folded (tail) versions.
split_rhat <- function(draws) {
  rhat_basic <- function(m) {
    half <- floor(nrow(m) / 2)
    m <- cbind(m[seq_len(half), , drop = FALSE],
               m[nrow(m) - half + seq_len(half), , drop = FALSE])
    s <- nrow(m)
    b <- s * var(colMeans(m))
    w <- mean(apply(m, 2, var))
    sqrt(((s - 1) / s * w + b / s) / w)
  }
  z_scale <- function(m) {
    r <- rank(m, ties.method = "average")
    matrix(qnorm((r - 3 / 8) / (length(m) + 1 / 4)), nrow(m))
  }
  bulk <- rhat_basic(z_scale(draws))
  tail <- rhat_basic(z_scale(abs(draws - median(draws))))
  max(bulk, tail)
}

# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------

set.seed(seed)

Sigma_x <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
X_raw   <- matrix(rnorm(n * p), n, p) %*% chol(Sigma_x)

# Standardize with the 1/n convention of the notes, so 1'X = 0 and X_j'X_j = n.
x_mean <- colMeans(X_raw)
x_sd   <- sqrt(colMeans(sweep(X_raw, 2, x_mean)^2))
X      <- sweep(sweep(X_raw, 2, x_mean), 2, x_sd, "/")

y <- alpha_true + drop(X %*% beta_true) + rnorm(n, sd = sigma_true)

stopifnot(max(abs(colSums(X))) < 1e-8,
          max(abs(colSums(X^2) - n)) < 1e-8)

# ---------------------------------------------------------------------------
# Gibbs sampler
# ---------------------------------------------------------------------------

horseshoe_gibbs <- function(y, X, n_iter, n_burn, a_sigma, b_sigma, s_tau,
                            init) {
  n <- nrow(X); p <- ncol(X)
  XtX <- crossprod(X)
  Xty <- drop(crossprod(X, y))
  Xt1 <- colSums(X)

  alpha <- init$alpha; beta <- init$beta; sigma2 <- init$sigma2
  lambda2 <- init$lambda2; tau2 <- init$tau2; nu <- init$nu; xi <- init$xi

  n_keep <- n_iter - n_burn
  out <- list(alpha   = numeric(n_keep),
              beta    = matrix(NA_real_, n_keep, p),
              sigma2  = numeric(n_keep),
              tau2    = numeric(n_keep),
              lambda2 = matrix(NA_real_, n_keep, p))

  for (it in seq_len(n_iter)) {
    # 1. Intercept: N(mean(y - X beta), sigma^2 / n)
    alpha <- rnorm(1, mean(y - X %*% beta), sqrt(sigma2 / n))

    # 2. Slopes as one block: N(m, sigma^2 A^{-1}), A = X'X + D^{-1}.
    #    chol() returns upper U with A = U'U; m = A^{-1} b via two triangular
    #    solves, and solving U v = z gives Cov(v) = A^{-1}.
    d <- tau2 * lambda2
    A <- XtX
    diag(A) <- diag(A) + 1 / d
    U <- chol(A)
    b <- Xty - alpha * Xt1
    m <- backsolve(U, forwardsolve(t(U), b))
    beta <- m + sqrt(sigma2) * backsolve(U, rnorm(p))

    # 3. Noise variance: IG(a + (n + p)/2, b + (R + Q)/2)
    resid  <- y - alpha - drop(X %*% beta)
    R      <- sum(resid^2)
    Q      <- sum(beta^2 / d)
    sigma2 <- rinvgamma(1, a_sigma + (n + p) / 2, b_sigma + (R + Q) / 2)

    # 4. Local variances: IG(1, 1/nu_j + beta_j^2 / (2 sigma^2 tau^2))
    lambda2 <- rinvgamma(p, 1, 1 / nu + beta^2 / (2 * sigma2 * tau2))

    # 5. Global variance: IG((p + 1)/2, 1/xi + sum(beta^2 / lambda2) / (2 sigma^2))
    tau2 <- rinvgamma(1, (p + 1) / 2,
                      1 / xi + sum(beta^2 / lambda2) / (2 * sigma2))

    # 6. Auxiliary variables
    nu <- rinvgamma(p, 1, 1 + 1 / lambda2)
    xi <- rinvgamma(1, 1, 1 / s_tau^2 + 1 / tau2)

    if (!all(is.finite(c(alpha, beta, sigma2, lambda2, tau2, nu, xi))) ||
        any(c(sigma2, lambda2, tau2, nu, xi) <= 0)) {
      stop("Non-finite or non-positive state at iteration ", it)
    }

    # 7. Save the completed state after burn-in
    if (it > n_burn) {
      k <- it - n_burn
      out$alpha[k]      <- alpha
      out$beta[k, ]     <- beta
      out$sigma2[k]     <- sigma2
      out$tau2[k]       <- tau2
      out$lambda2[k, ]  <- lambda2
    }
  }
  out
}

# Dispersed starting values, one set per chain.
chain_seeds <- seed + seq_len(n_chains)
t_start <- Sys.time()
chains <- lapply(seq_len(n_chains), function(ch) {
  set.seed(chain_seeds[ch])
  init <- list(alpha   = mean(y) + rnorm(1),
               beta    = rnorm(p, sd = 0.1),
               sigma2  = var(y) * exp(rnorm(1)),
               lambda2 = rep(1, p),
               tau2    = s_tau^2 * exp(rnorm(1)),
               nu      = rep(1, p),
               xi      = 1)
  horseshoe_gibbs(y, X, n_iter, n_burn, a_sigma, b_sigma, s_tau, init)
})
t_end <- Sys.time()

# ---------------------------------------------------------------------------
# Posterior summaries
# ---------------------------------------------------------------------------

pool <- function(name) do.call(rbind, lapply(chains, function(ch) {
  v <- ch[[name]]
  if (is.matrix(v)) v else matrix(v)
}))
beta_draws    <- pool("beta")
lambda2_draws <- pool("lambda2")
sigma2_draws  <- drop(pool("sigma2"))
tau2_draws    <- drop(pool("tau2"))
alpha_draws   <- drop(pool("alpha"))

probs <- c((1 - cred_level) / 2, 1 - (1 - cred_level) / 2)
beta_ci <- t(apply(beta_draws, 2, quantile, probs = probs))

# Shrinkage factor kappa_j = 1 / (1 + c_j tau^2 lambda_j^2) with c_j = n.
kappa_draws <- 1 / (1 + n * tau2_draws * lambda2_draws)

# Convergence screen: R-hat for each slope, log tau, log sigma, intercept.
per_chain <- function(f) sapply(chains, f)
rhat <- c(
  sapply(seq_len(p), function(j) split_rhat(per_chain(function(ch) ch$beta[, j]))),
  alpha     = split_rhat(per_chain(function(ch) ch$alpha)),
  log_sigma = split_rhat(per_chain(function(ch) 0.5 * log(ch$sigma2))),
  log_tau   = split_rhat(per_chain(function(ch) 0.5 * log(ch$tau2)))
)
names(rhat)[seq_len(p)] <- paste0("beta", seq_len(p))

summary_tab <- data.frame(
  j         = seq_len(p),
  true      = beta_true,
  post_mean = colMeans(beta_draws),
  lower     = beta_ci[, 1],
  upper     = beta_ci[, 2],
  q_rel     = colMeans(abs(beta_draws) > delta),   # Rule 2
  r_shrink  = colMeans(1 - kappa_draws),           # Rule 3
  rhat      = rhat[seq_len(p)]
)
summary_tab$sel_ci    <- summary_tab$lower > 0 | summary_tab$upper < 0
summary_tab$sel_q     <- summary_tab$q_rel > q_cutoff
summary_tab$sel_r     <- summary_tab$r_shrink > r_cutoff
summary_tab$covered   <- summary_tab$lower <= beta_true &
                         beta_true <= summary_tab$upper

ols <- coef(lm(y ~ X))[-1]
truly_nonzero <- beta_true != 0

decision_table <- function(sel) {
  c(true_pos  = sum(sel & truly_nonzero),
    false_pos = sum(sel & !truly_nonzero),
    false_neg = sum(!sel & truly_nonzero))
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

cat("Horseshoe regression, synthetic data\n")
cat(sprintf("n = %d, p = %d, nonzero slopes = %d, rho = %.2f, seed = %d\n",
            n, p, sum(truly_nonzero), rho, seed))
cat(sprintf("Prior: sigma^2 ~ IG(%g, %g), tau ~ HC(0, %.5f)\n",
            a_sigma, b_sigma, s_tau))
cat(sprintf("MCMC: %d chains x %d iterations (%d burn-in), %.1f s elapsed\n\n",
            n_chains, n_iter, n_burn,
            as.numeric(difftime(t_end, t_start, units = "secs"))))

cat("Nonzero true slopes and largest posterior means:\n")
show_j <- union(which(truly_nonzero),
                order(-abs(summary_tab$post_mean))[1:8])
print(cbind(summary_tab[show_j, c("j", "true", "post_mean", "lower", "upper",
                                  "q_rel", "r_shrink", "rhat")],
            ols = ols[show_j]),
      digits = 3, row.names = FALSE)

cat(sprintf("\nIntercept: mean %.3f  (true %.3f)\n",
            mean(alpha_draws), alpha_true))
cat(sprintf("sigma:     mean %.3f  (true %.3f)\n",
            mean(sqrt(sigma2_draws)), sigma_true))
cat(sprintf("tau:       median %.5f\n", median(sqrt(tau2_draws))))

cat(sprintf("\nEstimation error, sum of squared errors over all slopes:\n"))
cat(sprintf("  horseshoe posterior mean: %.3f\n",
            sum((summary_tab$post_mean - beta_true)^2)))
cat(sprintf("  OLS:                      %.3f\n", sum((ols - beta_true)^2)))
cat(sprintf("%.0f%% interval coverage of true slopes: %d / %d\n",
            100 * cred_level, sum(summary_tab$covered), p))

cat("\nVariable-selection decisions (counts against the known truth):\n")
print(rbind(
  `Rule 1: CI excludes 0`              = decision_table(summary_tab$sel_ci),
  `Rule 2: P(|beta| > delta) > cutoff` = decision_table(summary_tab$sel_q),
  `Rule 3: E(1 - kappa) > cutoff`      = decision_table(summary_tab$sel_r)
))

cat(sprintf("\nMax R-hat over slopes, alpha, log sigma, log tau: %.3f (%s)\n",
            max(rhat), names(rhat)[which.max(rhat)]))
if (max(rhat) > 1.01) {
  cat("Warning: R-hat above 1.01; consider longer chains before use.\n")
}

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------

script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
if (is.na(script_dir)) script_dir <- getwd()
run_id  <- format(t_start, "%Y%m%dT%H%M%S")
out_dir <- file.path(script_dir, "output", paste0("run_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(list(settings = list(seed = seed, chain_seeds = chain_seeds, n = n,
                             p = p, rho = rho, alpha_true = alpha_true,
                             sigma_true = sigma_true, beta_true = beta_true,
                             a_sigma = a_sigma, b_sigma = b_sigma,
                             s_tau = s_tau, n_chains = n_chains,
                             n_iter = n_iter, n_burn = n_burn),
             data = list(X_raw = X_raw, X = X, y = y),
             chains = chains, summary = summary_tab, rhat = rhat),
        file.path(out_dir, "horseshoe_example.rds"))
write.csv(summary_tab, file.path(out_dir, "slope_summary.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "session_info.txt"))

pdf(file.path(out_dir, "horseshoe_example.pdf"), width = 9, height = 5)
plot(seq_len(p), summary_tab$post_mean, pch = 19, cex = 0.7,
     ylim = range(beta_ci, beta_true, ols),
     xlab = "Predictor j", ylab = expression(beta[j]),
     main = sprintf("Horseshoe posterior means and %.0f%% intervals",
                    100 * cred_level))
segments(seq_len(p), beta_ci[, 1], seq_len(p), beta_ci[, 2])
points(seq_len(p), beta_true, pch = 4, col = "red")
points(seq_len(p), ols, pch = 1, col = "grey50", cex = 0.7)
abline(h = 0, lty = 3)
legend("topright", c("posterior mean + interval", "truth", "OLS"),
       pch = c(19, 4, 1), col = c("black", "red", "grey50"), bty = "n")
invisible(dev.off())

cat("\nOutputs written to", out_dir, "\n")
