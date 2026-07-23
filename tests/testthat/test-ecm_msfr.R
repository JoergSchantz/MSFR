# Tests for R/ecm_msfr.R (the MSFR model: common Phi + study-specific Lambda_s +
# covariate beta), covering the three checks TODO.md's P0 section called for once the
# .update_Lambda()/WoodburyMatrix::/step-halving fixes were applied:
#   1. brute-force numerical cross-check of the E-step (.exp_values())
#   2. a simulated multi-study recovery test for ecm_msfr()
#   3. confirmation that Phi/Lambda_s/beta/Psi_s are all genuinely updated
#      iteration-to-iteration (the bug class that hit ecm_fr()'s `Phi` earlier)

simulate_msfr_data <- function(seed, p = 10, S = 2, k = 2, j_s = c(1, 1),
                                p_b = 2, n_s = c(300, 300)) {
  set.seed(seed)
  true_Phi <- matrix(0, p, k)
  true_Phi[lower.tri(true_Phi, diag = TRUE)] <-
    runif(sum(lower.tri(true_Phi, diag = TRUE)), 0.4, 0.9)
  true_Lambda <- lapply(j_s, function(j) matrix(runif(p * j, 0.3, 0.8), p, j))
  true_beta <- matrix(runif(p * p_b, -1, 1), p, p_b)
  true_psi <- lapply(1:S, function(s) runif(p, 0.3, 0.8))

  B_s <- lapply(n_s, function(n) matrix(rnorm(n * p_b), n, p_b))
  X_s <- Map(function(B, Lam, psi, n) {
    Fc <- matrix(rnorm(n * k), n, k)
    Fs <- matrix(rnorm(n * ncol(Lam)), n, ncol(Lam))
    E  <- matrix(rnorm(n * p), n, p) %*% diag(sqrt(psi))
    B %*% t(true_beta) + Fc %*% t(true_Phi) + Fs %*% t(Lam) + E
  }, B_s, true_Lambda, true_psi, n_s)

  list(X_s = X_s, B_s = B_s, p = p, S = S, k = k, j_s = j_s, p_b = p_b, n_s = n_s,
       true_Phi = true_Phi, true_Lambda = true_Lambda, true_beta = true_beta,
       true_psi = true_psi)
}

test_that(
  ".exp_values() (MSFR E-step) matches a brute-force reference", 
  {
  d <- simulate_msfr_data(seed = 2024)
  Phi <- matrix(rnorm(d$p * d$k, sd = 0.5), d$p, d$k)
  Lambda_s <- lapply(d$j_s, function(j) matrix(rnorm(d$p * j, sd = 0.5), d$p, j))
  psi_vals <- lapply(1:d$S, function(s) runif(d$p, 0.3, 1.2))
  Psi_s <- lapply(psi_vals, diag)
  cov_s <- lapply(1:d$S, function(s) {
    M <- matrix(rnorm(d$p * d$p), d$p, d$p)
    crossprod(M) / d$p + diag(d$p)
  })
  X_s_tilde <- d$X_s

  ref <- lapply(1:d$S, function(s) {
    L <- Lambda_s[[s]]; Psi <- Psi_s[[s]]; C <- cov_s[[s]]
    Sig <- Phi %*% t(Phi) + L %*% t(L) + Psi
    Sig1 <- solve(Sig)
    delta_Phi <- t(Phi) %*% Sig1
    delta_L   <- t(L) %*% Sig1
    I_k <- diag(d$k); I_j <- diag(ncol(L))
    list(
      exp_xf = C %*% t(delta_Phi),
      exp_xl = C %*% t(delta_L),
      exp_ff = delta_Phi %*% C %*% t(delta_Phi) + (I_k - delta_Phi %*% Phi),
      exp_ll = delta_L %*% C %*% t(delta_L) + (I_j - delta_L %*% L),
      exp_fl = delta_Phi %*% C %*% t(delta_L) - delta_Phi %*% L,
      exp_f  = delta_Phi %*% t(X_s_tilde[[s]]),
      exp_l  = delta_L %*% t(X_s_tilde[[s]])
    )
  })

  out1 <- .exp_values(Phi, Lambda_s, Psi_s, CM_step = 1, cov_s = cov_s)
  out4 <- .exp_values(Phi, Lambda_s, Psi_s, CM_step = 4, X_s_tilde = X_s_tilde)

  for (s in 1:d$S) {
    expect_equal(out1$exp_xf[[s]], ref[[s]]$exp_xf, tolerance = 1e-8)
    expect_equal(out1$exp_xl[[s]], ref[[s]]$exp_xl, tolerance = 1e-8)
    expect_equal(out1$exp_ff[[s]], ref[[s]]$exp_ff, tolerance = 1e-8)
    expect_equal(out1$exp_ll[[s]], ref[[s]]$exp_ll, tolerance = 1e-8)
    expect_equal(out1$exp_fl[[s]], ref[[s]]$exp_fl, tolerance = 1e-8)
    expect_equal(out4$exp_f[[s]], ref[[s]]$exp_f, tolerance = 1e-8)
    expect_equal(out4$exp_l[[s]], ref[[s]]$exp_l, tolerance = 1e-8)
  }
})

test_that(
  "ecm_msfr() recovers simulated Phi/Lambda_s/beta/Psi_s reasonably well", 
  {
  d <- simulate_msfr_data(seed = 5)
  start <- start_msfa(d$X_s, d$B_s, p_b = d$p_b, k = d$k, j_s = d$j_s,
                       constraint = "block_lower2")

  res <- ecm_msfr(d$X_s, d$B_s, start, nIt = 500, tol = 1e-7, trace = FALSE)

  expect_true(is.finite(res$loglik))
  expect_true(all(is.finite(res$Phi)))
  expect_true(all(is.finite(res$beta)))
  expect_true(all(sapply(res$psi_s, function(x) all(is.finite(x)))))

  # beta is identified (unlike Phi/Lambda_s, which are only identified up to
  # rotation/sign under this constraint), so it's the most direct recovery check
  expect_lt(max(abs(res$beta - d$true_beta)), 0.3)

  # reconstructed covariance should track the (covariate-adjusted) sample covariance
  for (s in 1:d$S) {
    X_tilde_s <- d$X_s[[s]] - d$B_s[[s]] %*% t(res$beta)
    Sigma_hat <- res$Phi %*% t(res$Phi) + res$Lambda_s[[s]] %*% t(res$Lambda_s[[s]]) +
      diag(res$psi_s[[s]])
    expect_lt(max(abs(Sigma_hat - cov(X_tilde_s))), 0.3)
  }
})

test_that(
  "ecm_msfr() actually updates Phi/Lambda_s/beta/Psi_s across iterations",
  {
  # Regression test for the bug class already found in ecm_fr(): a parameter that's
  # computed as *_new every iteration but never carried forward into the variable the
  # next iteration's E-step actually reads, so it silently never changes.
  d <- simulate_msfr_data(seed = 7)
  start <- start_msfa(d$X_s, d$B_s, p_b = d$p_b, k = d$k, j_s = d$j_s,
                       constraint = "block_lower2")

  res <- ecm_msfr(d$X_s, d$B_s, start, nIt = 25, tol = 1e-12, trace = FALSE)

  expect_gt(res$iter, 1)  # sanity: make sure it didn't stop after the 1st iteration
  expect_gt(max(abs(res$Phi - start$Phi)), 1e-6)
  expect_gt(max(abs(res$beta - start$beta)), 1e-6)
  for (s in 1:d$S) {
    expect_gt(max(abs(res$Lambda_s[[s]] - start$Lambda_s[[s]])), 1e-6)
    expect_gt(max(abs(res$psi_s[[s]] - start$psi_s[[s]])), 1e-6)
  }
})
