# regular unit test for `start_msfa()`
test_that(
  "param initialisation behaves as intended for various inputs",
  {
    expect_error(start_msfa(NaN, NaN, NaN, NaN))
    expect_error(start_msfa(matrix(), matrix(), numeric(),numeric()))
    expect_error(start_msfa(matrix(), matrix(), character(),character()))
  }
)
test_that(
  "param initalisation of betas works for Scenario1_MSFR dataset",
  {
    data(Scenario1_MSFR)
    init <- start_msfa(X_s, B_s, k, j_s)
    beta_fm1 <- init$beta

    B <- Reduce(rbind, B_s)
    X <- Reduce(rbind, X_s)
    beta_mle = t(solve(t(B) %*% (B)) %*% t(B) %*% X)

    expect_true(all(abs(beta_fm1 - beta_mle) < 1e-10))
  }
)

test_that(
  "constraint blocks work",
  {
    data(Scenario1_MSFR)
    init_bl2  <- start_msfa(X_s, B_s, k ,j_s, constraint = "block_lower2")
    init_bl1  <- start_msfa(X_s, B_s, k ,j_s, constraint = "block_lower1")
    init_null <- start_msfa(X_s, B_s, k ,j_s, constraint = "null")
    
    for(params in list(init_bl2, init_bl1, init_null)){
      Phi    <- params$Phi
      Lambda <- params$Lambda
      psi    <- params$psi
      beta   <- params$beta
      # test for sensible values
      expect_true(all(Phi < Inf))
      expect_true(all(vapply(Lambda, function(L) all(L < Inf), logical(1))))
      expect_true(all(vapply(psi, function(p) all(p < Inf), logical(1))))
      expect_true(all(beta < Inf))
            
      expect_true(all(Phi > -Inf))
      expect_true(all(vapply(Lambda, function(L) all(L > -Inf), logical(1))))
      expect_true(all(vapply(psi, function(p) all(p > -Inf), logical(1))))
      expect_true(all(beta > -Inf))
    }
    
  }
)


test_that(
  "method 'fa' works",
  {
    data(Scenario1_MSFR)
    init <- start_msfa(X_s, B_s, k ,j_s, method = "fa")
    Phi    <- init$Phi
    Lambda <- init$Lambda
    psi    <- init$psi
    beta   <- init$beta

    # test for sensible values
    expect_true(all(Phi < Inf))
    expect_true(all(vapply(Lambda, function(L) all(L < Inf), logical(1))))
    expect_true(all(vapply(psi, function(p) all(p < Inf), logical(1))))
    expect_true(all(beta < Inf))
          
    expect_true(all(Phi > -Inf))
    expect_true(all(vapply(Lambda, function(L) all(L > -Inf), logical(1))))
    expect_true(all(vapply(psi, function(p) all(p > -Inf), logical(1))))
    expect_true(all(beta > -Inf))
  }
)
