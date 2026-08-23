test_that(
  "vcov_msfa does its job lol",
  {
    data(Scenario1_MSFR)
    init <- start_msfa(X_s, B_s, k, j_s)
    mle <- ecm_msfr(X_s, B_s, init, tol = 10^-7)
    expect_no_error(vcov_msfa(X_s, mle, getgrad = T))
    expect_no_warning(vcov_msfa(X_s, mle, getgrad = T))
  }
)
