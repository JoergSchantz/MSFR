test_that(
  "heat_plot behaves es intended for ecm outputs",
  {
    data(Scenario1_MSFR)
    init <- start_msfa(X_s, B_s, k, j_s)
    mat <- ecm_msfr(X_s, B_s, init, tol = 10^-4)$Phi
    expect_no_error(heat_plot(mat))
    expect_no_warning(heat_plot(mat))
  }
)
