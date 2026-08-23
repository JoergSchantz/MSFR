# TESTS ----
test_that(
  "ecm_msfr demonstrates reasonable behavior for unintended inputs",
  {
    expect_error(ecm_msfr(matrix(), matrix(), numeric(), numeric()))
    expect_error(ecm_msfr(NaN, NaN, NaN, NaN))
    expect_error(ecm_msfr(matrix(0), matrix(0), 1, 1))
  }
)

test_that(
  "ecm_msfr has finite outputs for reasonable inputs",
  {
    con <- list("block_lower1", "block_lower2")
    for (c in con) {
      data(Scenario1_MSFR)
        
      start <- start_msfa(X_s, B_s, k, j_s, constraint = c)
      res <- ecm_msfr(X_s, B_s, start, nIt = 100, constraint = c)

      expect_true(all(is.finite(res$Phi)))
      expect_true(all(is.finite(res$Beta)))
      expect_true(all(sapply(res$Lambda_s, is.finite)))
      expect_true(all(sapply(res$Psi_s, is.finite)))
    }
  }
)
