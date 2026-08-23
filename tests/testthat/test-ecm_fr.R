test_that(
  "ecm_fr() has finite output for reasonable inputs",
  {
    data(Scenario1_MSFR)
    cnst <- list("block_lower1", "block_lower2", "null")
    mtd  <- "adhoc"
    for ( c in cnst ) {
      init <- start_fr(X_s, B_s, k, j_s, constraint = c, method = mtd)
      out  <- ecm_fr(X_s, B_s, init)
      expect_true(all(sapply(out$Psi_s, is.finite)))
      expect_true(all(is.finite(out$Phi)))
    }

    # mtd <- "fa"
    # init <- start_fr(X_s, B_s, k, j_s, method = mtd)
    # out  <- ecm_fr(X_s, B_s, init)
    # expect_true(all(sapply(out$Psi_s, is.finite)))
    # expect_true(all(is.finite(out$Phi)))
  }
)
