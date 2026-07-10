## =============================================================================
## .exp_values_fr()  --  EM-hot-loop optimized version
##
## This assumes the earlier Woodbury refactor is correct and asks: given this
## function now runs once per EM iteration (possibly thousands of times over
## a run, S times each), what's left on the table?
##
## Three findings, in order of actual impact:
##
## (1) ALGORITHMIC: the previous version computed a Woodbury block from
##     Psi_s1 (for Sig_s1 / delta_Phi / Delta_Phi) AND a second, separate
##     Woodbury block from 1/Psi_s (for Woodbury_f / E_fis_x_is). Since
##     Psi_s1 IS Psi_s^{-1} (confirmed), those two blocks are, mathematically,
##     the exact same computation performed twice. This was already true of
##     the *original* code (`Sig_s1` used Psi_s1 directly; `Woodbury_f` used
##     `.inv_Psi(Psi_s[[s]])` to recompute the same inverse from scratch).
##     Merging them removes an entire redundant O(n k^2) block per group,
##     per call -- this dominates everything else below.
##
##     A consequence: Psi_s is no longer touched anywhere in this function.
##     det(Psi_s) is obtained from Psi_s1 via -sum(log(psi_s1)) instead.
##     `Psi_s` is kept as a parameter purely so existing call sites (which
##     pass it positionally) don't break -- because R evaluates arguments
##     lazily, an unused parameter that the body never references is never
##     forced, so keeping it costs ~nothing here even though it's unused.
##     (If you're willing to touch every call site, dropping it outright
##     also saves whatever `[[` your EM driver used to build the list.)
##
## (2) CONSTANT-FACTOR / BLAS: t(X) %*% Y patterns were replaced with
##     crossprod()/tcrossprod(), which do the transpose and multiply in one
##     fused BLAS call instead of materialising a transposed copy first.
##     The k x k inverse now uses chol2inv(chol(.)) instead of solve(),
##     which is faster for symmetric PD matrices (with a symmetrization
##     step + solve() fallback, since this now runs unattended for
##     thousands of iterations and a rare float-asymmetry chol() failure
##     would be worse than a slightly slower solve()).
##
## (3) ALLOCATION: get_diag_vec / woodbury_block / apply_Sig_s1 used to be
##     closures defined *inside* .exp_values_fr, so R allocated fresh
##     function objects (and closure environments) on every single call.
##     They're now module-level functions, created once. Sig_s1[[s]] no
##     longer stores a redundant copy of `Phi` or a per-group reference to
##     `apply_Sig_s1` (same object, S times over, every call) -- it stores
##     only what's actually group-specific. `dim(Phi)` is read once instead
##     of twice. `I_k` can be built once outside the EM loop (it depends
##     only on k, which is fixed for the run) and handed in, instead of
##     being reallocated on every call.
##
##     Note what's deliberately *not* touched: the `for (s in seq_len(S))`
##     loop itself. S is small by assumption, so R's per-iteration loop
##     overhead is negligible next to the O(n k^2) linear algebra inside
##     each iteration -- "vectorizing across groups" (e.g. block-diagonal
##     stacking) would add real complexity for no measurable win, and
##     could even hurt by inflating the size of intermediate matrices.
##     The vectorization that matters here is *within* each iteration,
##     handed off to BLAS via crossprod/tcrossprod/chol2inv.
##
##     Also new: `getdet` now always produces `logds_s` (the numerically
##     safe form) but only exponentiates to `ds_s` if `raw_det = TRUE`,
##     since an EM driver tracking log-likelihood has no use for S
##     `exp()` calls it's going to discard every iteration.
## =============================================================================

## ---- module-level helpers: created once, not re-created on every call ----

.get_diag_vec <- function(Psi, dense_warn_n = 3000L) {
  if (is.matrix(Psi)) {
    if (nrow(Psi) > dense_warn_n) {
      warning("A ", nrow(Psi), " x ", nrow(Psi), " dense matrix was passed ",
              "for a diagonal covariance; store it as a plain vector ",
              "upstream to avoid O(n^2) memory and repeated diag() ",
              "extraction on every EM iteration.", call. = FALSE)
    }
    diag(Psi)
  } else {
    Psi
  }
}

# One Woodbury block for Sig = diag(1/prec_vec) + Phi %*% t(Phi).
.woodbury_block <- function(Phi, prec_vec, I_k) {
  Y <- Phi * prec_vec                  # n x k,  == diag(prec_vec) %*% Phi
  A <- crossprod(Phi, Y)               # k x k,  == t(Phi) %*% diag(prec_vec) %*% Phi  (fused, no explicit t())
  A <- (A + t(A)) * 0.5                # guard against float asymmetry before chol()
  M <- tryCatch(
    chol2inv(chol(I_k + A)),
    error = function(e) solve(I_k + A) # robust fallback; should essentially never trigger
  )
  W <- tcrossprod(M, Y)                # k x n,  == M %*% t(Y)  (fused, no explicit t())
  list(A = A, M = M, Y = Y, W = W)
}

# Apply the *implicit* Sig_s^{-1} = diag(prec) - Y M t(Y) to an n x p matrix
# X without ever forming it densely: O(n k p) instead of O(n^2 p).
.apply_Sig_s1 <- function(compact, X) {
  YtX <- crossprod(compact$Y, X)                        # k x p, == t(Y) %*% X
  compact$diag * X - compact$Y %*% (compact$M %*% YtX)  # == diag(prec) %*% X  -  Y M t(Y) X
}

## ---- main function --------------------------------------------------------

.exp_values_fr <- function(Phi, Psi_s, Psi_s1, cov_s, X_s_tilde,
                            getdet = FALSE,
                            raw_det = FALSE,
                            dense_Sig_s1 = FALSE,
                            dense_warn_n = 3000L,
                            I_k = NULL)
{
  # `Psi_s` intentionally unused -- see note (1) above. Kept as a parameter
  # only so existing positional call sites keep working; R's lazy argument
  # evaluation means it's never forced, so this costs nothing here.

  d <- dim(Phi); n <- d[1L]; k <- d[2L]
  if (is.null(I_k)) I_k <- diag(1, k)   # pass this in from the EM driver to
                                         # avoid rebuilding it every iteration
  S <- length(Psi_s1)

  Sig_s1     <- vector("list", S)
  ds_s       <- vector("list", S)
  logds_s    <- vector("list", S)
  Txsfcs     <- vector("list", S)
  Tfcsfcs    <- vector("list", S)
  E_fis_x_is <- vector("list", S)

  for (s in seq_len(S)) {

    psi_s1 <- .get_diag_vec(Psi_s1[[s]], dense_warn_n)   # Psi_s^{-1}, length n

    ## single Woodbury block. Under the confirmed Psi_s1 == Psi_s^{-1},
    ## this replaces what used to be two identical blocks: it feeds
    ## Sig_s1 / delta_Phi / Delta_Phi *and* the old separate
    ## Woodbury_f / E_fis_x_is computation, because they were always the
    ## same object.
    wb <- .woodbury_block(Phi, psi_s1, I_k)
    delta_Phi_s <- wb$W    # t(Phi) %*% Sig_s^{-1}                 (k x n)
    Delta_Phi_s <- wb$M    # I_k - t(Phi) %*% Sig_s^{-1} %*% Phi   (k x k)

    Sig_s1[[s]] <- list(diag = psi_s1, Y = wb$Y, M = wb$M)   # apply via .apply_Sig_s1()

    if (dense_Sig_s1) {
      if (n > dense_warn_n) {
        warning("Densifying a ", n, " x ", n, " matrix for group ", s,
                 "; this defeats the memory savings of the Woodbury ",
                 "representation and should not be enabled inside an EM ",
                 "loop.", call. = FALSE)
      }
      Sig_s1[[s]]$dense <- diag(psi_s1) - tcrossprod(wb$Y %*% wb$M, wb$Y)
    }

    if (getdet) {
      # det(Psi_s) = 1 / det(Psi_s1) = 1 / prod(psi_s1); on the log scale,
      # so this stays finite for large n regardless of raw_det.
      logdet_Ik_A  <- as.numeric(determinant(I_k + wb$A, logarithm = TRUE)$modulus)
      logds_s[[s]] <- -sum(log(psi_s1)) + logdet_Ik_A
      if (raw_det) ds_s[[s]] <- exp(logds_s[[s]])
    }

    Txsfcs[[s]]     <- tcrossprod(cov_s[[s]], delta_Phi_s)          # n x k
    Tfcsfcs[[s]]    <- delta_Phi_s %*% Txsfcs[[s]] + Delta_Phi_s    # k x k, Txsfcs reused
    E_fis_x_is[[s]] <- tcrossprod(delta_Phi_s, X_s_tilde[[s]])      # k x N_s
  }

  list(Txsfcs = Txsfcs, Tfcsfcs = Tfcsfcs, E_fis_x_is = E_fis_x_is,
       ds_s = ds_s, logds_s = logds_s, Sig_s1 = Sig_s1)
}
