321qcvbnm,.-# MSFR package — roadmap to a releasable R package

This is a working audit of the whole package (everything except `R/vcov_msfa.R` and
`R/heat_plot.R`, per instruction) plus a mapping of the non-coding steps still open in
[`CRAN_checklist.md`](CRAN_checklist.md). Everything below was checked against the actual
code — either by reading it closely or by running it — not guessed. File:line references
point at the current state of the repo.

**How this is organized:** priority tiers, not file-by-file. P0 blocks everything else
(the package can't be trusted or tested until these are fixed); P1–P2 are what "the package
builds and has a safety net" requires; P3+ is what turns a working package into a
publishable/citable one. Checkboxes match the style of `CRAN_checklist.md` so the two files
can eventually be read together.

---

## P0 — Correctness blockers

These make specific functions wrong or unusable. Nothing else in this document is worth
doing until these are addressed, because right now roughly half the package's model-fitting
functions don't run at all.

### `ecm_msfr()` (`R/ecm_msfr.R`) — now runs end-to-end, verified by  tests

This is the file `file-structure.md` calls "the main function for the multi-study factor
regression algorithm." As of this pass, all four bu gs found in it are fixed and verified —
`tests/testthat/test-ecm_msfr.R` runs clean: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 28 ]`, via
`devtools::test()`, i.e. through the real package namespace, not ad hoc `source()`.

- [x] **`.update_Lambda()` (`R/ecm_msfr.R:38-46`).** Two rounds: the original bug was `p`
  referenced but never defined (this top-level function can't see its caller's local `p` via
  lexical scoping); the first fix attempt replaced it with `dim(Phi)[i]`, which swapped one
  undefined-variable bug for another (`i` isn't in scope here either — confirmed by
  `test-ecm_msfr.R` failing with `Error in ... : object 'i' not found` at the exact call
  site). Now fixed to `p <- dim(Phi)[1]` and reverified — the recovery and
  parameter-movement tests both get past CM3 and pass.
- [x] **`kronecker(WoodburyMatrix::t(exp_ff), nPsi1)` (`R/ecm_msfr.R:177`)** — the
  `WoodburyMatrix::` prefix (almost certainly a typo for plain `t()`) was loading that
  package's namespace and corrupting `kronecker()`'s S4 method dispatch for the very next
  call on plain matrices. Fixed to plain `t(exp_ff)`; `WoodburyMatrix` was never declared as
  a dependency anywhere, so this only "worked" during earlier ad hoc testing because the
  package happened to be installed on this dev machine.
- [x] **The step-halving / rank-repair branch (`R/ecm_msfr.R:249-267`)** — all three bugs
  (unused shrunk step, wrong field name from `.vect2param()`, wrong variable name for the
  inverse-Ψ update) are fixed on inspection. **Still not exercised by any test** — nothing in
  the current suite forces `Omega` to be rank-deficient, so this fix is verified by reading,
  not by execution. Still worth a dedicated test that deliberately triggers it (see P2).
- [x] **New bug found while verifying the `.update_Lambda()` fix, same class as `ecm_fr()`'s
  earlier `Phi` bug: `res$psi_s` (the returned uniquenesses) was set once from
  `start$psi_s` and never reassigned in the loop.** The internal `Psi_s`/`Psi_s1` (matrix
  form, used for the E-step) *did* update correctly every iteration via `Psi_s <- Psi_new` —
  but the lowercase `psi_s` (vector form, the thing actually returned to the caller) has
  exactly one assignment in the whole function, at initialization. Confirmed empirically:
  after the `.update_Lambda()` fix, `test-ecm_msfr.R`'s parameter-movement test failed with
  `res$psi_s` differing from `start$psi_s` by exactly `0` after 25 iterations, while
  `Phi`/`Lambda_s`/`beta` all moved substantially — and `grep -n "\bpsi_s\b"` confirmed there
  is genuinely no reassignment anywhere in the loop body. This also explains why the recovery
  test's `Sigma_hat` reconstruction was off (it's built from `res$psi_s`). Fixed by adding
  `psi_s <- psi_new` alongside the other "assign vars for new cycle" carry-overs — `psi_new`
  (vector form) was already being correctly computed each iteration for the identifiability
  constraint, just never carried into the variable actually returned.

This is worth flagging as a pattern, not just two isolated bugs: this session has now found
the *same* bug shape three times independently — `ecm_fr()`'s `Phi`, and now `ecm_msfr()`'s
`Lambda`-reshape and `psi_s` — a parameter gets a fresh `*_new` value computed correctly every
iteration, but the "carry it into next iteration / into the return value" step is missing or
targets the wrong variable. Worth specifically checking `ecm_fa()` and any future modules for
this exact shape rather than assuming it's been ruled out just because those functions'
existing tests pass.

### Naming: the `ecm_msfa` → `ecm_msfr` rename is half-done, and the package is currently broken because of it

The function itself has been renamed: `R/ecm_msfr.R:107` now defines `ecm_msfr <- function(...)`,
not `ecm_msfa`. But nothing downstream of that was regenerated:

- [x] **`NAMESPACE` still has `export(ecm_msfa)`, and `ecm_msfa` no longer exists anywhere in
  `R/`.** This isn't a documentation nicety — it means `library(MSFR)` would currently export
  a name that resolves to nothing. Confirmed directly: running the new test suite via
  `devtools::test()` prints `Warning message: Objects listed as exports, but not present in
  namespace: * ecm_msfa`. Fix by re-running `devtools::document()` once `.update_Lambda()` is
  fixed (roxygen will pick up the new name from the `@export`ed function and regenerate
  `NAMESPACE` — no manual edit needed).
- [x] The roxygen title above the function (`R/ecm_msfr.R:63`) still reads *"Estimates the
  parameters of a MSFA model"* — needs updating to describe the actual (covariate-adjusted,
  MSFR) model; `devtools::document()` won't fix prose, only structure.
- [x] `file-structure.md:2` still describes `R/ecm_msfr.R` only in terms of "the multi-study
  factor regression algorithm" without naming the function, so it doesn't contradict the new
  name, but it's worth a pass once the rename is fully settled.
- [x] `README.md:31` calls `ecm_msfa(...)` — needs updating to `ecm_msfr(...)` as part of the
  broader README rewrite already tracked in P3.

### `R/helpers.R` — verify the fix already made this session holds

Earlier this session, a dead, pre-refactor duplicate of `.exp_values_fr()` living in
`helpers.R` was found to silently shadow the correct one in `exp_values_fr.R` under R's
default (alphabetical, no `Collate:` field) package-load order — meaning `ecm_fr()` would
work in ad hoc interactive testing but silently break in a real package build. It's been
removed.

- [ ] `tests/testthat/` now exists (see P2) — add a regression test that just asserts
  `formals(.exp_values_fr)` has the expected names (`Phi, Psi_s, cov_s, X_s_tilde, getdet` —
  `Psi_s1` was dropped from the signature this session, since `.exp_values_fr()` now computes
  the inverse internally; see `test-ecm_msfr.R` for the pattern to follow, or add it as its
  own `test-exp_values_fr.R`). Cheap, and it's exactly the kind of silent, load-order-dependent
  breakage that's easy to reintroduce without noticing.

---

## P1 — Package infrastructure (the package doesn't currently build)

Verified directly: `devtools::check()` on the real `DESCRIPTION` fails immediately at the
build stage, before any code is even checked.

- [x] **`Authors@R` has no one with role `"cre"` (maintainer) and a valid email.**
  `R CMD build` refuses outright: *"Authors@R field gives no person with maintainer role,
  valid email address and non-empty name."* Currently: Roberta De Vito has `"aut"` + email,
  Alejandra Avalos-Pacheco has `"aut"` with no email, Jörg Schantz has `"ctb"` with no email.
  Someone needs a `"cre"` role and a real, working email address — presumably Jörg Schantz,
  given the project's current maintenance activity, but that's a call for the actual authors
  to make, not something to guess at here.
- [x] **`DESCRIPTION` has no `Imports:` field at all**, despite the code depending on four
  external packages. Confirmed by `R CMD check` on a patched copy: *"Namespace dependencies
  missing from DESCRIPTION Imports/Depends entries: 'pracma', 'psych', 'robust', 'statmod'."*
  This isn't just a CRAN nicety — without it, `install.packages()`/`remotes::install_github()`
  won't pull in the packages this code actually needs. Add via
  `usethis::use_package("robust")` / `usethis::use_package("psych")` /
  `usethis::use_package("statmod")` / `usethis::use_package("pracma")` (the last one is used
  by `vcov_msfa.R`, out of scope for code review here, but the dependency still needs
  declaring at the package level).
- [x] **`NAMESPACE` has three whole-package `import(...)` lines** (`psych`, `robust`,
  `statmod`) where only a handful of specific functions are actually used
  (`covRob`/`fa`/`vecmat`/`matvec`/etc.). Whole-package `@import` pulls the entire exported
  surface of each dependency into this package's namespace, which is both unnecessary and a
  latent symbol-collision risk (exactly the kind of thing that made the accidental
  `WoodburyMatrix::t()` interaction above so easy to miss). Convert to `@importFrom` for the
  specific functions actually called, matching how `stats::cor/cov/factanal/prcomp` are
  already (correctly) imported.
- [x] **The `Description:` field in `DESCRIPTION` is the paper's abstract, pasted with its
  original line-wrap hyphenation intact** — "impor- tant", "dif- ferent", "covari- ate",
  "His- panic" all appear as broken words. This is the text CRAN (and every user running
  `install.packages()`) would see. Needs a rewrite as a proper package description (what the
  package *does*, in software terms — fits MSFR/FR/FA models via ECM — not the paper's
  scientific abstract).
- [ ] Add a package-level help page. `man/MSFR-package.Rd` exists but there's no
  corresponding `"_PACKAGE"` roxygen source anywhere in `R/` — confirmed:
  `devtools::document()` deletes `MSFR-package.Rd` as an orphan rather than regenerating it.
  Re-run `usethis::use_package_doc()` and re-document.
- [ ] `Data/Scenario1_MSFR.rda` is a bundled dataset with **no documentation page**
  (`man/Scenario1_MSFR.Rd` doesn't exist). CRAN requires every bundled dataset to be
  documented (`@format`, `@source`, description of each column/list element). Also worth
  double-checking it's actually current — it's referenced by `README.md`'s example (see
  below), which itself looks unmaintained.
- [x] `R CMD build` also emits a NOTE that the package now depends on R ≥ 3.5.0 because of
  how `Scenario1_MSFR.rda` was serialized. Either re-save it with an older serialization
  version or just declare `Depends: R (>= 3.5.0)` explicitly in `DESCRIPTION` so this is a
  documented decision rather than an implicit one.
- [ ] Bump `Version:` off the `usethis::create_package()` default (`0.0.0.9000`) once the
  above are addressed and there's a first coherent state worth tagging.

---

## P2 — Testing (infrastructure now started)

- [x] `usethis::use_testthat(3)` has been run: `tests/testthat/` and `tests/testthat.R` now
  exist, and `DESCRIPTION` declares `Suggests: testthat (>= 3.0.0)` +
  `Config/testthat/edition: 3`.
- [x] `tests/testthat/test-ecm_msfr.R` exists and covers the MSFR model (E-step brute-force
  cross-check, recovery test, parameter-movement test — see P0 for current pass/fail status).
  Run with `devtools::test(filter = "ecm_msfr")` or `devtools::test()` for the whole suite.
- [ ] Everything else is still using disposable scratch scripts, run and then deleted —
  nothing beyond the one file above is left behind to catch a regression yet. Given how many
  silent, load-order- or scoping-dependent bugs have turned up in this codebase (the
  `helpers.R` shadowing issue, `ecm_fr()`'s `Phi` never being carried forward, `ecm_msfr.R`'s
  `.update_Lambda()` bugs), filling this in for the rest of the package is not optional
  polish — it's the difference between catching the next one of these in five seconds and
  re-discovering it by accident months later. Port the verification patterns already proven
  out this session into permanent tests, following `test-ecm_msfr.R`'s pattern:
  - Brute-force numerical cross-checks: compare `.exp_values_fa()`/`.exp_values_fr()`
    against a direct `solve()`-based reference implementation (this was done ad hoc multiple
    times this session — turn it into `test-exp_values_fa.R` / `test-exp_values_fr.R`).
  - End-to-end recovery tests: simulate data from a known `Λ`/`Φ`/`β`/`Ψ`, fit with
    `ecm_fa()`/`ecm_fr()`, assert the fit converges and reconstructs the simulated
    covariance/coefficients within a reasonable tolerance (`ecm_msfr()`'s version of this is
    now in place, but currently failing — see P0).
  - Regression tests for the specific bug classes already found: `Φ`/`Λₛ`/`β` must all
    differ from their initial values after fitting (this pattern is in
    `test-ecm_msfr.R`; port the same check to `ecm_fa()`/`ecm_fr()`); `.exp_values_fr` must
    have exactly the current argument list (catches the `helpers.R` shadowing bug class —
    still not written, see the `helpers.R` item above).
  - `q_s = 1` / `k = 1` edge cases for all three model-fitting functions (single-factor
    models) — spot-checked ad hoc for FA and FR earlier this session (both passed) and never
    turned into a permanent test; `test-ecm_msfr.R`'s default scenario uses `j_s = c(1, 1)`
    (single study-specific factor per study) but `k = 2`, so MSFR's `k = 1` case specifically
    is still untested.
  - A rank-deficiency case that actually forces `ecm_msfr()`'s step-halving branch, once
    that branch has a test forcing it to trigger at all (see P0).
- [ ] `Local_testing/` (currently git-ignored, not part of the package) already has real,
  substantive scratch scripts — `simple_scenario.R`, `WB_testing.R`, `runtimes.txt` — that
  look like a natural starting point to mine for real `testthat` cases rather than writing
  everything from scratch. Worth a read before starting P2 in earnest. It also contains its
  own copy of `start_msfa.R`, separate from `R/start_msfa.R` — check whether it's a stale
  fork or something that was deliberately being tested; if stale, it's just noise.
- [ ] Once tests exist: `devtools::test_coverage()` to see what's actually exercised, then
  `usethis::use_github_action('test-coverage')` + coverage badge per the CRAN checklist.
- [ ] Add `devtools::check()` as a routine step from here on — every substantial change in
  this session was verified with hand-rolled scripts using `devtools::load_all()`, which is
  the right instinct, but `devtools::check()` (once P1 unblocks it) catches a much wider
  class of issues automatically (as this audit's use of it already demonstrated).

---

## P3 — Documentation

- [ ] **Re-run `devtools::document()`.** `man/` is stale relative to `R/`: this session added
  `.exp_values_fa()` and rewrote `.exp_values_fr()`, both fully roxygen-documented in source,
  but neither has a rendered `.Rd` page yet, and `ecm_fa.Rd`/`ecm_fr.Rd` in `man/` predate the
  bug fixes made to those functions' implementations. (`RoxygenNote`/`Config/roxygen2/version`
  in `DESCRIPTION` should also be reconciled with whatever `roxygen2` version is actually
  used to render — a check-copy render bumped `RoxygenNote` to 7.3.3 while `DESCRIPTION`
  currently only pins `Config/roxygen2/version: 8.0.0`.)
- [ ] **`README.md` needs a full rewrite, not a touch-up.** Specific problems found by
  reading it closely:
  - It's hand-written with RMarkdown-style code fences (` ```{r ...} `) but saved as a plain
    `.md`, so none of the "examples" are ever actually executed/validated — exactly the
    failure mode `usethis::use_readme_rmd()` (unchecked in `CRAN_checklist.md`) exists to
    prevent.
  - From the "Fitting the model via ECM" section onward, code fences are opened but never
    closed — every subsequent section is nested inside one giant unterminated code block when
    rendered.
  - `data(Scenario1_MSFR.rda)` is invalid usage — `data()` takes the object name
    (`Scenario1_MSFR`), not the filename with extension.
  - The final line, `` beta <- ECM_MSFR$$beta ``, has a stray extra `$` — not valid R.
  - It calls `ecm_msfa(...)`, correctly matching the actual current export name — worth
    noting given the naming confusion flagged in P0; whichever name is settled on there needs
    to be reflected here too.
  - It only documents the MSFR flow — nothing about `ecm_fa()`/`ecm_fr()`, which are
    both exported, user-facing entry points.
  - Fix: `usethis::use_readme_rmd()`, write one working example per exported model-fitting
    function, `devtools::build_readme()` to actually knit and validate it.
- [ ] Add a vignette (`usethis::use_vignette(...)`) once the README example is solid — a
  natural one would walk through simulating data, fitting all three models
  (`ecm_fa`/`ecm_fr`/`ecm_msfa`), and comparing them, mirroring the paper's own simulation
  study structure.
- [ ] `pkgdown` site (`usethis::use_pkgdown()`, `pkgdown::build_site()`,
  `usethis::use_pkgdown_github_pages()`) — low effort once docs/README/vignette are in good
  shape, and it's explicitly on the `CRAN_checklist.md` list.
- [ ] `usethis::use_news_md()` — start tracking user-visible changes from here forward,
  given how much has already changed this session (E-step consolidation across FA/FR,
  multiple correctness fixes).
- [ ] `cffr::cff_write()` for a citation file, given this package exists specifically to
  accompany published research (the MSFR/MSFA/FR papers).

---

## P4 — Repository hygiene

None of these block functionality, but they affect what actually ends up in a built package
tarball or a fresh clone, and a couple are close to becoming real problems.

- [ ] **Two CSV files at the package root, `likelihoods_253.csv` (2.8MB) and
  `likelihoods_no253.csv` (1.5MB), are not excluded by `.Rbuildignore`.** Together they're
  ~4.3MB, a large fraction of CRAN's soft package-size expectations, and they look like
  debug output — `ecm_msfr.R:311` has a commented-out
  `write.csv2(l.df, "likelihoods_no253.csv")`. If they're not meant to ship with the
  package, add them to `.Rbuildignore`; if they're reference data for a vignette or test,
  they belong under `inst/extdata/` with a comment explaining what they are, not sitting
  unexplained at the root.
- [ ] `.Rbuildignore` currently only excludes `MSFR.Rproj`, `.Rproj.user`, and `LICENSE.md`.
  It should also cover: `Local_testing/` (already git-ignored, but `.Rbuildignore` is
  independent of git and should agree), `.DS_Store`, and the project-internal planning docs
  at the root (`CRAN_checklist.md`, `file-structure.md`, `initial.md`,
  `project-requirements.md`, and this file) — `R CMD build`/`check` will otherwise flag these
  as non-standard top-level files/directories.
- [ ] **`R/v_Matrix/` (`exp_msfa_sparse.R`, `start_msfa_sparse.R`) is unfinished exploratory
  code for a sparse/`Matrix`-and-`irlba`-backed variant of MSFA.** Confirmed it is *not*
  currently picked up by the package build (R only sources top-level files directly under
  `R/`, not subdirectories, and this was verified directly — the live `.exp_values` in a
  loaded package has none of this file's content) — so it's inert today, not a live risk.
  But it redefines `.exp_values`, `.wb_identity`, `.wb_identity2`, `.get_exp_xl`, etc. under
  the *exact same names* as the real, working versions, references an undefined global
  `test` object (clearly leftover from interactive development), and depends on `Matrix`
  and `irlba`, neither declared anywhere. If anyone ever "helpfully" flattens `R/` or adds a
  `Collate:` field that pulls subdirectories in, this becomes the same silent-shadowing
  hazard already found and fixed in `helpers.R`. Either finish and properly integrate it
  (own file names, no naming collisions, declared dependencies), or move it out of `R/`
  entirely (e.g. into `Local_testing/` or a `dev/` folder) until it's ready.
- [ ] `file-structure.md` is already slightly out of sync with the real tree (it references
  `R/exp_values_msfr.R` and `R/start_msfr.R`, which don't exist — the real files are
  `R/exp_values.R` and `R/start_msfa.R`) and needs a pass to also document
  `R/exp_values_fa.R`/`R/exp_values_fr.R` accurately (their entries were added this session
  but the pre-existing MSFR-related entries were never corrected). Worth fixing in the same
  pass as the naming cleanup in P0.

---

## P5 — Known design debt (not blocking, worth tracking)

Lower-priority items already identified and deliberately deferred earlier this session,
recorded here so they aren't lost:

- [ ] `.exp_values()` (the general MSFR E-step, `exp_values.R`) still passes a dense `p × p`
  matrix for `Ψₛ⁻¹` in its first-level Woodbury reduction (`wb_f`/`wb_l`). `.wb_identity()`/
  `.wb_identity2()` now support a vector fast-path for exactly this case (added this session
  for FA/FR, ~1.6–2× measured speedup), but MSFR's first-level reduction was deliberately
  left on the dense path since it wasn't in scope. The second-level reduction genuinely can't
  use the vector path (its `W` argument isn't diagonal), but the first-level one could.
- [ ] `ecm_fr()`'s `npar` (and therefore `AIC`/`BIC`) reuses `ecm_fa()`'s formula
  (`p*S + sum(tot_s*(p-(tot_s-1)/2))`), which counts loading parameters once per study — but
  FR has a single *shared* `Φ` across studies, not a per-study copy, and the formula omits
  `β`'s `p*p_b` parameters entirely. `AIC`/`BIC` from `ecm_fr()` are currently not
  trustworthy for model comparison.
- [ ] `ecm_fr()` accepts `block_lower`, `robust`, `corr`, `mcd` parameters but doesn't use
  any of them: `Φ`'s upper triangle is never zeroed (no identifiability constraint enforced,
  unlike `ecm_fa()` and `ecm_msfa()`, both of which do enforce it), and `cov_s` is always
  the plain sample covariance regardless of `robust`/`corr`/`mcd`.
- [ ] `ecm_fa()` accepts a `traceIT` parameter (documented as controlling trace-print
  frequency) but the loop hardcodes `i %% 100` instead of `i %% traceIT`.
- [ ] Minor style/consistency cleanup once the above is settled: `@import psych` in
  `ecm_fa.R`/`ecm_fr.R`'s roxygen block is unused in both files (should be dropped once P1's
  `importFrom` cleanup happens); a couple of `solve(A) %*% B` patterns could be
  `solve(A, B)` for a small numerical-stability/performance win, same reasoning already
  applied elsewhere this session.

---

## Suggested order of attack

1. **P0** — nothing else is trustworthy until `ecm_msfa()` actually runs and every
   model-fitting function has been confirmed working via `devtools::load_all()`, not ad hoc
   `source()`.
2. **P1** — get `devtools::check()` to run past the build stage at all; this is also what
   unblocks meaningfully testing P0's fixes in a real package context.
3. **P2** — lock in everything fixed so far with a real test suite before doing anything else,
   given how many of the bugs found this session were silent (wrong answers or load-order-
   dependent breakage, not crashes).
4. **P3** — once the package is correct and tested, make it legible to someone who isn't the
   author.
5. **P4** — cheap, do in idle moments; won't conflict with anything above.
6. **P5** — genuine design decisions, worth a deliberate conversation with the paper authors
   (especially the FR identifiability/robust-covariance gaps) rather than a unilateral fix.
