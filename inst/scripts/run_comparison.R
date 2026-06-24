# ============================================================
# run_comparison.R  --  produces manuscript Table 3
#
# Three arms scored on identical scenarios, data generation, and
# decision rule, varying ONLY the dose-response model / criterion:
#   (1) BART utility   : your validated engine (reproduces Table 2)
#   (2) Parametric util: same decision layer, BART swapped for a
#                        Bayesian probit dose-response (utility design
#                        without BART)  -- the apples-to-apples foil
#   (3) EffTox         : canonical trade-off-contour design via trialr
#                        /escalation  -- the external anchor
#
# Place this in bartEffTox/inst/scripts/ (next to run_simulations.R)
# and run after the package is built/loaded.
#
# --------------------------- READ FIRST ---------------------------
# * Authored without an R session here, so treat as a script to run
#   and check. The first thing to verify is that arm (1) reproduces
#   Table 2 within Monte Carlo error; if it does not, the wiring is
#   off, not the science.
# * Arm (2) reuses your bart_efftox_interim decision logic VERBATIM
#   (PAVA, admissibility eq 2, utility eq 3, no-skip, tie-break);
#   only the model-fit block is replaced. The parametric model is an
#   independent Bayesian probit in [1, dose, dose^2] per outcome
#   (Albert-Chib), regularized by a N(0, b0 I) prior, with covariates
#   entered linearly and predictions marginalized over the reference
#   covariates exactly as your interim does. Document this choice in
#   the paper; if you prefer the original 6-parameter EffTox model as
#   the parametric foil, arm (3) already supplies the canonical one.
# * Arm (3) uses its own validated conduct (escalation/trialr). Its
#   acceptability is matched to your admissibility (p_t = 1 - c_T,
#   p_e = 1 - c_E); its start-low/escalate-by-recommendation conduct
#   is intrinsic to EffTox and is not forced to match your run-in.
#   Stan is slow: keep its R modest and report Monte Carlo SEs.
# * All arms are scored against the SAME target, scenario$true_opt
#   (the utility-optimal admissible dose), with P_none reported for
#   the no-acceptable-dose scenarios.
# ============================================================

# ---- load your engine (run_trial, simulate_design, summarise_ocs, etc.) ----
# Tries, in order: an already-loaded package, devtools::load_all on the repo
# root found by walking up from this script or the working directory, then an
# installed library(). Stops with a clear message if none expose the engine.
.need <- c("run_trial","simulate_design","summarise_ocs","bart_efftox_interim",
           "pava_project","utility_fn","runin_admissibility")
.have <- function() all(vapply(.need, exists, logical(1), mode = "function"))

if (!.have()) {
  # locate the package root (the dir containing DESCRIPTION for bartEffTox)
  cand <- c(".", "..", "../..", "../../..",
            tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA))
  root <- NA_character_
  for (p in cand) {
    if (is.na(p)) next
    d <- file.path(p, "DESCRIPTION")
    if (file.exists(d) && any(grepl("bartEffTox", readLines(d), fixed = TRUE))) {
      root <- normalizePath(p); break
    }
  }
  if (!is.na(root) && requireNamespace("devtools", quietly = TRUE)) {
    message("Loading engine via devtools::load_all('", root, "')")
    devtools::load_all(root, quiet = TRUE)
  } else if (requireNamespace("bartEffTox", quietly = TRUE)) {
    message("Loading engine via library(bartEffTox)")
    library(bartEffTox)
  }
}
if (!.have())
  stop("bartEffTox engine not found. Run from the repo (so DESCRIPTION is ",
       "reachable) with devtools installed, or install the package first:\n",
       "  devtools::load_all('/path/to/BART-EffTox/bartEffTox')\n",
       "then source this script. Missing: ",
       paste(.need[!vapply(.need, exists, logical(1), mode='function')], collapse=", "))

# Scenarios embedded verbatim from inst/scripts/run_simulations.R, so the
# script does not depend on how the packaged data object is named or loaded.
scenarios <- list(
  S1 = list(name="Parametric benchmark",
            pT=c(0.05,0.10,0.18,0.25,0.45), pE=c(0.10,0.25,0.40,0.55,0.65),
            covars=FALSE, N=36L, true_opt=4L,
            true_adm=c(FALSE,TRUE,TRUE,TRUE,FALSE)),
  S2 = list(name="Plateau efficacy",
            pT=c(0.05,0.11,0.30,0.45,0.62), pE=c(0.15,0.42,0.48,0.50,0.51),
            covars=FALSE, N=36L, true_opt=2L,
            true_adm=c(FALSE,TRUE,TRUE,FALSE,FALSE)),
  S3 = list(name="Non-monotone efficacy",
            pT=c(0.04,0.09,0.17,0.28,0.48), pE=c(0.10,0.28,0.50,0.38,0.22),
            covars=FALSE, N=36L, true_opt=3L,
            true_adm=c(FALSE,TRUE,TRUE,TRUE,FALSE)),
  S4 = list(name="Biomarker-modified efficacy", pB=0.40,
            pT_xb0=c(0.05,0.10,0.18,0.32,0.50), pT_xb1=c(0.05,0.10,0.18,0.32,0.50),
            pE_xb0=c(0.08,0.15,0.22,0.25,0.25), pE_xb1=c(0.15,0.35,0.60,0.70,0.72),
            covars=TRUE, N=36L, true_opt=3L,
            true_adm=c(FALSE,TRUE,TRUE,FALSE,FALSE),
            pT=c(0.05,0.10,0.18,0.32,0.50),
            pE=0.60*c(0.08,0.15,0.22,0.25,0.25)+0.40*c(0.15,0.35,0.60,0.70,0.72)),
  S5 = list(name="Frailty-modified toxicity", pF=0.30,
            pT_xf0=c(0.02,0.05,0.10,0.19,0.33), pT_xf1=c(0.08,0.20,0.38,0.58,0.78),
            pE_xf0=c(0.10,0.28,0.46,0.55,0.56), pE_xf1=c(0.10,0.28,0.46,0.55,0.56),
            covars=TRUE, N=36L, true_opt=3L,
            true_adm=c(FALSE,TRUE,TRUE,FALSE,FALSE),
            pT=0.70*c(0.02,0.05,0.10,0.19,0.33)+0.30*c(0.08,0.20,0.38,0.58,0.78),
            pE=c(0.10,0.28,0.46,0.55,0.56)),
  S6 = list(name="Small-sample stress test",
            pT=c(0.06,0.14,0.28,0.42,0.58), pE=c(0.14,0.38,0.42,0.40,0.32),
            covars=FALSE, N=24L, true_opt=2L,
            true_adm=c(FALSE,TRUE,TRUE,FALSE,FALSE)),
  S7 = list(name="Combined modern oncology", pB=0.40, pF=0.30,
            pT_xf0=c(0.03,0.07,0.14,0.24,0.42), pT_xf1=c(0.07,0.16,0.30,0.50,0.70),
            pE_xb0=c(0.08,0.18,0.28,0.25,0.18), pE_xb1=c(0.15,0.40,0.65,0.55,0.38),
            covars=TRUE, N=36L, true_opt=3L,
            true_adm=c(FALSE,TRUE,TRUE,FALSE,FALSE),
            pT=0.70*c(0.03,0.07,0.14,0.24,0.42)+0.30*c(0.07,0.16,0.30,0.50,0.70),
            pE=0.60*c(0.08,0.18,0.28,0.25,0.18)+0.40*c(0.15,0.40,0.65,0.55,0.38))
)

R_COMP      <- 2000L
R_EFFTOX    <- 1000L        # Stan is slow; smaller R, report MCSE
N_CORES     <- 4L
MASTER_SEED <- 2025L
DESIGN <- list(cohort = 3L, phi_T = 0.30, phi_E = 0.20, c_T_pre = 0.70,
               c_T = 0.40, c_E = 0.50, a_T = 1L, b_T = 1L,
               M = 20L, n_iter = 400L, n_burn = 80L)

## ============================================================
## Arm (2) building blocks: parametric probit interim
## ============================================================

# Vectorized truncated-normal sampler via inverse CDF.
.rtn <- function(mu, y) {           # y=1 -> (0, Inf); y=0 -> (-Inf, 0)
  lo <- ifelse(y == 1L, 0, -Inf); hi <- ifelse(y == 1L, Inf, 0)
  Plo <- pnorm(lo - mu); Phi <- pnorm(hi - mu)
  u <- runif(length(mu), Plo, Phi)
  mu + qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
}

# Albert-Chib probit Gibbs. Returns n_iter x nrow(Xpred) prob draws.
.ac_probit <- function(Xtr, y, Xpred, n_iter, n_burn, b0 = 4) {
  k  <- ncol(Xtr)
  V  <- solve(crossprod(Xtr) + diag(1 / b0, k))
  cV <- chol(V)
  beta <- rep(0, k)
  out  <- matrix(NA_real_, n_iter, nrow(Xpred))
  for (it in seq_len(n_burn + n_iter)) {
    z    <- .rtn(as.vector(Xtr %*% beta), y)
    bhat <- V %*% crossprod(Xtr, z)
    beta <- as.vector(bhat + t(cV) %*% rnorm(k))
    if (it > n_burn) out[it - n_burn, ] <- pnorm(as.vector(Xpred %*% beta))
  }
  out
}

# Drop-in replacement for bart_efftox_interim: SAME signature and return,
# parametric probit in place of BART. Decision logic copied verbatim.
parametric_interim <- function(dose_obs, yT, yE, X_cov = NULL,
                               dose_grid = 1:5, X_ref = NULL,
                               current_dose = 1L,
                               phi_T = 0.30, phi_E = 0.20,
                               c_T = 0.40, c_E = 0.50,
                               M = 20L, n_iter = 400L, n_burn = 80L) {
  J             <- length(dose_grid)
  dose_std      <- 2 * (dose_obs  - 1L) / (J - 1L) - 1
  dose_std_grid <- 2 * (dose_grid - 1L) / (J - 1L) - 1
  basis <- function(d, Xc = NULL) {
    B <- cbind(1, d, d^2)
    if (!is.null(Xc)) B <- cbind(B, Xc)
    B
  }
  Xtr <- basis(dose_std, X_cov)
  if (!is.null(X_cov) && !is.null(X_ref)) {
    n_ref  <- nrow(X_ref)
    Xpred  <- do.call(rbind, lapply(dose_std_grid,
                function(d) basis(rep(d, n_ref), X_ref)))
    marginal <- TRUE
  } else {
    Xpred <- basis(dose_std_grid); n_ref <- 1L; marginal <- FALSE
  }

  pT_raw <- .ac_probit(Xtr, yT, Xpred, n_iter, n_burn)
  pE_raw <- .ac_probit(Xtr, yE, Xpred, n_iter, n_burn)

  if (marginal) {
    pT_samp <- pE_samp <- matrix(NA_real_, n_iter, J)
    for (j in seq_len(J)) {
      idx          <- ((j - 1L) * n_ref + 1L):(j * n_ref)
      pT_samp[, j] <- rowMeans(pT_raw[, idx, drop = FALSE])
      pE_samp[, j] <- rowMeans(pE_raw[, idx, drop = FALSE])
    }
  } else { pT_samp <- pT_raw; pE_samp <- pE_raw }

  # ---- from here, IDENTICAL to bart_efftox_interim ----
  n_j <- tabulate(dose_obs, nbins = J); w_j <- n_j + 1L
  pT_proj <- matrix(NA_real_, n_iter, J)
  for (s in seq_len(n_iter)) pT_proj[s, ] <- pava_project(pT_samp[s, ], w_j)

  pr_toxic   <- colMeans(pT_proj > phi_T)
  pr_ineffic <- colMeans(pE_samp < phi_E)
  adm <- (pr_toxic < c_T) & (pr_ineffic < c_E)

  safety_flag <- mean(pT_proj[, current_dose] > phi_T) > c_T
  eligible    <- which(adm & (seq_len(J) <= current_dose + 1L))
  if (safety_flag) eligible <- eligible[eligible <= current_dose]

  if (length(eligible) == 0L)
    return(list(next_dose = NA_integer_, stop_trial = TRUE,
                admissible_set = integer(0L),
                mean_utility = rep(NA_real_, J),
                post_pT = pT_proj, post_pE = pE_samp,
                pr_toxic = round(pr_toxic, 3L),
                pr_ineffic = round(pr_ineffic, 3L)))

  mean_util <- numeric(J)
  for (j in seq_len(J)) mean_util[j] <- mean(utility_fn(pE_samp[, j], pT_proj[, j]))
  best_util <- max(mean_util[eligible])
  tied      <- eligible[abs(mean_util[eligible] - best_util) < 1e-6]
  best      <- if (current_dose %in% tied) current_dose else min(tied)

  list(next_dose = best, stop_trial = FALSE, admissible_set = eligible,
       mean_utility = round(mean_util, 3L),
       post_pT = pT_proj, post_pE = pE_samp,
       pr_toxic = round(pr_toxic, 3L), pr_ineffic = round(pr_ineffic, 3L))
}

# run_trial with the interim generalized (body copied from run_trial).
run_trial_arm <- function(scenario, interim_fn, max_n = 36L, cohort = 3L,
                          phi_T = 0.30, phi_E = 0.20, c_T_pre = 0.70,
                          c_T = 0.40, c_E = 0.50, a_T = 1L, b_T = 1L,
                          M = 20L, n_iter = 400L, n_burn = 80L,
                          use_covars = FALSE) {
  J <- length(scenario$pT)
  dose_obs <- integer(0); yT <- integer(0); yE <- integer(0)
  X_cov_obs <- NULL; current <- 1L; bart_active <- FALSE
  gen_cohort <- function(d, m) {
    if (isTRUE(scenario$covars) && use_covars) {
      xb <- if (!is.null(scenario$pB)) rbinom(m, 1L, scenario$pB) else NULL
      xf <- if (!is.null(scenario$pF)) rbinom(m, 1L, scenario$pF) else NULL
      pT_i <- if (!is.null(xf)) ifelse(xf == 1L, scenario$pT_xf1[d], scenario$pT_xf0[d]) else scenario$pT[d]
      pE_i <- if (!is.null(xb)) ifelse(xb == 1L, scenario$pE_xb1[d], scenario$pE_xb0[d]) else scenario$pE[d]
      cols <- Filter(Negate(is.null), list(xb = xb, xf = xf))
      X_new <- if (length(cols) > 0) do.call(cbind, cols) else NULL
    } else { pT_i <- scenario$pT[d]; pE_i <- scenario$pE[d]; X_new <- NULL }
    list(yT = rbinom(m, 1L, pT_i), yE = rbinom(m, 1L, pE_i),
         d = rep(d, m), X = X_new)
  }
  dose_safe <- function(d) {
    nj <- sum(dose_obs == d); tj <- sum(yT[dose_obs == d])
    if (nj == 0L) return(TRUE)
    !runin_admissibility(tj, nj, phi_T, c_T_pre, a_T, b_T)
  }
  for (d in seq_len(J)) {
    if (length(dose_obs) >= max_n) break
    nxt <- gen_cohort(d, cohort)
    dose_obs <- c(dose_obs, nxt$d); yT <- c(yT, nxt$yT); yE <- c(yE, nxt$yE)
    if (!is.null(nxt$X)) X_cov_obs <- rbind(X_cov_obs, nxt$X)
    current <- d
    if (!dose_safe(d)) { current <- max(1L, d - 1L); bart_active <- TRUE; break }
    if (d == J) bart_active <- TRUE
  }
  if (length(dose_obs) >= max_n) bart_active <- TRUE
  rec_dose <- NA_integer_; stopped_early <- FALSE
  if (bart_active) {
    while (length(dose_obs) < max_n) {
      X_ref <- if (use_covars && !is.null(X_cov_obs)) X_cov_obs else NULL
      res <- interim_fn(dose_obs = dose_obs, yT = yT, yE = yE, X_cov = X_cov_obs,
                        dose_grid = seq_len(J), X_ref = X_ref, current_dose = current,
                        phi_T = phi_T, phi_E = phi_E, c_T = c_T, c_E = c_E,
                        M = M, n_iter = n_iter, n_burn = n_burn)
      if (res$stop_trial) { stopped_early <- TRUE; break }
      current <- res$next_dose
      nxt <- gen_cohort(current, cohort)
      dose_obs <- c(dose_obs, nxt$d); yT <- c(yT, nxt$yT); yE <- c(yE, nxt$yE)
      if (!is.null(nxt$X)) X_cov_obs <- rbind(X_cov_obs, nxt$X)
    }
  }
  if (!stopped_early && length(dose_obs) >= max_n) {
    X_ref <- if (use_covars && !is.null(X_cov_obs)) X_cov_obs else NULL
    final <- interim_fn(dose_obs = dose_obs, yT = yT, yE = yE, X_cov = X_cov_obs,
                        dose_grid = seq_len(J), X_ref = X_ref, current_dose = current,
                        phi_T = phi_T, phi_E = phi_E, c_T = c_T, c_E = c_E,
                        M = M, n_iter = n_iter, n_burn = n_burn)
    if (!final$stop_trial) rec_dose <- final$next_dose
  }
  list(rec_dose = rec_dose, stopped_early = stopped_early,
       n = length(dose_obs), n_at_dose = tabulate(dose_obs, nbins = J),
       mean_tox = mean(yT), mean_eff = mean(yE))
}

# simulate_design with a pluggable trial fn (body copied from simulate_design).
simulate_arm <- function(scenario, interim_fn, R = 2000L, n_cores = 2L,
                         seed = 42L, chunk_size = 50L, ...) {
  set.seed(seed); seeds <- sample.int(1e6L, R)
  chunks <- split(seq_len(R), ceiling(seq_len(R) / chunk_size))
  out <- vector("list", R)
  for (ch in chunks) {
    ch_res <- parallel::mclapply(ch, function(r) {
      set.seed(seeds[r])
      tryCatch(run_trial_arm(scenario, interim_fn, ...),
               error = function(e) structure(list(message = conditionMessage(e)),
                                             class = "trial_error"))
    }, mc.cores = n_cores)
    for (i in seq_along(ch)) out[[ch[i]]] <- ch_res[[i]]
  }
  out
}

## ============================================================
## Arm (3): canonical EffTox via escalation/trialr (optional, slow)
## ============================================================
efftox_ocs <- function(sc, R = R_EFFTOX, cohort = 3L, seed = MASTER_SEED) {
  if (!requireNamespace("escalation", quietly = TRUE) ||
      !requireNamespace("trialr", quietly = TRUE))
    return(list(PCS = NA, P_none = NA))
  J  <- length(sc$pT); rd <- seq_len(J)
  ## TODO: elicit the contour points and ESS-calibrate the priors to your
  ## trial. Values below are the trialr tutorial example -- placeholders.
  design <- escalation::get_trialr_efftox(
    real_doses = rd, efficacy_hurdle = 0.20, toxicity_hurdle = 0.30,
    p_e = 0.50, p_t = 0.60,                       # matched to c_E, c_T
    eff0 = 0.5, tox1 = 0.65, eff_star = 0.7, tox_star = 0.25,
    alpha_mean = -7.9593, alpha_sd = 3.5487, beta_mean = 1.5482, beta_sd = 3.5018,
    gamma_mean = 0.7367, gamma_sd = 2.5423, zeta_mean = 3.4181, zeta_sd = 2.4406,
    eta_mean = 0, eta_sd = 0.2, psi_mean = 0, psi_sd = 1)
  design <- escalation::dont_skip_doses(design)
  sims <- escalation::simulate_trials(
    design, num_sims = R, true_prob_eff = sc$pE, true_prob_tox = sc$pT,
    sample_sizes = rep(cohort, floor(sc$N / cohort)))
  pr  <- escalation::prob_recommend(sims)
  none <- as.numeric(pr["NoDose"]); if (is.na(none)) none <- as.numeric(pr["0"])
  pcs <- if (!is.na(sc$true_opt)) as.numeric(pr[as.character(sc$true_opt)]) else NA
  list(PCS = 100 * pcs, P_none = 100 * none)
}

## ============================================================
## DRIVER: build Table 3
## ============================================================
J <- 5L; rows <- list()
for (nm in names(scenarios)) {
  sc <- scenarios[[nm]]

  ## (1) BART arm -- your validated engine; should reproduce Table 2
  bart <- summarise_ocs(simulate_design(
            sc, R = R_COMP, n_cores = N_CORES, seed = MASTER_SEED, chunk_size = 50L,
            max_n = sc$N, cohort = DESIGN$cohort, phi_T = DESIGN$phi_T,
            phi_E = DESIGN$phi_E, c_T_pre = DESIGN$c_T_pre, c_T = DESIGN$c_T,
            c_E = DESIGN$c_E, a_T = DESIGN$a_T, b_T = DESIGN$b_T, M = DESIGN$M,
            n_iter = DESIGN$n_iter, n_burn = DESIGN$n_burn, use_covars = sc$covars),
            sc$true_opt, sc$true_adm, J)

  ## (2) Parametric-utility arm -- same conduct, BART removed
  para <- summarise_ocs(simulate_arm(
            sc, parametric_interim, R = R_COMP, n_cores = N_CORES,
            seed = MASTER_SEED, chunk_size = 50L,
            max_n = sc$N, cohort = DESIGN$cohort, phi_T = DESIGN$phi_T,
            phi_E = DESIGN$phi_E, c_T_pre = DESIGN$c_T_pre, c_T = DESIGN$c_T,
            c_E = DESIGN$c_E, a_T = DESIGN$a_T, b_T = DESIGN$b_T,
            n_iter = DESIGN$n_iter, n_burn = DESIGN$n_burn, use_covars = sc$covars),
            sc$true_opt, sc$true_adm, J)

  ## (3) EffTox arm -- comment in once trialr/escalation + contour set
  eff <- list(PCS = NA, P_none = NA)
  # eff <- efftox_ocs(sc)

  has_opt <- !is.na(sc$true_opt)
  mcse <- function(p) round(100 * sqrt((p/100) * (1 - p/100) / R_COMP), 1)
  rows[[nm]] <- data.frame(
    Scenario = nm, Name = sc$name,
    BART      = round(100 * bart$PCS, 1),
    Parametric= round(100 * para$PCS, 1),
    EffTox    = round(eff$PCS, 1),
    BART_none = round(100 * bart$P_none, 1),
    Para_none = round(100 * para$P_none, 1),
    metric    = if (has_opt) "PCS" else "P_none",
    mcse_BART = mcse(100 * bart$PCS),
    stringsAsFactors = FALSE)
}
table3 <- do.call(rbind, rows)
cat("\n==================== Table 3 ====================\n"); print(table3)
write.csv(table3, "table3_comparison.csv", row.names = FALSE)

## ---- validation check: BART column vs published Table 2 ----
published <- c(S1 = 12.1, S2 = 39.0, S3 = 36.4, S4 = 24.5,
               S5 = 33.0, S6 = 30.0, S7 = 30.9)
chk <- merge(data.frame(Scenario = names(published), Table2 = published),
             table3[, c("Scenario", "BART")], by = "Scenario")
chk$diff <- round(chk$BART - chk$Table2, 1)
cat("\nBART arm vs published Table 2 (|diff| should be within MC error):\n")
print(chk)
if (any(abs(chk$diff) > 4, na.rm = TRUE))
  cat("\nWARNING: BART arm does not reproduce Table 2; fix wiring before trusting Table 3.\n")

## ---- LaTeX for manuscript Table 3 (PCS columns; fills the dots) ----
con <- file("table3_comparison.tex", "w")
cat("% replace the placeholder body of Table~\\ref{tab:comparison} with these rows\n",
    file = con)
for (i in seq_len(nrow(table3))) with(table3[i, ],
  cat(sprintf("%s & %s & %.1f & %.1f & %s \\\\\n", Scenario, Name, BART, Parametric,
              ifelse(is.na(EffTox), "$\\cdot$", sprintf("%.1f", EffTox))), file = con))
close(con)
cat("\nWrote table3_comparison.csv and table3_comparison.tex.\n")
