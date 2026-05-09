# ============================================================
# R/bart_efftox.R
# Core functions for the BART-EffTox Phase I/II dose-finding design.
# All functions correspond directly to the manuscript specification
# (Sections 3-4). See vignette("bart-efftox-intro") for usage.
# ============================================================

#' Bilinear efficacy-toxicity utility
#'
#' Computes the utility score \eqn{U(p_E, p_T)} used for dose selection
#' in the BART-EffTox design (manuscript equation~1). The utility is
#' bilinear in the outcome probabilities with four corner values.
#'
#' @param pE Numeric. Efficacy probability (scalar or vector).
#' @param pT Numeric. Toxicity probability (same length as \code{pE}).
#' @param u10 Numeric. Corner value: efficacy without toxicity (default 100).
#' @param u11 Numeric. Corner value: efficacy with toxicity (default 40).
#' @param u00 Numeric. Corner value: neither outcome (default 20).
#' @param u01 Numeric. Corner value: toxicity without efficacy (default 0).
#'
#' @return Numeric vector of utility values.
#'
#' @references
#' Thall PF, Cook JD (2004). Dose-finding based on efficacy-toxicity
#' trade-offs. \emph{Biometrics} \strong{60}(3):684--693.
#'
#' @export
#' @examples
#' utility_fn(pE = 0.55, pT = 0.25)   # d4 in Scenario S1
#' utility_fn(pE = c(0.10, 0.25, 0.40, 0.55, 0.65),
#'            pT = c(0.05, 0.10, 0.18, 0.25, 0.45))
utility_fn <- function(pE, pT, u10 = 100, u11 = 40, u00 = 20, u01 = 0) {
  u10 * pE * (1 - pT) +
  u11 * pE * pT +
  u00 * (1 - pE) * (1 - pT) +
  u01 * (1 - pE) * pT
}


#' Weighted PAVA isotonic projection
#'
#' Projects a vector onto the monotone cone
#' \eqn{\mathcal{C} = \{q \in [0,1]^J : q_1 \le \cdots \le q_J\}}
#' using the weighted pool-adjacent-violators algorithm (PAVA).
#' Applied to each posterior toxicity draw before admissibility and
#' utility calculations (Proposition~3, Section~4.3 of the manuscript).
#'
#' @param p Numeric vector of length \eqn{J}. Raw posterior toxicity
#'   probabilities (or any vector to be isotonically projected).
#' @param w Numeric vector of weights of length \eqn{J}.
#'   Default \code{NULL} gives equal weights of 1.
#'   In the design, \code{w[j] = n_j + 1} where \code{n_j} is the
#'   number of patients observed at dose \eqn{d_j}.
#'
#' @return Numeric vector of the same length as \code{p}, non-decreasing,
#'   with values in \eqn{[0,1]}.
#'
#' @details
#' Properties (Proposition~3):
#' \enumerate{
#'   \item \strong{Order enforcement}: output is non-decreasing.
#'   \item \strong{Idempotence}: already-monotone inputs are unchanged.
#'   \item \strong{Nonexpansiveness}: \eqn{\|\tilde{p} - p^0\|_W \le \|p - p^0\|_W}
#'     for any monotone \eqn{p^0}.
#'   \item \strong{Blockwise averaging}: violations resolved by merging
#'     adjacent blocks and replacing with the weighted mean.
#' }
#' No isotonic constraint is applied to efficacy probabilities.
#'
#' @references
#' Robertson T, Wright FT, Dykstra RL (1988).
#' \emph{Order Restricted Statistical Inference}. Wiley, New York.
#'
#' @export
#' @examples
#' # Non-monotone input gets projected
#' pava_project(c(0.05, 0.20, 0.15, 0.30, 0.45))
#' # Already monotone: returned unchanged
#' pava_project(c(0.05, 0.10, 0.18, 0.25, 0.45))
#' # With data-adaptive weights
#' pava_project(c(0.05, 0.20, 0.15, 0.30, 0.45), w = c(4, 8, 6, 3, 3))
pava_project <- function(p, w = NULL) {
  J <- length(p)
  if (is.null(w)) w <- rep(1, J)
  r <- p
  for (pass in seq_len(100L)) {
    changed <- FALSE
    for (i in seq_len(J - 1L)) {
      if (r[i] > r[i + 1L] + 1e-12) {
        avg      <- (w[i] * r[i] + w[i + 1L] * r[i + 1L]) / (w[i] + w[i + 1L])
        r[i]     <- avg
        r[i + 1L]<- avg
        changed  <- TRUE
      }
    }
    if (!changed) break
  }
  pmin(pmax(r, 0), 1)
}


#' Beta-binomial run-in safety rule
#'
#' Computes the posterior tail probability
#' \eqn{\Pr(p_j > \phi_T \mid y_j, n_j)} under a
#' \eqn{\mathrm{Beta}(a_T, b_T)} conjugate prior, and returns
#' \code{TRUE} if escalation from dose \eqn{d_j} should be blocked
#' (i.e., the tail exceeds the cutoff \eqn{c_{T,\mathrm{pre}}}).
#' This implements the run-in safety rule of Proposition~4.
#'
#' @param y_j Integer. Observed toxicity count at dose \eqn{d_j}.
#' @param n_j Integer. Total patients at dose \eqn{d_j}.
#' @param phi_T Numeric. Toxicity threshold (default 0.30).
#' @param c_T_pre Numeric. Run-in blocking cutoff (default 0.70).
#' @param a_T Numeric. Beta prior shape parameter (default 1, uniform).
#' @param b_T Numeric. Beta prior shape parameter (default 1, uniform).
#'
#' @return Logical. \code{TRUE} if escalation is blocked.
#'
#' @references
#' Manuscript Proposition~4 and Table~1.
#'
#' @export
#' @examples
#' runin_admissibility(y_j = 0, n_j = 3)   # 0/3 tox -> not blocked
#' runin_admissibility(y_j = 2, n_j = 3)   # 2/3 tox -> blocked (c_T_pre=0.70)
#' # Posterior tail value
#' 1 - pbeta(0.30, 1 + 2, 1 + 3 - 2)       # ~0.916
runin_admissibility <- function(y_j, n_j,
                                phi_T   = 0.30,
                                c_T_pre = 0.70,
                                a_T     = 1L,
                                b_T     = 1L) {
  tail_prob <- 1 - pbeta(phi_T, a_T + y_j, b_T + n_j - y_j)
  tail_prob > c_T_pre
}


#' Single BART-EffTox interim analysis
#'
#' Fits probit BART models for toxicity and efficacy, applies the PAVA
#' isotonic projection to the toxicity posterior, computes posterior
#' admissibility and posterior mean utility for each dose, and selects
#' the next dose (or stops the trial). Implements the BART-adaptive
#' phase decision rule described in Sections 3.3--3.4 of the manuscript.
#'
#' @param dose_obs Integer vector. Dose levels (1 to \code{J}) for each
#'   enrolled patient.
#' @param yT Integer vector. Binary toxicity outcomes (0/1), same length
#'   as \code{dose_obs}.
#' @param yE Integer vector. Binary efficacy outcomes (0/1), same length
#'   as \code{dose_obs}.
#' @param X_cov Numeric matrix or \code{NULL}. Patient-level covariates,
#'   one row per patient. If \code{NULL}, no covariate adjustment.
#' @param dose_grid Integer vector. Dose levels to evaluate (default 1:5).
#' @param X_ref Numeric matrix or \code{NULL}. Reference covariate values
#'   for marginalizing predictions. If \code{NULL}, uses \code{X_cov}.
#' @param current_dose Integer. Current dose level (for no-skip rule).
#' @param phi_T Numeric. Toxicity threshold (default 0.30).
#' @param phi_E Numeric. Efficacy threshold (default 0.20).
#' @param c_T Numeric. Toxicity admissibility cutoff (default 0.40).
#' @param c_E Numeric. Efficacy admissibility cutoff (default 0.50).
#' @param M Integer. Number of BART trees (default 20).
#' @param n_iter Integer. Posterior draws (default 400).
#' @param n_burn Integer. Burn-in draws (default 80).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{next_dose}}{Integer or \code{NA}. Selected next dose.}
#'   \item{\code{stop_trial}}{Logical. \code{TRUE} if no admissible dose found.}
#'   \item{\code{admissible_set}}{Integer vector. Posterior-admissible dose indices.}
#'   \item{\code{mean_utility}}{Numeric vector length \eqn{J}. Posterior mean utilities.}
#'   \item{\code{pr_toxic}}{Numeric vector. \eqn{\Pr(\tilde p_T > \phi_T \mid \mathcal{D}_n)}.}
#'   \item{\code{pr_ineffic}}{Numeric vector. \eqn{\Pr(p_E < \phi_E \mid \mathcal{D}_n)}.}
#'   \item{\code{post_pT}}{Matrix \eqn{S \times J}. PAVA-projected toxicity draws.}
#'   \item{\code{post_pE}}{Matrix \eqn{S \times J}. Efficacy posterior draws.}
#' }
#'
#' @importFrom BART pbart
#' @export
bart_efftox_interim <- function(dose_obs,
                                yT,
                                yE,
                                X_cov        = NULL,
                                dose_grid    = 1:5,
                                X_ref        = NULL,
                                current_dose = 1L,
                                phi_T        = 0.30,
                                phi_E        = 0.20,
                                c_T          = 0.40,
                                c_E          = 0.50,
                                M            = 20L,
                                n_iter       = 400L,
                                n_burn       = 80L) {

  J            <- length(dose_grid)
  dose_std     <- 2 * (dose_obs  - 1L) / (J - 1L) - 1
  dose_std_grid<- 2 * (dose_grid - 1L) / (J - 1L) - 1

  # Build predictor matrices
  X_train <- if (!is.null(X_cov)) cbind(dose_std, X_cov) else
               matrix(dose_std, ncol = 1L)

  if (!is.null(X_cov) && !is.null(X_ref)) {
    n_ref  <- nrow(X_ref)
    X_pred <- do.call(rbind,
                lapply(dose_std_grid, function(d) cbind(d, X_ref)))
    marginal <- TRUE
  } else {
    X_pred   <- matrix(dose_std_grid, ncol = 1L)
    n_ref    <- 1L
    marginal <- FALSE
  }

  # Fit probit BART for toxicity and efficacy (Albert-Chib augmentation)
  fit_T <- BART::pbart(X_train, yT, X_pred,
                       ntree    = M,
                       ndpost   = n_iter,
                       nskip    = n_burn,
                       printevery = 0L)
  fit_E <- BART::pbart(X_train, yE, X_pred,
                       ntree    = M,
                       ndpost   = n_iter,
                       nskip    = n_burn,
                       printevery = 0L)

  pT_raw <- fit_T$prob.test   # S x (J * n_ref)
  pE_raw <- fit_E$prob.test

  # Marginalize over reference covariates if present
  if (marginal) {
    pT_samp <- pE_samp <- matrix(NA_real_, n_iter, J)
    for (j in seq_len(J)) {
      idx         <- ((j - 1L) * n_ref + 1L):(j * n_ref)
      pT_samp[, j]<- rowMeans(pT_raw[, idx, drop = FALSE])
      pE_samp[, j]<- rowMeans(pE_raw[, idx, drop = FALSE])
    }
  } else {
    pT_samp <- pT_raw
    pE_samp <- pE_raw
  }

  # PAVA isotonic projection on each toxicity draw (Proposition 3)
  n_j    <- tabulate(dose_obs, nbins = J)
  w_j    <- n_j + 1L
  pT_proj<- matrix(NA_real_, n_iter, J)
  for (s in seq_len(n_iter))
    pT_proj[s, ] <- pava_project(pT_samp[s, ], w_j)

  # Posterior admissibility (equation 2)
  pr_toxic  <- colMeans(pT_proj > phi_T)
  pr_ineffic<- colMeans(pE_samp < phi_E)
  adm <- (pr_toxic < c_T) & (pr_ineffic < c_E)

  # No-skip escalation: restrict to current dose +/- 1
  safety_flag <- mean(pT_proj[, current_dose] > phi_T) > c_T
  eligible    <- which(adm & (seq_len(J) <= current_dose + 1L))
  if (safety_flag) eligible <- eligible[eligible <= current_dose]

  if (length(eligible) == 0L)
    return(list(next_dose = NA_integer_, stop_trial = TRUE,
                admissible_set = integer(0L),
                mean_utility   = rep(NA_real_, J),
                post_pT = pT_proj, post_pE = pE_samp,
                pr_toxic = round(pr_toxic, 3L),
                pr_ineffic = round(pr_ineffic, 3L)))

  # Posterior mean utility (plug-in marginal; equation 3)
  mean_util <- numeric(J)
  for (j in seq_len(J))
    mean_util[j] <- mean(utility_fn(pE_samp[, j], pT_proj[, j]))

  # Select dose with highest utility; tie-break: prefer current dose
  best_util <- max(mean_util[eligible])
  tied      <- eligible[abs(mean_util[eligible] - best_util) < 1e-6]
  best      <- if (current_dose %in% tied) current_dose else min(tied)

  list(next_dose      = best,
       stop_trial     = FALSE,
       admissible_set = eligible,
       mean_utility   = round(mean_util, 3L),
       post_pT        = pT_proj,
       post_pE        = pE_samp,
       pr_toxic       = round(pr_toxic, 3L),
       pr_ineffic     = round(pr_ineffic, 3L))
}


#' Simulate a single BART-EffTox trial
#'
#' Runs one complete two-phase trial (balanced run-in followed by
#' BART-adaptive phase) under a given scenario. Implements the trial
#' conduct described in Section~3.5 of the manuscript.
#'
#' @param scenario A list with fields:
#'   \code{pT} (length-J toxicity probabilities),
#'   \code{pE} (length-J efficacy probabilities),
#'   \code{covars} (logical, whether patient covariates are simulated),
#'   and optional covariate-specific fields (\code{pB}, \code{pF},
#'   \code{pT_xb0}, \code{pT_xb1}, \code{pE_xb0}, \code{pE_xb1},
#'   \code{pT_xf0}, \code{pT_xf1}, \code{pE_xf0}, \code{pE_xf1}).
#' @param max_n Integer. Maximum sample size \eqn{N}.
#' @param cohort Integer. Cohort size \eqn{m} (default 3).
#' @param phi_T,phi_E,c_T_pre,c_T,c_E Numeric. Design parameters.
#' @param a_T,b_T Numeric. Beta prior shape parameters (default 1).
#' @param M,n_iter,n_burn Integer. BART tuning parameters.
#' @param use_covars Logical. Whether to use covariate information.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{rec_dose}}{Integer or \code{NA}. Final recommended dose.}
#'   \item{\code{stopped_early}}{Logical.}
#'   \item{\code{n}}{Integer. Total patients enrolled.}
#'   \item{\code{n_at_dose}}{Integer vector length \eqn{J}.}
#'   \item{\code{mean_tox}}{Numeric. Observed toxicity rate.}
#'   \item{\code{mean_eff}}{Numeric. Observed efficacy rate.}
#' }
#'
#' @export
run_trial <- function(scenario,
                      max_n    = 36L,
                      cohort   = 3L,
                      phi_T    = 0.30,
                      phi_E    = 0.20,
                      c_T_pre  = 0.70,
                      c_T      = 0.40,
                      c_E      = 0.50,
                      a_T      = 1L,
                      b_T      = 1L,
                      M        = 20L,
                      n_iter   = 400L,
                      n_burn   = 80L,
                      use_covars = FALSE) {

  J           <- length(scenario$pT)
  dose_obs    <- integer(0)
  yT          <- integer(0)
  yE          <- integer(0)
  X_cov_obs   <- NULL
  current     <- 1L
  bart_active <- FALSE

  # Helper: generate one cohort
  gen_cohort <- function(d, m) {
    if (isTRUE(scenario$covars) && use_covars) {
      xb   <- if (!is.null(scenario$pB)) rbinom(m, 1L, scenario$pB) else NULL
      xf   <- if (!is.null(scenario$pF)) rbinom(m, 1L, scenario$pF) else NULL
      pT_i <- if (!is.null(xf))
                ifelse(xf == 1L, scenario$pT_xf1[d], scenario$pT_xf0[d])
              else scenario$pT[d]
      pE_i <- if (!is.null(xb))
                ifelse(xb == 1L, scenario$pE_xb1[d], scenario$pE_xb0[d])
              else scenario$pE[d]
      cols  <- Filter(Negate(is.null), list(xb = xb, xf = xf))
      X_new <- if (length(cols) > 0) do.call(cbind, cols) else NULL
    } else {
      pT_i  <- scenario$pT[d]
      pE_i  <- scenario$pE[d]
      X_new <- NULL
    }
    list(yT = rbinom(m, 1L, pT_i),
         yE = rbinom(m, 1L, pE_i),
         d  = rep(d, m),
         X  = X_new)
  }

  # Helper: beta-binomial safety check
  dose_safe <- function(d) {
    nj <- sum(dose_obs == d)
    tj <- sum(yT[dose_obs == d])
    if (nj == 0L) return(TRUE)
    !runin_admissibility(tj, nj, phi_T, c_T_pre, a_T, b_T)
  }

  # Phase 1: balanced run-in
  for (d in seq_len(J)) {
    if (length(dose_obs) >= max_n) break
    nxt      <- gen_cohort(d, cohort)
    dose_obs <- c(dose_obs, nxt$d)
    yT       <- c(yT, nxt$yT)
    yE       <- c(yE, nxt$yE)
    if (!is.null(nxt$X)) X_cov_obs <- rbind(X_cov_obs, nxt$X)
    current  <- d
    if (!dose_safe(d)) { current <- max(1L, d - 1L); bart_active <- TRUE; break }
    if (d == J) bart_active <- TRUE
  }
  if (length(dose_obs) >= max_n) bart_active <- TRUE

  rec_dose     <- NA_integer_
  stopped_early<- FALSE

  # Phase 2: BART-adaptive
  if (bart_active) {
    while (length(dose_obs) < max_n) {
      X_ref <- if (use_covars && !is.null(X_cov_obs)) X_cov_obs else NULL
      res   <- bart_efftox_interim(
        dose_obs = dose_obs, yT = yT, yE = yE,
        X_cov = X_cov_obs, dose_grid = seq_len(J),
        X_ref = X_ref, current_dose = current,
        phi_T = phi_T, phi_E = phi_E,
        c_T = c_T, c_E = c_E,
        M = M, n_iter = n_iter, n_burn = n_burn)
      if (res$stop_trial) { stopped_early <- TRUE; break }
      current  <- res$next_dose
      nxt      <- gen_cohort(current, cohort)
      dose_obs <- c(dose_obs, nxt$d)
      yT       <- c(yT, nxt$yT)
      yE       <- c(yE, nxt$yE)
      if (!is.null(nxt$X)) X_cov_obs <- rbind(X_cov_obs, nxt$X)
    }
  }

  # Final recommendation
  if (!stopped_early && length(dose_obs) >= max_n) {
    X_ref <- if (use_covars && !is.null(X_cov_obs)) X_cov_obs else NULL
    final <- bart_efftox_interim(
      dose_obs = dose_obs, yT = yT, yE = yE,
      X_cov = X_cov_obs, dose_grid = seq_len(J),
      X_ref = X_ref, current_dose = current,
      phi_T = phi_T, phi_E = phi_E,
      c_T = c_T, c_E = c_E,
      M = M, n_iter = n_iter, n_burn = n_burn)
    if (!final$stop_trial) rec_dose <- final$next_dose
  }

  list(rec_dose     = rec_dose,
       stopped_early= stopped_early,
       n            = length(dose_obs),
       n_at_dose    = tabulate(dose_obs, nbins = J),
       mean_tox     = mean(yT),
       mean_eff     = mean(yE))
}


#' Simulate multiple BART-EffTox trials in parallel
#'
#' Runs \code{R} independent trials using \code{parallel::mclapply}
#' with chunked execution for macOS stability. Each chunk uses a
#' separate random seed derived from \code{seed}.
#'
#' @param scenario A scenario list as described in \code{\link{run_trial}}.
#' @param R Integer. Number of trials (default 2000).
#' @param n_cores Integer. Parallel workers (default 2; macOS: keep <=4).
#' @param seed Integer. Master seed for reproducibility.
#' @param chunk_size Integer. Trials per mclapply chunk (default 50).
#' @param ... Additional arguments passed to \code{\link{run_trial}}.
#'
#' @return A list of length \code{R}. Each element is either the return
#'   value of \code{\link{run_trial}} or a \code{trial_error} object if
#'   the trial threw an error.
#'
#' @export
simulate_design <- function(scenario,
                            R          = 2000L,
                            n_cores    = 2L,
                            seed       = 42L,
                            chunk_size = 50L,
                            ...) {
  set.seed(seed)
  seeds  <- sample.int(1e6L, R)
  chunks <- split(seq_len(R), ceiling(seq_len(R) / chunk_size))
  all_results <- vector("list", R)

  for (ch in chunks) {
    ch_res <- parallel::mclapply(ch, function(r) {
      set.seed(seeds[r])
      tryCatch(
        run_trial(scenario, ...),
        error = function(e) structure(
          list(message = conditionMessage(e), seed = seeds[r]),
          class = "trial_error"))
    }, mc.cores = n_cores)
    for (i in seq_along(ch)) all_results[[ch[i]]] <- ch_res[[i]]
  }
  all_results
}


#' Summarize operating characteristics from simulation results
#'
#' Computes operating characteristics from a list of trial simulation
#' results, as reported in Table~2 and Table~3 of the manuscript.
#' All probabilities are marginal over all \code{R} trials.
#'
#' @param results List. Output of \code{\link{simulate_design}}.
#' @param true_opt Integer. Index of the true optimal dose.
#' @param true_adm Logical vector of length \eqn{J}. \code{TRUE} for
#'   truly admissible doses.
#' @param J Integer. Number of dose levels (default 5).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{R}}{Number of valid (non-failed) trials.}
#'   \item{\code{n_fail}}{Number of failed trials.}
#'   \item{\code{PCS}}{Marginal probability of correct selection.}
#'   \item{\code{P_inadm}}{Marginal probability of inadmissible recommendation.}
#'   \item{\code{P_early}}{Probability of early stopping before \eqn{N}.}
#'   \item{\code{P_none}}{Probability of no final recommendation.}
#'   \item{\code{mean_N}}{Mean patients enrolled.}
#'   \item{\code{sel_pct}}{Marginal dose-specific recommendation percentages.}
#'   \item{\code{n_dose}}{Mean patients per dose.}
#'   \item{\code{mean_tox}}{Mean observed toxicity rate.}
#'   \item{\code{mean_eff}}{Mean observed efficacy rate.}
#' }
#'
#' @export
summarise_ocs <- function(results, true_opt, true_adm, J = 5L) {
  bad    <- sapply(results,
                   function(x) is.null(x) || inherits(x, "trial_error"))
  n_fail <- sum(bad)
  if (n_fail > 0)
    warning(sprintf("%d/%d trials failed or returned NULL",
                    n_fail, length(results)))
  results <- results[!bad]

  if (length(results) == 0L) {
    warning("No valid trials after filtering.")
    return(list(R = 0L, n_fail = n_fail,
                PCS = NA, P_inadm = NA,
                P_early = NA, P_none = NA,
                mean_N = NA,
                sel_pct = rep(NA, J), n_dose = rep(NA, J),
                mean_tox = NA, mean_eff = NA))
  }

  rec     <- as.integer(sapply(results, function(x) {
    v <- x[["rec_dose"]]
    if (is.null(v)) NA_integer_ else as.integer(v) }))
  stopped <- as.logical(sapply(results, `[[`, "stopped_early"))
  n_tot   <- as.numeric(sapply(results, `[[`, "n"))
  n_at    <- do.call(rbind, lapply(results, `[[`, "n_at_dose"))
  tox_r   <- as.numeric(sapply(results, `[[`, "mean_tox"))
  eff_r   <- as.numeric(sapply(results, `[[`, "mean_eff"))

  R_valid <- length(results)

  list(
    R        = R_valid,
    n_fail   = n_fail,
    PCS      = mean(rec == true_opt,              na.rm = TRUE),
    P_inadm  = mean(rec %in% which(!true_adm),    na.rm = TRUE),
    P_early  = mean(stopped,                      na.rm = TRUE),
    P_none   = mean(is.na(rec)),
    mean_N   = mean(n_tot, na.rm = TRUE),
    sel_pct  = round(100 * tabulate(rec[!is.na(rec)], nbins = J) / R_valid, 1),
    n_dose   = if (is.matrix(n_at)) round(colMeans(n_at), 1) else rep(NA, J),
    mean_tox = round(mean(tox_r, na.rm = TRUE), 3),
    mean_eff = round(mean(eff_r, na.rm = TRUE), 3)
  )
}
