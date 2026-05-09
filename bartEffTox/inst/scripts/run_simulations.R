# ============================================================
# run_sims_100.R
# BART-EffTox: 100-trial pilot simulation across all 7 scenarios
#
# Usage: Rscript run_sims_100.R
#   or:  source("run_sims_100.R")
#
# Output: sims_100_results.rds  -- raw trial results (list of lists)
#         sims_100_summary.csv  -- OC table ready for manuscript
#         sims_100_summary.txt  -- formatted console output
#
# Runtime estimate (BART package, 4 cores):
#   ~0.5 sec/trial x 100 trials x 7 scenarios = ~6 min
# ============================================================

# ---- DEPENDENCY CHECK ----
if (!requireNamespace("BART", quietly=TRUE)) {
  stop(paste0(
    "\n\nBART package required. Install with: install.packages('BART')\n",
    "The self-contained stump sampler produces INVALID OCs (80-100% stopping).\n"
  ))
}

# ---- ARCHITECTURE FILE CHECK ----
# Ensure bart_efftox_final_architecture_v2.R is the CURRENT version
# with the pre-BART balanced exploration phase.
# If you see 80-100% stopping rates, you have an old version of this file.
arch_file <- "bart_efftox_final_architecture_v2.R"
if (!file.exists(arch_file)) {
  stop(paste0("Cannot find ", arch_file,
              " -- place it in the same directory as this script."))
}

# Check version before sourcing
arch_lines <- readLines(arch_file)
has_prebart <- any(grepl("bart_active|dose_safe|c_T_pre", arch_lines))
if (!has_prebart) {
  stop(paste0(
    "\n\nOLD version of ", arch_file, " detected!\n",
    "The pre-BART balanced exploration phase is missing.\n",
    "Download the current version from your outputs folder.\n",
    "Key indicator: the file must contain 'bart_active' and 'c_T_pre'.\n"
  ))
}
source(arch_file)

# Final check: version sentinel
if (!exists("BART_EFFTOX_ARCH_VERSION")) {
  warning("Architecture file loaded but version sentinel not found.",
          " Ensure you have the latest version.")
} else {
  cat(sprintf("Architecture version: %s\n", BART_EFFTOX_ARCH_VERSION))
}

# ============================================================
# SCENARIO DEFINITIONS (Supplementary Table S1, v3)
# ============================================================
scenarios <- list(

  S1 = list(
    name   = "Parametric benchmark",
    pT     = c(0.05, 0.10, 0.18, 0.25, 0.45),
    pE     = c(0.10, 0.25, 0.40, 0.55, 0.65),
    covars = FALSE, N = 36L, true_opt = 4L,
    true_adm = c(FALSE, TRUE, TRUE, TRUE, FALSE)
  ),

  S2 = list(
    name   = "Plateau efficacy",
    pT     = c(0.05, 0.11, 0.30, 0.45, 0.62),
    pE     = c(0.15, 0.42, 0.48, 0.50, 0.51),
    covars = FALSE, N = 36L, true_opt = 2L,
    true_adm = c(FALSE, TRUE, TRUE, FALSE, FALSE)
  ),

  S3 = list(
    name   = "Non-monotone efficacy",
    pT     = c(0.04, 0.09, 0.17, 0.28, 0.48),
    pE     = c(0.10, 0.28, 0.50, 0.38, 0.22),
    covars = FALSE, N = 36L, true_opt = 3L,
    true_adm = c(FALSE, TRUE, TRUE, TRUE, FALSE)
  ),

  S4 = list(
    name    = "Biomarker-modified efficacy",
    pB      = 0.40,
    pT_xb0  = c(0.05, 0.10, 0.18, 0.32, 0.50),
    pT_xb1  = c(0.05, 0.10, 0.18, 0.32, 0.50),
    pE_xb0  = c(0.08, 0.15, 0.22, 0.25, 0.25),
    pE_xb1  = c(0.15, 0.35, 0.60, 0.70, 0.72),
    covars  = TRUE, N = 36L, true_opt = 3L,
    true_adm = c(FALSE, TRUE, TRUE, FALSE, FALSE),
    # Marginal pT and pE for admissibility check
    pT = c(0.05, 0.10, 0.18, 0.32, 0.50),
    pE = 0.60*c(0.08,0.15,0.22,0.25,0.25) + 0.40*c(0.15,0.35,0.60,0.70,0.72)
  ),

  S5 = list(
    name    = "Frailty-modified toxicity",
    pF      = 0.30,
    pT_xf0  = c(0.02, 0.05, 0.10, 0.19, 0.33),
    pT_xf1  = c(0.08, 0.20, 0.38, 0.58, 0.78),
    pE_xf0  = c(0.10, 0.28, 0.46, 0.55, 0.56),
    pE_xf1  = c(0.10, 0.28, 0.46, 0.55, 0.56),
    covars  = TRUE, N = 36L, true_opt = 3L,
    true_adm = c(FALSE, TRUE, TRUE, FALSE, FALSE),
    pT = 0.70*c(0.02,0.05,0.10,0.19,0.33) + 0.30*c(0.08,0.20,0.38,0.58,0.78),
    pE = c(0.10, 0.28, 0.46, 0.55, 0.56)
  ),

  S6 = list(
    name   = "Small-sample stress test",
    pT     = c(0.06, 0.14, 0.28, 0.42, 0.58),
    pE     = c(0.14, 0.38, 0.42, 0.40, 0.32),
    covars = FALSE, N = 24L, true_opt = 2L,
    true_adm = c(FALSE, TRUE, TRUE, FALSE, FALSE)
  ),

  S7 = list(
    name    = "Combined modern oncology",
    pB      = 0.40, pF = 0.30,
    pT_xf0  = c(0.03, 0.07, 0.14, 0.24, 0.42),
    pT_xf1  = c(0.07, 0.16, 0.30, 0.50, 0.70),
    pE_xb0  = c(0.08, 0.18, 0.28, 0.25, 0.18),
    pE_xb1  = c(0.15, 0.40, 0.65, 0.55, 0.38),
    covars  = TRUE, N = 36L, true_opt = 3L,
    true_adm = c(FALSE, TRUE, TRUE, FALSE, FALSE),
    pT = 0.70*c(0.03,0.07,0.14,0.24,0.42) + 0.30*c(0.07,0.16,0.30,0.50,0.70),
    pE = 0.60*c(0.08,0.18,0.28,0.25,0.18) + 0.40*c(0.15,0.40,0.65,0.55,0.38)
  )
)

# ============================================================
# SIMULATION SETTINGS
# ============================================================
R_pilot  <- 2000L         # full R=2000 run
N_CORES  <- 4L            # 4 cores; adjust to your machine
MASTER_SEED <- 2025L      # master seed for reproducibility

# BART-EffTox design parameters (Section 3, matching Section 4 notation)
# Phase 1 — Pre-BART balanced exploration (new design):
#   c_T_pre=0.70: permissive safety check during pre-BART phase
#   Requires 2+ tox events in 3 patients to block a dose (Table 1 updated)
#   Ensures BART always receives data at 3-5 dose levels before activating
# Phase 2 — BART-EffTox adaptive selection:
#   c_T=0.40, c_E=0.50: relaxed BART-phase admissibility cutoffs
DESIGN <- list(
  cohort  = 3L,
  phi_T   = 0.30,   # toxicity threshold (both phases)
  phi_E   = 0.20,   # efficacy threshold (both phases)
  c_T_pre = 0.70,   # pre-BART safety cutoff (permissive, Prop 4 updated)
  c_T     = 0.40,   # BART-phase toxicity admissibility cutoff
  c_E     = 0.50,   # BART-phase efficacy admissibility cutoff
  a_T     = 1L,     # beta prior shape (uniform)
  b_T     = 1L,     # beta prior shape (uniform)
  M       = 20L,    # BART trees
  n_iter  = 400L,   # 400 posterior draws (validated against n_iter=1250)
  n_burn  = 80L     # burn-in (validated)
)

# ============================================================
# RUN SIMULATIONS
# ============================================================
cat("============================================================\n")
cat(" BART-EffTox simulation: R=2000, 7 scenarios\n")
cat("============================================================\n")
cat(sprintf(" Cores: %d | Seed: %d | M: %d\n",
    N_CORES, MASTER_SEED, DESIGN$M))
cat(sprintf(" phi_T=%.2f phi_E=%.2f c_T_pre=%.2f c_T=%.2f c_E=%.2f\n\n",
    DESIGN$phi_T, DESIGN$phi_E, DESIGN$c_T_pre, DESIGN$c_T, DESIGN$c_E))

all_results <- list()
timing      <- numeric(length(scenarios))

for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]

  cat(sprintf("[%s] %s (N=%d, true opt=d%d)... ",
      sc_name, sc$name, sc$N, sc$true_opt))
  flush.console()

  t0 <- proc.time()

  res <- simulate_design(
    scenario   = sc,
    R          = R_pilot,
    n_cores    = N_CORES,
    seed       = MASTER_SEED,
    chunk_size = 50L,          # macOS: keep <=50; Linux: can use 200+
    # design parameters
    max_n      = sc$N,
    cohort     = DESIGN$cohort,
    phi_T      = DESIGN$phi_T,
    phi_E      = DESIGN$phi_E,
    c_T_pre    = DESIGN$c_T_pre,
    c_T        = DESIGN$c_T,
    c_E        = DESIGN$c_E,
    a_T        = DESIGN$a_T,
    b_T        = DESIGN$b_T,
    M          = DESIGN$M,
    n_iter     = DESIGN$n_iter,
    n_burn     = DESIGN$n_burn,
    use_covars = sc$covars
  )

  elapsed      <- (proc.time() - t0)["elapsed"]
  timing[sc_name] <- elapsed
  all_results[[sc_name]] <- res

  cat(sprintf("done (%.0f sec, %.1f sec/trial)\n",
      elapsed, elapsed / R_pilot))
}

cat(sprintf("\nTotal runtime: %.1f sec (%.1f min)\n\n",
    sum(timing), sum(timing)/60))

# ============================================================
# SUMMARISE OPERATING CHARACTERISTICS
# ============================================================
cat("============================================================\n")
cat(" OPERATING CHARACTERISTICS SUMMARY\n")
cat("============================================================\n\n")

J <- 5L
oc_list <- list()

for (sc_name in names(scenarios)) {
  sc  <- scenarios[[sc_name]]
  res <- all_results[[sc_name]]
  oc  <- summarise_ocs(res,
                        true_opt = sc$true_opt,
                        true_adm = sc$true_adm,
                        J        = J)
  oc_list[[sc_name]] <- oc

  cat(sprintf("--- %s: %s ---\n", sc_name, sc$name))
  cat(sprintf("  R (completed):    %5d\n",             oc$R))
  if (oc$n_fail > 0)
    cat(sprintf("  Trials failed:    %5d  <-- CHECK\n", oc$n_fail))
  cat(sprintf("  PCS (d%d):         %5.1f%%\n", sc$true_opt, 100*oc$PCS))
  cat(sprintf("  Pr(inadmissible): %5.1f%%\n",             100*oc$P_inadm))
  cat(sprintf("  Pr(early stop):   %5.1f%%\n",             100*oc$P_stop))
  cat(sprintf("  Pr(no rec):       %5.1f%%\n",             100*oc$P_norec))
  cat(sprintf("  Mean N:           %5.1f\n",               oc$mean_N))
  cat(sprintf("  Mean tox rate:    %5.3f\n",               oc$mean_tox))
  cat(sprintf("  Mean eff rate:    %5.3f\n",               oc$mean_eff))
  cat(sprintf("  Dose selection:   %s\n",
      paste(sprintf("d%d=%.0f%%", 1:J, oc$sel_pct), collapse="  ")))
  cat(sprintf("  Mean N per dose:  %s\n\n",
      paste(sprintf("d%d=%.1f", 1:J, oc$n_dose), collapse="  ")))
}

# ============================================================
# SAVE OUTPUTS
# ============================================================

# 1. Raw results
saveRDS(all_results, "sims_100_results.rds")
cat("Saved: sims_100_results.rds\n")

# 2. Summary CSV
summary_rows <- lapply(names(scenarios), function(sc_name) {
  sc <- scenarios[[sc_name]]
  oc <- oc_list[[sc_name]]
  data.frame(
    Scenario     = sc_name,
    Name         = sc$name,
    N            = sc$N,
    TrueOpt      = sc$true_opt,
    R            = oc$R,
    N_fail       = oc$n_fail,
    PCS          = round(100*oc$PCS,    1),
    P_inadm      = round(100*oc$P_inadm,1),
    P_stop       = round(100*oc$P_stop, 1),
    P_norec      = round(100*oc$P_norec,1),
    MeanN        = round(oc$mean_N,     1),
    MeanTox      = round(oc$mean_tox,   3),
    MeanEff      = round(oc$mean_eff,   3),
    SelD1        = oc$sel_pct[1],
    SelD2        = oc$sel_pct[2],
    SelD3        = oc$sel_pct[3],
    SelD4        = oc$sel_pct[4],
    SelD5        = oc$sel_pct[5],
    NdoseD1      = oc$n_dose[1],
    NdoseD2      = oc$n_dose[2],
    NdoseD3      = oc$n_dose[3],
    NdoseD4      = oc$n_dose[4],
    NdoseD5      = oc$n_dose[5],
    Runtime_sec  = round(timing[sc_name], 1),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, "sims_100_summary.csv", row.names=FALSE)
cat("Saved: sims_100_summary.csv\n")

# 3. Formatted text table
sink("sims_100_summary.txt")
cat(sprintf("BART-EffTox simulation: R=%d trials per scenario\n", R_pilot))
cat(sprintf("Seed=%d  M=%d  n_iter=%d  n_burn=%d\n",
    MASTER_SEED, DESIGN$M, DESIGN$n_iter, DESIGN$n_burn))
cat(sprintf("phi_T=%.2f  phi_E=%.2f  c_T_pre=%.2f  c_T=%.2f  c_E=%.2f\n\n",
    DESIGN$phi_T, DESIGN$phi_E, DESIGN$c_T_pre, DESIGN$c_T, DESIGN$c_E))
cat(sprintf("%-4s %-32s %4s %5s %6s %6s %6s %6s\n",
    "Scen", "Name", "N", "PCS%", "Inadm%", "Stop%", "NoRec%", "MeanN"))
cat(strrep("-", 75), "\n")
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  oc <- oc_list[[sc_name]]
  cat(sprintf("%-4s %-32s %4d %5.1f %6.1f %6.1f %6.1f %6.1f\n",
      sc_name, substr(sc$name, 1, 32), sc$N,
      100*oc$PCS, 100*oc$P_inadm, 100*oc$P_stop, 100*oc$P_norec, oc$mean_N))
}
sink()
cat("Saved: sims_100_summary.txt\n")

cat("\nDone.\n")
