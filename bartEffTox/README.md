# bartEffTox

**Bayesian Nonparametric EffTox Design for Phase I/II Oncology Dose Finding**

[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://cran.r-project.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This package implements the **BART-EffTox** design, a Bayesian Phase I/II
dose-finding design that extends the EffTox decision framework of
[Thall and Cook (2004)](https://doi.org/10.1111/j.0006-341X.2004.00218.x)
by replacing its parametric dose-response model with separate probit
Bayesian Additive Regression Trees (BART) models for efficacy and
toxicity. The two posterior outcome models are combined at the decision
stage through a common efficacy-toxicity utility function.

A pool-adjacent-violators (PAVA) isotonic projection enforces monotone
toxicity at each decision stage without constraining the efficacy surface,
allowing the design to accommodate non-monotone efficacy profiles and
patient-level covariate effects.

## Installation

```r
# Install from GitHub (requires the remotes package)
# install.packages("remotes")
remotes::install_github("[username]/bartEffTox")

# The BART package is required for the BART fitting step
install.packages("BART")
```

## Quick start

```r
library(bartEffTox)

# Define a simple scenario (Scenario S1 from the manuscript)
sc <- list(
  name     = "Monotone upper-dose optimum",
  pT       = c(0.05, 0.10, 0.18, 0.25, 0.45),
  pE       = c(0.10, 0.25, 0.40, 0.55, 0.65),
  covars   = FALSE,
  N        = 36L,
  true_opt = 4L,
  true_adm = c(FALSE, TRUE, TRUE, TRUE, FALSE)
)

# Simulate R = 100 trials (set n_cores to your machine's core count)
set.seed(2025)
results <- simulate_design(
  scenario   = sc,
  R          = 100L,
  n_cores    = 2L,
  max_n      = 36L,
  cohort     = 3L,
  phi_T      = 0.30, phi_E = 0.20,
  c_T_pre    = 0.70, c_T   = 0.40, c_E = 0.50,
  M          = 20L,  n_iter = 400L, n_burn = 80L
)

# Summarize operating characteristics
oc <- summarise_ocs(results, true_opt = 4L,
                    true_adm = c(FALSE, TRUE, TRUE, TRUE, FALSE))
cat(sprintf("PCS = %.1f%%  P_inadm = %.1f%%  P_none = %.1f%%\n",
    100 * oc$PCS, 100 * oc$P_inadm, 100 * oc$P_none))
```

## Functions

| Function | Description |
|---|---|
| `utility_fn()` | Bilinear efficacy-toxicity utility (eq. 1) |
| `pava_project()` | Weighted PAVA isotonic projection (Proposition 3) |
| `runin_admissibility()` | Beta-binomial run-in safety rule (Proposition 4) |
| `bart_efftox_interim()` | Single BART-phase interim analysis (Section 3.3-3.4) |
| `run_trial()` | Simulate one complete two-phase trial (Section 3.5) |
| `simulate_design()` | Parallel simulation over R trials |
| `summarise_ocs()` | Operating characteristic summary (Tables 2-3) |

## Datasets

- `bart_efftox_scenarios`: Seven simulation scenarios from the manuscript
  (Table 1, Section 5.1), including non-monotone efficacy, biomarker-modified
  efficacy, frailty-modified toxicity, and combined covariate scenarios.

## Reproducing manuscript results

```r
# Load the seven manuscript scenarios
data(bart_efftox_scenarios)

# Run the full simulation (requires BART package; ~2h on 4 cores)
# See inst/scripts/run_simulations.R for the complete pipeline
```

See `inst/scripts/` for the full simulation pipeline and figure-generation
scripts that reproduce Tables 1-3 and Figures 1-6 of the manuscript.

## Design parameters

The default parameters match the manuscript calibration:

| Parameter | Value | Description |
|---|---|---|
| `phi_T` | 0.30 | Toxicity threshold |
| `phi_E` | 0.20 | Efficacy threshold |
| `c_T_pre` | 0.70 | Run-in blocking cutoff (Proposition 4) |
| `c_T` | 0.40 | BART-phase toxicity admissibility cutoff |
| `c_E` | 0.50 | BART-phase efficacy admissibility cutoff |
| `M` | 20 | BART trees |
| `n_iter` | 400 | Posterior draws (1000 for final analysis) |
| `n_burn` | 80 | Burn-in draws |

## Requirements

- R >= 4.1.0
- [BART](https://cran.r-project.org/package=BART) >= 2.9
  (Sparapani, Spanbauer & McCulloch 2021)
- parallel (base R)

**macOS note:** Use `n_cores <= 4` and `chunk_size = 50` to avoid
silent worker crashes with `mclapply` on Apple Silicon.

## Citation

If you use this package, please cite:

> [Author names redacted for review] (2025). A Bayesian Nonparametric
> EffTox Design for Phase I/II Oncology Dose Finding with Monotone
> Toxicity Control and Covariate Adaptation. *Statistics in Medicine*,
> [under review].

and the underlying BART implementation:

> Sparapani R, Spanbauer C, McCulloch R (2021). Nonparametric machine
> learning and efficient computation with Bayesian additive regression
> trees: the BART R package. *Journal of Statistical Software*
> **97**(1):1-66. https://doi.org/10.18637/jss.v097.i01

## References

- Thall PF, Cook JD (2004). Dose-finding based on efficacy-toxicity
  trade-offs. *Biometrics* **60**(3):684-693.
- Chipman HA, George EI, McCulloch RE (2010). BART: Bayesian additive
  regression trees. *Annals of Applied Statistics* **4**(1):266-298.
- Robertson T, Wright FT, Dykstra RL (1988). *Order Restricted
  Statistical Inference*. Wiley, New York.
- Albert JH, Chib S (1993). Bayesian analysis of binary and
  polychotomous response data. *JASA* **88**(422):669-679.
