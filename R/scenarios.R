#' Seven simulation scenarios for BART-EffTox calibration
#'
#' A list of seven scenario specifications matching Table~1 and
#' Section~5.1 of the manuscript. Each scenario defines the true
#' dose-response surface for \eqn{J = 5} dose levels, along with the
#' true optimal dose and admissible set.
#'
#' @format A named list with elements \code{S1} through \code{S7}.
#' Each element is a list with fields:
#' \describe{
#'   \item{\code{name}}{Character. Scenario description.}
#'   \item{\code{N}}{Integer. Maximum sample size.}
#'   \item{\code{pT}}{Numeric vector length 5. Marginal toxicity probabilities.}
#'   \item{\code{pE}}{Numeric vector length 5. Marginal efficacy probabilities.}
#'   \item{\code{true_opt}}{Integer. True optimal dose index.}
#'   \item{\code{true_adm}}{Logical vector. Truly admissible doses.}
#'   \item{\code{covars}}{Logical. Whether patient covariates are simulated.}
#'   \item{Covariate fields}{Optional: \code{pB}, \code{pF},
#'     \code{pT_xb0}, \code{pT_xb1}, \code{pE_xb0}, \code{pE_xb1},
#'     \code{pT_xf0}, \code{pT_xf1}, \code{pE_xf0}, \code{pE_xf1}.}
#' }
#'
#' @references
#' Manuscript Table~1, Section~5.1.
#'
#' @examples
#' data(bart_efftox_scenarios)
#' bart_efftox_scenarios$S1$pT
#' bart_efftox_scenarios$S2$true_opt
"bart_efftox_scenarios"
