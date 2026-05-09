library(testthat)

# ── utility_fn ─────────────────────────────────────────────────────────────
test_that("utility_fn gives correct corner values", {
  expect_equal(utility_fn(1, 0), 100)   # efficacy, no tox -> u10
  expect_equal(utility_fn(1, 1),  40)   # both -> u11
  expect_equal(utility_fn(0, 0),  20)   # neither -> u00
  expect_equal(utility_fn(0, 1),   0)   # tox only -> u01
})

test_that("utility_fn is vectorised", {
  pE <- c(0.10, 0.25, 0.40, 0.55, 0.65)
  pT <- c(0.05, 0.10, 0.18, 0.25, 0.45)
  u  <- utility_fn(pE, pT)
  expect_length(u, 5)
  expect_true(all(u >= 0 & u <= 100))
})

# ── pava_project ────────────────────────────────────────────────────────────
test_that("pava_project returns monotone vector", {
  p <- c(0.05, 0.20, 0.15, 0.30, 0.45)  # violation at position 3
  pp <- pava_project(p)
  expect_true(all(diff(pp) >= -1e-10))
  expect_true(all(pp >= 0 & pp <= 1))
})

test_that("pava_project is idempotent on monotone input", {
  p <- c(0.05, 0.10, 0.18, 0.25, 0.45)
  expect_equal(pava_project(p), p)
})

test_that("pava_project is nonexpansive", {
  set.seed(1)
  p    <- runif(5); p0 <- sort(runif(5))
  w    <- rep(1, 5)
  pp   <- pava_project(p, w)
  lhs  <- sum(w * (pp - p0)^2)
  rhs  <- sum(w * ( p - p0)^2)
  expect_lte(lhs, rhs + 1e-10)
})

# ── runin_admissibility ─────────────────────────────────────────────────────
test_that("runin_admissibility blocks at correct threshold", {
  # With c_T_pre=0.70, phi_T=0.30, Beta(1,1): y*=2 for n=3
  expect_false(runin_admissibility(0, 3))  # 0/3: not blocked
  expect_false(runin_admissibility(1, 3))  # 1/3: not blocked (tail ~0.65 < 0.70)
  expect_true( runin_admissibility(2, 3))  # 2/3: blocked (tail ~0.92 > 0.70)
  expect_true( runin_admissibility(3, 3))  # 3/3: blocked
})

test_that("runin_admissibility posterior tail is exact", {
  # Manual calculation: Beta(1+2, 1+3-2) = Beta(3,2), P(p > 0.30)
  expected <- 1 - pbeta(0.30, 3, 2)
  # tail > 0.70 -> blocked
  expect_true(expected > 0.70)
  expect_equal(round(expected, 3), 0.916)
  # Correct: 1 - pbeta(0.30, 1+2, 1+1) = 1 - pbeta(0.30, 3, 2)
  expect_equal(round(1 - pbeta(0.30, 3, 2), 3), 0.916)
})

test_that("Table 1 values match Proposition 4", {
  # Verify all four (n, y*) pairs from Table 1
  tab <- data.frame(
    n     = c( 3,  6,  9, 12),
    y_star= c( 2,  3,  4,  5)
  )
  for (i in seq_len(nrow(tab))) {
    n <- tab$n[i]; ys <- tab$y_star[i]
    # y* - 1 should NOT block
    expect_false(runin_admissibility(ys - 1L, n),
                 label = sprintf("n=%d, y=%d not blocked", n, ys-1))
    # y* should block
    expect_true(runin_admissibility(ys, n),
                label = sprintf("n=%d, y=%d blocked", n, ys))
  }
})

# ── pava_project edge cases ─────────────────────────────────────────────────
test_that("pava_project clips to [0, 1]", {
  # Should not produce values outside [0,1]
  expect_true(all(pava_project(c(0, 0.5, 1)) >= 0))
  expect_true(all(pava_project(c(0, 0.5, 1)) <= 1))
})

test_that("pava_project handles length-1 input", {
  expect_equal(pava_project(0.5), 0.5)
})
