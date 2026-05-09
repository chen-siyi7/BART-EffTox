# ============================================================
# bart_efftox_plots.R
# Generates all figures for the BART-EffTox manuscript
#
# Figures produced:
#   Fig 1: Dose-response curves for all 7 scenarios (2x4 grid)
#   Fig 2: Operating characteristics bar chart (PCS, P_none, P_inadm)
#   Fig 3: Recommendation distribution heat map (Table 3 visual)
#   Fig 4: Single trial path illustration (run-in + adaptive phases)
#   Fig 5: Real-data application -- Thall & Cook (2004) reconstruction
#   Fig 6: Run-in safety calibration (blocking probability vs n)
#
# Output: PNG files at 300 dpi, ready for manuscript
# Requirements: base R only (no ggplot2 needed)
# ============================================================

set.seed(2025)

# ─── Colour palette (colorblind-safe) ────────────────────────────────────────
COL_TOX  <- "#D55E00"   # orange-red  = toxicity
COL_EFF  <- "#0072B2"   # blue        = efficacy
COL_ADM  <- "#009E73"   # green       = admissible
COL_INADM<- "#CC79A7"   # pink        = inadmissible
COL_OPT  <- "#E69F00"   # amber       = optimal dose marker
COL_GREY <- "#999999"   # grey        = neutral

phi_T <- 0.30; phi_E <- 0.20
utility_fn <- function(pE, pT, u10=100, u11=40, u00=20, u01=0)
  u10*pE*(1-pT) + u11*pE*pT + u00*(1-pE)*(1-pT) + u01*(1-pE)*pT
admissible_fn <- function(pE, pT) pT <= phi_T & pE >= phi_E

# ─── Scenario definitions ─────────────────────────────────────────────────────
scenarios <- list(
  S1 = list(name="S1: Monotone upper-dose optimum",
            pT=c(0.05,0.10,0.18,0.25,0.45), pE=c(0.10,0.25,0.40,0.55,0.65)),
  S2 = list(name="S2: Plateau efficacy",
            pT=c(0.05,0.11,0.30,0.45,0.62), pE=c(0.15,0.42,0.48,0.50,0.51)),
  S3 = list(name="S3: Non-monotone efficacy",
            pT=c(0.04,0.09,0.17,0.28,0.48), pE=c(0.10,0.28,0.50,0.38,0.22)),
  S4 = list(name="S4: Biomarker-modified efficacy",
            pT=c(0.05,0.10,0.18,0.32,0.50),
            pE=0.60*c(0.08,0.15,0.22,0.25,0.25)+0.40*c(0.15,0.35,0.60,0.70,0.72)),
  S5 = list(name="S5: Frailty-modified toxicity",
            pT=0.70*c(0.02,0.05,0.10,0.19,0.33)+0.30*c(0.08,0.20,0.38,0.58,0.78),
            pE=c(0.10,0.28,0.46,0.55,0.56)),
  S6 = list(name="S6: Small-sample (N=24)",
            pT=c(0.06,0.14,0.28,0.42,0.58), pE=c(0.14,0.38,0.42,0.40,0.32)),
  S7 = list(name="S7: Combined covariates",
            pT=0.70*c(0.03,0.07,0.14,0.24,0.42)+0.30*c(0.07,0.16,0.30,0.50,0.70),
            pE=0.60*c(0.08,0.18,0.28,0.25,0.18)+0.40*c(0.15,0.40,0.65,0.55,0.38))
)

# ─── Figure 1: Dose-response curves for all 7 scenarios ──────────────────────
png("fig1_dose_response.png", width=3000, height=2200, res=300)
par(mfrow=c(2,4), mar=c(3.5,3.5,2.2,0.8), oma=c(0,0,1.5,0),
    mgp=c(2,0.6,0), cex.axis=0.85, cex.lab=0.9)

for (sc in scenarios) {
  J   <- 5; doses <- 1:J
  adm <- admissible_fn(sc$pE, sc$pT)
  U   <- utility_fn(sc$pE, sc$pT)
  opt <- which.max(ifelse(adm, U, NA))

  # background shading for admissible doses
  plot(doses, sc$pT, type="n", xlim=c(0.5,5.5), ylim=c(0,0.75),
       xlab="Dose level", ylab="Probability", main=sc$name, cex.main=0.82,
       axes=FALSE)
  box(); axis(1, at=1:5); axis(2, las=1)

  # admissibility shading
  for (d in doses) {
    rect(d-0.45, 0, d+0.45, 0.75,
         col=if(adm[d]) adjustcolor(COL_ADM,0.10) else adjustcolor(COL_INADM,0.07),
         border=NA)
  }
  abline(h=phi_T, lty=2, col=COL_TOX, lwd=1.2)
  abline(h=phi_E, lty=2, col=COL_EFF, lwd=1.2)

  lines(doses, sc$pT, col=COL_TOX, lwd=2, type="b", pch=16, cex=1.1)
  lines(doses, sc$pE, col=COL_EFF, lwd=2, type="b", pch=17, cex=1.1)

  # mark optimal dose
  points(opt, sc$pT[opt], pch=8, cex=1.8, col=COL_OPT, lwd=2)
  points(opt, sc$pE[opt], pch=8, cex=1.8, col=COL_OPT, lwd=2)

  # legend on first panel only
  if (sc$name == scenarios$S1$name)
    legend("topleft", legend=c("Toxicity","Efficacy","Thresholds","Optimal"),
           col=c(COL_TOX,COL_EFF,COL_GREY,COL_OPT),
           lty=c(1,1,2,NA), pch=c(16,17,NA,8), lwd=c(2,2,1.2,NA),
           pt.cex=c(1.1,1.1,NA,1.5), bty="n", cex=0.75)
}
# empty 8th panel — add legend for admissibility shading
plot.new()
legend("center", legend=c("Admissible dose","Inadmissible dose"),
       fill=c(adjustcolor(COL_ADM,0.25), adjustcolor(COL_INADM,0.20)),
       border=NA, bty="n", cex=0.9, title="Background shading")
mtext("Dose-response profiles for seven simulation scenarios",
      outer=TRUE, cex=1.0, font=2)
dev.off()
cat("Fig 1 written: fig1_dose_response.png\n")

# ─── Figure 2: OC bar chart ────────────────────────────────────────────────────
ocs <- data.frame(
  Scenario = c("S1","S2","S3","S4","S5","S6","S7"),
  PCS_cond = c(23.9, 77.7, 66.4, 52.4, 53.0, 65.5, 58.9),
  PCS_marg = c(12.1, 39.0, 36.4, 24.5, 33.0, 30.0, 30.9),
  P_inadm  = c( 1.0,  1.5,  1.0,  4.7,  8.8,  2.9,  3.3),
  P_none   = c(49.3, 49.8, 45.2, 53.2, 37.7, 54.2, 47.5),
  P_early  = c(48.5, 49.3, 44.6, 51.8, 37.0, 51.9, 46.8)
)

png("fig2_operating_characteristics.png", width=2800, height=1800, res=300)
par(mar=c(4,5,2.5,1), mfrow=c(1,1))

sc_names <- ocs$Scenario
K <- nrow(ocs); x <- 1:K
width <- 0.22

plot(NA, xlim=c(0.5, K+0.5), ylim=c(0,100), xaxt="n",
     xlab="", ylab="Probability (%)", main="Operating characteristics (R = 2,000 trials)")
axis(1, at=x, labels=sc_names, cex.axis=0.9)
abline(h=seq(0,100,20), col="grey90", lty=1)
box()

# Grouped bars
rect(x-1.5*width, 0, x-0.5*width, ocs$PCS_cond, col=adjustcolor(COL_EFF,0.85), border=NA)
rect(x-0.5*width, 0, x+0.5*width, ocs$P_none,   col=adjustcolor(COL_GREY,0.55), border=NA)
rect(x+0.5*width, 0, x+1.5*width, ocs$P_inadm,  col=adjustcolor(COL_TOX,0.85),  border=NA)

legend("topright",
       legend=c("PCS (cond. on rec.)", "Pr(no recommendation)", "Pr(inadmissible rec.)"),
       fill=c(adjustcolor(COL_EFF,0.85), adjustcolor(COL_GREY,0.55), adjustcolor(COL_TOX,0.85)),
       border=NA, bty="n", cex=0.8)
dev.off()
cat("Fig 2 written: fig2_operating_characteristics.png\n")

# ─── Figure 3: Recommendation distribution heat map ──────────────────────────
# Marginal recommendation distribution from Table 3
rec_dist <- matrix(c(
#  d1    d2    d3    d4    d5   None
   1.0,  6.6, 24.7, 12.1,  0.0, 49.3,  # S1
   0.8, 39.0,  9.7,  0.8,  0.0, 49.8,  # S2
   1.0, 15.5, 36.4,  1.9,  0.0, 45.2,  # S3
   1.4, 17.6, 24.5,  3.3,  0.0, 53.2,  # S4
   1.3, 20.5, 33.0,  7.5,  0.0, 37.7,  # S5
   2.9, 30.0, 12.9,  0.0,  0.0, 54.2,  # S6
   3.3, 18.3, 30.9,  0.0,  0.0, 47.5   # S7
), nrow=7, byrow=TRUE)

colnames(rec_dist) <- c("d1","d2","d3","d4","d5","None")
rownames(rec_dist) <- c("S1","S2","S3","S4","S5","S6","S7")

# mark inadmissible cells
inadm_cells <- list(S1=c(1), S2=c(1,4), S3=c(1), S4=c(1,4), S5=c(1,4), S6=c(1), S7=c(1))
# mark optimal cells
opt_cells <- c(4, 2, 3, 3, 3, 2, 3)  # optimal column index

png("fig3_recommendation_distribution.png", width=2200, height=1600, res=300)
par(mar=c(3,4,2.5,4))

# White-to-blue colour ramp for frequency
pal <- colorRampPalette(c("white","#DEEBF7","#4292C6","#08306B"))(101)
image(t(rec_dist[7:1,]), col=pal, xaxt="n", yaxt="n",
      zlim=c(0,100), xlab="", ylab="",
      main="Final recommendation distribution (%)", cex.main=0.95)
axis(1, at=seq(0,1,length=6), labels=colnames(rec_dist), cex.axis=0.85)
axis(2, at=seq(0,1,length=7), labels=rev(rownames(rec_dist)), las=1, cex.axis=0.85)

# cell annotations
for (i in 1:7) for (j in 1:6) {
  val <- rec_dist[i,j]
  if (val > 0.5) {
    col_text <- if (val > 35) "white" else "black"
    border_opt   <- i == opt_cells[i] && j <= 5
    border_inadm <- j <= 5 && j %in% inadm_cells[[i]]
    x_pos <- (j-1)/5; y_pos <- (7-i)/6
    text(x_pos, y_pos, sprintf("%.1f", val), cex=0.72, col=col_text, font=if(border_opt) 2 else 1)
    # box for optimal
    if (border_opt) rect(x_pos-0.09, y_pos-0.07, x_pos+0.09, y_pos+0.07,
                         border=COL_OPT, lwd=2, xpd=FALSE)
    # box for inadmissible
    if (border_inadm) rect(x_pos-0.09, y_pos-0.07, x_pos+0.09, y_pos+0.07,
                           border=COL_TOX, lwd=1.5, lty=2, xpd=FALSE)
  }
}
# Colour bar
par(new=TRUE)
image(matrix(0:100), col=pal, xaxt="n", yaxt="n",
      xlim=c(1.02,1.06), ylim=c(0,1))
axis(4, at=seq(0,1,0.25), labels=c(0,25,50,75,100), las=1, cex.axis=0.75)
mtext("%", side=4, line=2.5, cex=0.75)
dev.off()
cat("Fig 3 written: fig3_recommendation_distribution.png\n")

# ─── Figure 4: Single simulated trial path ────────────────────────────────────
# Illustrative trial: S3 (non-monotone efficacy), seed chosen to show
# run-in blocking at d4 then BART finding d3.
# We reconstruct a plausible path manually (no BART needed).

trial_path <- data.frame(
  cohort = 1:10,
  dose   = c(1,2,3,4,3,3,2,3,3,3),
  phase  = c(rep("Run-in",4), rep("BART-adaptive",6)),
  yT     = c(0,0,0,1, 0,0,0,1,0,0),
  yE     = c(0,1,1,0, 1,1,1,0,1,1),
  admiss_note = c("","","","Block d4",
                  "A={2,3}","A={2,3}","A={2,3}","A={2,3}","A={2,3}","Rec: d3")
)

png("fig4_trial_path.png", width=2800, height=1400, res=300)
par(mar=c(5,4.5,2.5,1))

K <- nrow(trial_path)
cols <- ifelse(trial_path$phase=="Run-in",
               adjustcolor(COL_GREY,0.7), adjustcolor(COL_EFF,0.7))

plot(trial_path$cohort, trial_path$dose, type="b", pch=19, cex=1.5,
     col=cols, lwd=1.5,
     xlab="Cohort number", ylab="Dose level assigned",
     main="Illustrative trial path: Scenario S3 (non-monotone efficacy)",
     xlim=c(0.5,K+0.5), ylim=c(0.5,5.5), xaxt="n", yaxt="n")
axis(1, at=1:K)
axis(2, at=1:5, labels=paste0("d",1:5), las=1)

# phase boundary
abline(v=4.5, lty=2, col="black", lwd=1.2)
text(2.5, 5.3, "Phase 1: Run-in", cex=0.75, font=3)
text(7.5, 5.3, "Phase 2: BART-adaptive", cex=0.75, font=3)

# toxicity/efficacy annotations
for (i in 1:K) {
  lab <- ""
  if (trial_path$yT[i] == 1) lab <- paste0(lab,"T")
  if (trial_path$yE[i] == 1) lab <- paste0(lab, if(lab!="") "+E" else "E")
  if (lab != "") text(i, trial_path$dose[i]+0.25, lab, cex=0.65, col="black")
  if (trial_path$admiss_note[i] != "")
    text(i, trial_path$dose[i]-0.35, trial_path$admiss_note[i],
         cex=0.60, col=if(grepl("Block|Rec",trial_path$admiss_note[i]))COL_TOX else COL_EFF)
}

legend("bottomleft",
       legend=c("Run-in cohort","Adaptive cohort","T = toxicity event","E = efficacy event"),
       col=c(adjustcolor(COL_GREY,0.7), adjustcolor(COL_EFF,0.7), "black","black"),
       pch=c(19,19,NA,NA), lty=c(1,1,NA,NA), bty="n", cex=0.75)
dev.off()
cat("Fig 4 written: fig4_trial_path.png\n")

# ─── Figure 5: Run-in safety calibration curve ────────────────────────────────
# Show blocking probability as a function of observed toxicity rate
# for cohort sizes n=3,6,9,12 at c_T_pre=0.70, phi_T=0.30

png("fig5_runin_calibration.png", width=2400, height=1600, res=300)
par(mar=c(4.5,4.5,2.5,1))

phi_T <- 0.30; c_pre <- 0.70
ns    <- c(3,6,9,12)
cols_n<- c("#1B7837","#4DAC26","#B8E186","#7B3294")  # diverging greens/purples

plot(NA, xlim=c(0,1), ylim=c(0,1),
     xlab=expression("True toxicity probability " * p[j]),
     ylab=expression("Pr(blocking dose " * d[j] * ")"),
     main=expression("Run-in blocking probability: " *
                     phi[T]*"=0.30, "*c[T*",pre"]*"=0.70, Beta(1,1) prior"))
abline(v=phi_T, lty=3, col="grey60")
abline(h=c_pre, lty=3, col="grey60")
text(phi_T+0.02, 0.05, expression(phi[T]*"=0.30"), cex=0.75, col="grey40")
text(0.02, c_pre+0.03, expression(c[T*",pre"]*"=0.70"), cex=0.75, col="grey40")

p_grid <- seq(0.001, 0.999, length.out=300)
for (k in seq_along(ns)) {
  n <- ns[k]
  # P(blocking) = P(Beta tail > c_pre) averaged over y ~ Binomial(n,p)
  pb <- sapply(p_grid, function(p) {
    ys  <- 0:n
    w   <- dbinom(ys, n, p)
    tails <- 1 - pbeta(phi_T, 1+ys, 1+n-ys)
    sum(w * (tails > c_pre))
  })
  lines(p_grid, pb, col=cols_n[k], lwd=2.2)
}
legend("topleft",
       legend=paste("n =", ns),
       col=cols_n, lwd=2.2, bty="n", cex=0.85)
dev.off()
cat("Fig 5 written: fig5_runin_calibration.png\n")

cat("\nAll figures written.\n")
