###############################################################################
# Paper 3 - spatial_M4c.R   |   seed 20260826
#   JOB 1  SDEM theta / indirect-effect robustness across WELL-SPECIFIED graphs
#          only (fully connected, no isolates, non-pathological density).
#          Ill-posed graphs reported separately. queen_aug stays PRIMARY -
#          stability check, NOT a re-selection.
#   JOB 2  rebuild the stale M4 figures from the M4b (SDEM) result:
#          fig_M4b_impacts.png, fig_M4b_model_comparison.png
#          (fig_M4_*.png are left in place for before/after).
#
# Reads java_merged.gpkg, M4b_model_comparison.csv, M4b_impacts_with_CI.csv.
# Outputs (code/): M4c_theta_robustness.csv, M4c_illposed_graphs.csv,
#   fig_M4b_impacts.png, fig_M4b_model_comparison.png, M4c_results.txt
###############################################################################

suppressPackageStartupMessages({ library(sf); library(spdep); library(spatialreg) })
sf::sf_use_s2(FALSE)
SEED <- 20260826L
set.seed(SEED)
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else normalizePath(getwd())
})
setwd(here)
if (!all(file.exists("java_merged.gpkg", "M4b_impacts_with_CI.csv")))
  stop("run spatial_M1_M3.R and spatial_M4b.R first")

OUT_LOG <- file.path(here, "M4c_results.txt")
if (file.exists(OUT_LOG)) invisible(file.remove(OUT_LOG))
log_con <- file(OUT_LOG, open = "a")
log <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, log_con); flush(log_con) }
logcap <- function(x) { o <- capture.output(print(x)); cat(o, sep = "\n"); cat("\n")
  writeLines(o, log_con); flush(log_con) }
rule <- function(t = "") log("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78))

rule("PAPER 3  -  M4c  (SDEM theta robustness; M4b figures)")
log("run ", as.character(Sys.time()), " | seed ", SEED,
    " | chosen model (M4b) = SDEM | queen_aug PRIMARY")

## ---- data + graph builders ------------------------------------------
java <- st_read(file.path(here, "java_merged.gpkg"), quiet = TRUE)
java <- java[order(java$bps), ]
jdf  <- st_drop_geometry(java); jdf$logY <- base::log(jdf$Y)
FRM  <- logY ~ X1 + X2 + X3
aeqd <- "+proj=aeqd +lat_0=-7.5 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
jp <- st_transform(java, aeqd)
coords <- st_coordinates(st_point_on_surface(st_geometry(jp)))
add_sym_link <- function(nb, i, j) {
  put <- function(v, x) sort(unique(as.integer(c(v[v != 0L], x))))
  nb[[i]] <- put(nb[[i]], j); nb[[j]] <- put(nb[[j]], i); nb
}
ix  <- function(b) which(java$bps == b)
AUG <- rbind(c(3526L, 3578L), c(3101L, 3175L))
augment <- function(nb) {
  for (r in seq_len(nrow(AUG))) nb <- add_sym_link(nb, ix(AUG[r, 1]), ix(AUG[r, 2]))
  sym.attr.nb(nb)
}

## snap-tolerance probe: <=100 m no change; 500 m adds 4 near-touch links but
## Madura/Kep.Seribu stay disconnected; 5000 m over-connects across sea gaps.
log("Queen snap probe (m -> links / components):")
for (s in c(0, 100, 500, 1000, 5000)) {
  qq <- suppressWarnings(poly2nb(jp, queen = TRUE, snap = s))
  log("  snap ", sprintf("%5d", s), " -> ", sum(card(qq)), " / ", n.comp.nb(qq)$nc)
}

q119   <- poly2nb(jp, queen = TRUE, row.names = as.character(java$bps))
qa119  <- augment(q119)
r119   <- poly2nb(jp, queen = FALSE, row.names = as.character(java$bps))
ra119  <- augment(r119)
qs119  <- suppressWarnings(poly2nb(jp, queen = TRUE, snap = 500,
                                   row.names = as.character(java$bps)))
qsa119 <- augment(qs119)

keep   <- java$bps != 3101L                # drop Kepulauan Seribu isolate
jdf118 <- jdf[keep, ]
jp118  <- jp[keep, ]
q118   <- poly2nb(jp118, queen = TRUE, row.names = as.character(java$bps[keep]))
qa118  <- add_sym_link(q118, which(java$bps[keep] == 3526L),
                             which(java$bps[keep] == 3578L))   # Suramadu only (3101 gone)
qa118  <- sym.attr.nb(qa118)

klw <- function(k, cc = coords, rn = as.character(java$bps))
  nb2listw(knn2nb(knearneigh(cc, k = k), row.names = rn, sym = FALSE),
           style = "W", zero.policy = TRUE)

## inverse-distance band (ill-posed - carried only for the excluded table)
nn1  <- knn2nb(knearneigh(coords, k = 1)); dmax <- max(unlist(nbdists(nn1, coords)))
db   <- dnearneigh(coords, 0, dmax, row.names = as.character(java$bps))
while (n.comp.nb(db)$nc > 1) { dmax <- dmax * 1.15
  db <- dnearneigh(coords, 0, dmax, row.names = as.character(java$bps)) }
idw_lw <- nb2listw(db, glist = lapply(nbdists(db, coords), function(x) 1 / (x / 1000)),
                   style = "W", zero.policy = TRUE)

meta <- function(nb_or_lw) {
  nb <- if (inherits(nb_or_lw, "listw")) nb_or_lw$neighbours else nb_or_lw
  c(n = length(nb), ncomp = n.comp.nb(nb)$nc, isolates = sum(card(nb) == 0),
    mean_links = round(mean(card(nb)), 2))
}

## note/reason strings kept verbatim from the original run (they are CSV columns)
WELL <- list(
  queen_aug         = list(lw = nb2listw(qa119,  style = "W", zero.policy = TRUE), data = jdf,    note = "PRIMARY"),
  queen_aug_n118    = list(lw = nb2listw(qa118,  style = "W", zero.policy = TRUE), data = jdf118, note = "primary minus Kep. Seribu, still fully connected"),
  queen_bare_n118   = list(lw = nb2listw(q118,   style = "W", zero.policy = TRUE), data = jdf118, note = "isolate dropped but NOT augmented: Madura still a separate component (2 subgraphs)"),
  rook_aug          = list(lw = nb2listw(ra119,  style = "W", zero.policy = TRUE), data = jdf,    note = "Rook contiguity + 2 documented links"),
  queen_snap500_aug = list(lw = nb2listw(qsa119, style = "W", zero.policy = TRUE), data = jdf,    note = "Queen(snap 500 m) + 2 documented links; +4 near-touch links vs queen_aug"),
  knn4              = list(lw = klw(4), data = jdf, note = "k-nearest, k=4"),
  knn5              = list(lw = klw(5), data = jdf, note = "k-nearest, k=5"),
  knn6              = list(lw = klw(6), data = jdf, note = "k-nearest, k=6")
)
ILL <- list(
  queen_bare_n119 = list(lw = nb2listw(q119, style = "W", zero.policy = TRUE), data = jdf,
    reason = "3 disconnected subgraphs: Kepulauan Seribu isolated + Madura a closed 4-unit block; a spatial model on a disconnected graph is not estimating the same system"),
  idw_band_0_141km = list(lw = idw_lw, data = jdf,
    reason = paste0("inverse-distance band 0-", round(dmax / 1000, 0),
    " km, mean ", round(mean(card(db)), 0),
    " neighbours: over-dense, known wide-band pathology (X1 indirect -> -21.9 in M4b)"))
)

## refit SDEM on a graph; theta (= INDIRECT for SDEM) + normal-approx 95% CI
sdem_row <- function(tag, lw, data, note, R = 2000) {
  m  <- errorsarlm(FRM, data = data, listw = lw, Durbin = TRUE, zero.policy = TRUE)
  Co <- summary(m)$Coef
  lagr <- grep("^lag\\.", rownames(Co))
  tr <- trW(as(as_dgRMatrix_listw(lw), "CsparseMatrix"), type = "mult")
  set.seed(SEED)
  s <- summary(impacts(m, tr = tr, R = R), zstats = TRUE, short = TRUE)
  ind_est <- s$mat[, "Indirect"]; ind_se <- s$semat[, "Indirect"]; ind_p <- s$pzmat[, "Indirect"]
  z <- qnorm(0.975); md <- meta(lw)
  do.call(rbind, lapply(1:3, function(k) {
    lr <- lagr[k]
    data.frame(graph = tag, note = note, n = md["n"], components = md["ncomp"],
      isolates = md["isolates"], mean_links = md["mean_links"], lambda = round(m$lambda, 4),
      variable = c("X1", "X2", "X3")[k],
      theta_est = round(Co[lr, 1], 5), theta_se = round(Co[lr, 2], 5),
      theta_z = round(Co[lr, 3], 3), theta_p = signif(Co[lr, 4], 4),
      indirect_pct = round(100 * ind_est[k], 4),
      indirect_ci_lo_pct = round(100 * (ind_est[k] - z * ind_se[k]), 4),
      indirect_ci_hi_pct = round(100 * (ind_est[k] + z * ind_se[k]), 4),
      indirect_p = signif(ind_p[k], 4), row.names = NULL)
  }))
}

## ---- JOB 1a: well-specified graphs -------------------------------
rule("JOB 1  -  SDEM theta / indirect robustness, WELL-SPECIFIED graphs")
inv <- do.call(rbind, lapply(names(WELL), function(g) {
  md <- meta(WELL[[g]]$lw)
  data.frame(graph = g, n = md["n"], components = md["ncomp"], isolates = md["isolates"],
             mean_links = md["mean_links"], note = WELL[[g]]$note, row.names = NULL)
}))
logcap(inv)
log("queen_bare_n118 has 2 components (Madura detached) - the literal 'isolate")
log("dropped' graph; the robustness CLAIM rests on the fully-connected graphs.")
well <- do.call(rbind, lapply(names(WELL), function(g)
  sdem_row(g, WELL[[g]]$lw, WELL[[g]]$data, WELL[[g]]$note)))
write.csv(well, file.path(here, "M4c_theta_robustness.csv"), row.names = FALSE)
for (v in c("X1", "X2", "X3")) {
  log("W.", v, ":")
  logcap(well[well$variable == v, c("graph", "theta_est", "theta_p", "indirect_pct",
              "indirect_ci_lo_pct", "indirect_ci_hi_pct", "indirect_p")])
}

## ---- JOB 1b: ill-posed graphs (excluded) ------------------------
rule("JOB 1  -  ILL-POSED graphs (EXCLUDED from the robustness claim)")
ill <- do.call(rbind, lapply(names(ILL), function(g) {
  r <- sdem_row(g, ILL[[g]]$lw, ILL[[g]]$data, ILL[[g]]$reason)
  r$reason_excluded <- ILL[[g]]$reason; r
}))
write.csv(ill, file.path(here, "M4c_illposed_graphs.csv"), row.names = FALSE)
for (g in names(ILL)) log("  ", g, " : ", ILL[[g]]$reason)
logcap(ill[, c("graph", "variable", "theta_est", "theta_p", "indirect_pct",
               "indirect_ci_lo_pct", "indirect_ci_hi_pct", "indirect_p")])

## ---- JOB 1: verdict -------------------------------------------
rule("JOB 1  -  VERDICT (a sign/magnitude, b p=0.05 crossing, c CI-width)")
fc   <- well[well$graph != "queen_bare_n118", ]        # fully-connected subset
cont <- c("queen_aug", "queen_aug_n118", "rook_aug", "queen_snap500_aug")
sm <- function(v, sub = fc) {
  d <- sub[sub$variable == v, ]
  ciw <- d$indirect_ci_hi_pct - d$indirect_ci_lo_pct
  list(allpos = all(d$theta_est > 0), allneg = all(d$theta_est < 0),
       imin = min(d$indirect_pct), imax = max(d$indirect_pct),
       pmin = min(d$indirect_p), pmax = max(d$indirect_p),
       ncross = sum(d$indirect_p < 0.05), which_cross = d$graph[d$indirect_p < 0.05],
       cwmin = min(ciw), cwmax = max(ciw), cwmean = mean(ciw), cwsd = sd(ciw))
}
for (v in c("X1", "X2", "X3")) {
  s <- sm(v)
  log("W.", v, " (7 fully-connected graphs): sign ",
      if (s$allpos) "+ all 7" else if (s$allneg) "- all 7" else "MIXED",
      " | indirect ", round(s$imin, 2), "% .. ", round(s$imax, 2), "%",
      " | p ", signif(s$pmin, 3), " .. ", signif(s$pmax, 3), " (<0.05 in ", s$ncross, "/7",
      if (s$ncross) paste0(": ", paste(s$which_cross, collapse = ",")) else "", ")",
      " | 95% CI width ", round(s$cwmin, 1), "-", round(s$cwmax, 1), " pp (sd ",
      round(s$cwsd, 2), ")")
}
sx2 <- sm("X2"); sx3 <- sm("X3"); sx3c <- sm("X3", fc[fc$graph %in% cont, ])
log("(a) W.X2 & W.X3 positive in all 7, indirect ~ +1 to +2.5% - stable to W. W.X1 negative, never near sig.")
log("(b) W.X2 crosses p=0.05 in ", sx2$ncross, "/7 (~0.045); W.X3 in ", sx3$ncross,
    "/7 (min ", signif(sx3$pmin, 3), "). Neither robustly significant.")
log("(c) CI width ~flat: W.X2 ", round(sx2$cwmin, 1), "-", round(sx2$cwmax, 1),
    " pp; W.X3 contiguity family ~", round(sx3c$cwmean, 1),
    " pp, widens only under kNN. => power ceiling at n=119, not a W-definition effect.")

## ---- JOB 2: figures -----------------------------------------
rule("JOB 2  -  rebuild M4 figures from the M4b (SDEM) result")
impCI <- read.csv(file.path(here, "M4b_impacts_with_CI.csv"), stringsAsFactors = FALSE)
mc6   <- read.csv(file.path(here, "M4b_model_comparison.csv"), stringsAsFactors = FALSE)
PAL_eff <- c(DIRECT = "#1f2a44", INDIRECT = "#c2410c", TOTAL = "#2166ac")

## fig_M4b_impacts.png : forest plot with 95% CIs
png(file.path(here, "fig_M4b_impacts.png"), width = 2000, height = 1250, res = 175)
par(mar = c(4.6, 10.5, 4.4, 5))
impCI$variable <- factor(impCI$variable, levels = c("X3", "X2", "X1"))
impCI$effect   <- factor(impCI$effect,   levels = c("TOTAL", "INDIRECT", "DIRECT"))
impCI <- impCI[order(impCI$variable, impCI$effect), ]
vlab <- c(X1 = "X1 unemp.", X2 = "X2 HDI", X3 = "X3 LFPR")
rowlab <- sprintf("%-9s %s", vlab[as.character(impCI$variable)], impCI$effect)
yy <- seq_len(nrow(impCI))
xr <- range(c(impCI$ci_lo_pct, impCI$ci_hi_pct, 0)) + c(-0.5, 0.5)
pfmt <- function(p) ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))
plot(NA, xlim = xr, ylim = c(0.5, nrow(impCI) + 0.5), yaxt = "n",
     xlab = "semi-elasticity: % change in poverty rate per one-unit change in predictor", ylab = "")
abline(v = 0, col = "grey55", lwd = 1.4)
abline(h = c(3.5, 6.5), col = "grey88")
segments(impCI$ci_lo_pct, yy, impCI$ci_hi_pct, yy, col = PAL_eff[as.character(impCI$effect)], lwd = 3)
points(impCI$pct_per_unit, yy, pch = 19, cex = 1.5, col = PAL_eff[as.character(impCI$effect)])
axis(2, at = yy, labels = rowlab, las = 1, cex.axis = 0.92, font = 1)
text(par("usr")[2], yy, pfmt(impCI$p_value), pos = 4, cex = 0.72, col = "grey30", xpd = NA)
title(main = "M4b  Impact decomposition - SDEM  (points = estimate, bars = 95% CI)",
      adj = 0, font.main = 1, cex.main = 1.2)
mtext(paste0("Under SDEM there is no global multiplier: INDIRECT = the theta (W.X) coefficient exactly.  ",
             "The W.X2 and W.X3 indirect CIs just include zero."),
      side = 3, line = 0.3, adj = 0, cex = 0.72, col = "grey30")
legend("bottomleft", bty = "n", cex = 0.9, lwd = 3, col = PAL_eff, legend = names(PAL_eff), horiz = TRUE)
dev.off(); log("wrote fig_M4b_impacts.png")

## fig_M4b_model_comparison.png : AIC + residual Moran's I, 6 models
ord <- c("OLS", "SLX", "SAR", "SEM", "SDM", "SDEM")
mc6 <- mc6[match(ord, mc6$model), ]
absorb <- mc6$resid_Moran_p >= 0.05
png(file.path(here, "fig_M4b_model_comparison.png"), width = 2000, height = 1150, res = 175)
par(mfrow = c(1, 2), oma = c(0, 0, 4, 0), mar = c(4, 4.5, 2, 1))
b1 <- barplot(mc6$AIC, names.arg = ord, border = NA, ylab = "AIC",
              col = ifelse(ord == "SDEM", "#b2182b", "#7f9cc0"), ylim = c(0, max(mc6$AIC) * 1.15))
text(b1, mc6$AIC, sprintf("%.1f", mc6$AIC), pos = 3, cex = 0.8, col = "grey20")
title(main = "(1)  AIC", adj = 0, font.main = 1)
mtext("AIC is a TIEBREAK only under our rule.  SDM 20.04 vs SDEM 20.54 = 0.5 apart.",
      side = 1, line = 2.6, adj = 0, cex = 0.68, col = "grey30")
b2 <- barplot(mc6$resid_MoranI, names.arg = ord, border = NA,
              ylab = "residual Moran's I  (queen_aug, 9999 perm)",
              col = ifelse(absorb, "#177245", "#c2410c"),
              ylim = range(c(mc6$resid_MoranI, 0)) + c(-0.05, 0.05))
abline(h = 0, col = "grey55")
text(b2, mc6$resid_MoranI,
     ifelse(absorb, sprintf("p=%.2f", mc6$resid_Moran_p),
            sprintf("p=%s", format(mc6$resid_Moran_p, scientific = TRUE, digits = 1))),
     pos = ifelse(mc6$resid_MoranI >= 0, 3, 1), cex = 0.72, col = "grey20")
title(main = "(2)  Residual spatial dependence", adj = 0, font.main = 1)
mtext("green = absorbed (p >= 0.05);  orange = NOT absorbed.  OLS, SLX, SAR fail.",
      side = 1, line = 2.6, adj = 0, cex = 0.68, col = "grey30")
mtext(paste0("M4b  Six-candidate comparison.  CHOSEN = SDEM (red) - not on AIC: ",
             "robust LMlag n.s. (0.62), short-range M3 correlogram, Hausman rejects SEM (p=0.002)."),
      outer = TRUE, cex = 0.92, font = 2, line = 1.2)
dev.off(); log("wrote fig_M4b_model_comparison.png")
log("fig_M4_impacts.png / fig_M4_model_comparison.png (M4/SDM) left in place for before/after")

rule("M4c COMPLETE - STOP")
close(log_con)
cat("\nDONE - see ", OUT_LOG, "\n")
