###############################################################################
# Paper 3 - spatial_M5_M6.R : M5 typology + M6 spillover-aware targeting
#   seed 20260826  |  reads java_merged.gpkg, M4b_impacts_with_CI.csv,
#                       M3_lisaY_all_schemes.csv
#
# M5  Continuous typology, all 119 units, under augmented Queen:
#     (a) Moran-scatterplot quadrant on Y (signs only)
#     (b) HDI differential X2 - W.X2 ; shadow needs W.X2 > mean(X2) (absolute)
#     -> 4 mechanism-named types, each >= 5.
# M6  Net regional benefit of +1 HDI in i, homogeneous SDEM impacts (no rho):
#       net_i = b_dir*Y_i  -  sum_j theta*Y_j*w_ji ,  w_ji = 1/k_j  (NOT w_ij)
#     R0 = rank by Y_i ; R1 = rank by net_i. Finding = the DIFFERENCE.
#     Sensitivity: theta = 0 and the two 95% CI bounds. RATES only (no headcount).
#
# Outputs (code/): M5_typology.csv, M5_crosstab.csv, fig_M5_typology.png,
#   M6_rankings.csv, M6_sensitivity.csv, M5_M6_results.txt
###############################################################################

suppressPackageStartupMessages({ library(sf); library(spdep) })
sf::sf_use_s2(FALSE)
SEED <- 20260826L
set.seed(SEED)
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else normalizePath(getwd())
})
setwd(here)
if (!all(file.exists("java_merged.gpkg", "M4b_impacts_with_CI.csv",
                     "M3_lisaY_all_schemes.csv")))
  stop("run spatial_M1_M3.R, spatial_M4b.R first")

OUT_LOG <- file.path(here, "M5_M6_results.txt")
if (file.exists(OUT_LOG)) invisible(file.remove(OUT_LOG))
log_con <- file(OUT_LOG, open = "a")
log <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, log_con); flush(log_con) }
logcap <- function(x) { o <- capture.output(print(x)); cat(o, sep = "\n"); cat("\n")
  writeLines(o, log_con); flush(log_con) }
rule <- function(t = "") log("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78))

rule("PAPER 3  -  M5 TYPOLOGY  +  M6 SPILLOVER-AWARE TARGETING")
log("run ", as.character(Sys.time()), " | seed ", SEED)

## ---- data + primary weights (augmented Queen) --------------------------
java <- st_read(file.path(here, "java_merged.gpkg"), quiet = TRUE)
java <- java[order(java$bps), ]
d <- st_drop_geometry(java)
aeqd <- "+proj=aeqd +lat_0=-7.5 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
jp <- st_transform(java, aeqd)
qnb <- poly2nb(jp, queen = TRUE, row.names = as.character(java$bps))
add_sym_link <- function(nb, i, j) {
  put <- function(v, x) sort(unique(as.integer(c(v[v != 0L], x))))
  nb[[i]] <- put(nb[[i]], j); nb[[j]] <- put(nb[[j]], i); nb
}
ix <- function(b) which(java$bps == b)
AUG <- rbind(c(3526L, 3578L), c(3101L, 3175L))
for (r in seq_len(nrow(AUG))) qnb <- add_sym_link(qnb, ix(AUG[r, 1]), ix(AUG[r, 2]))
qnb <- sym.attr.nb(qnb)
lw  <- nb2listw(qnb, style = "W", zero.policy = TRUE)
Wm  <- listw2mat(lw)                    # Wm[i,j] = w_ij (row i sums to 1)
kdeg <- card(qnb)
stopifnot(nrow(Wm) == 119, isTRUE(all.equal(unname(rowSums(Wm)), rep(1, 119))))

###############################################################################
# M5.  CONTINUOUS TYPOLOGY
###############################################################################
rule("M5  TYPOLOGY (continuous; all 119 units)")

## (a) Moran quadrant on Y - signs only
z <- d$Y - mean(d$Y)
lag_z <- as.vector(Wm %*% z)
quad <- ifelse(z > 0 & lag_z > 0, "HH", ifelse(z < 0 & lag_z < 0, "LL",
        ifelse(z > 0 & lag_z < 0, "HL", "LH")))
logcap(as.data.frame(table(quadrant = quad)))

## (b) HDI differential ; neighbourhood-mean HDI
lag_X2 <- as.vector(Wm %*% d$X2)
hdi_diff <- d$X2 - lag_X2
hdi_pos <- ifelse(hdi_diff >= 0, "own_HDI_above_nbrs", "in_higher_HDI_nbrhd")
log("HDI differential X2 - W.X2: range ", round(min(hdi_diff), 2), " .. ",
    round(max(hdi_diff), 2), " ; median ", round(median(hdi_diff), 2))
logcap(as.data.frame(table(hdi_position = hdi_pos)))
logcap(addmargins(table(quadrant = factor(quad, c("HH", "HL", "LH", "LL")),
  hdi_position = factor(hdi_pos, c("in_higher_HDI_nbrhd", "own_HDI_above_nbrs")))))

## Types:
##   poor_own  = z > 0 ; poor_nbrs = lag_z > 0
##   shadow    = hdi_diff < 0  AND  W.X2 > mean(X2)   (needs a genuinely
##               prosperous neighbour; the sign-only rule was too weak - it
##               caught e.g. Sampang, poorest regency, no prosperous neighbour)
##   Prosperous core  : !poor_own & !poor_nbrs         (LL)
##   Agglomeration shadow : poor_own & shadow          (SDEM-predicted type)
##   Lagging region   : poor_own & !shadow
##   Sheltered fringe : !poor_own & poor_nbrs          (LH)
X2_mean <- mean(d$X2)
poor_own  <- z > 0
poor_nbrs <- lag_z > 0
shadow    <- hdi_diff < 0 & lag_X2 > X2_mean
n_reclass <- sum(poor_own & hdi_diff < 0 & !(lag_X2 > X2_mean))
log("shadow requires W.X2 > mean(X2) = ", round(X2_mean, 2),
    " (median(X2) = ", round(median(d$X2), 2),
    " would be more permissive; pop-weighted mean not computable). ",
    n_reclass, " units reclassified shadow -> lagging vs the old sign-only rule.")
type0 <- ifelse(!poor_own & !poor_nbrs, "Prosperous core",
        ifelse( poor_own &  shadow,     "Agglomeration shadow",
        ifelse( poor_own & !shadow,     "Lagging region",
        ifelse(!poor_own &  poor_nbrs,  "Sheltered fringe", "Prosperous core"))))
log("initial type counts:")
logcap(as.data.frame(table(type = type0)))

MIN <- 5L
type <- type0
tc <- table(type)
small <- names(tc)[tc < MIN]
if (length(small)) for (s in small) {
  tgt <- if (s == "Sheltered fringe") "Prosperous core" else "Lagging region"
  log("MERGE '", s, "' (", tc[s], " < ", MIN, ") -> '", tgt, "'")
  type[type == s] <- tgt
} else log("no type has < ", MIN, " members; no merge")
log("final typology:")
logcap(as.data.frame(table(type = type)))

d$z <- z; d$lag_z <- lag_z; d$quad <- quad
d$lag_X2 <- lag_X2; d$hdi_diff <- hdi_diff; d$hdi_pos <- hdi_pos; d$type <- type
tsum <- do.call(rbind, lapply(sort(unique(type)), function(tp) {
  s <- d[d$type == tp, ]
  data.frame(type = tp, n = nrow(s), Y = round(mean(s$Y), 2), X1 = round(mean(s$X1), 2),
             X2 = round(mean(s$X2), 2), X3 = round(mean(s$X3), 2),
             WX2 = round(mean(s$lag_X2), 2), hdi_diff = round(mean(s$hdi_diff), 2),
             row.names = NULL)
}))
log("per-type means (WX2 = mean neighbourhood HDI; Java mean HDI ",
    round(X2_mean, 2), " - the shadow type's WX2 must be above this):")
logcap(tsum)

ct2 <- table(quadrant = factor(quad, c("HH", "HL", "LH", "LL")), type = type)
write.csv(cbind(quadrant = rownames(ct2), as.data.frame.matrix(ct2)),
          file.path(here, "M5_crosstab.csv"), row.names = FALSE)
logcap(addmargins(ct2))

d$WX2 <- lag_X2
d$nbrhd_prosperous <- lag_X2 > X2_mean
m5 <- d[, c("bps", "name", "prov", "Y", "X1", "X2", "X3", "z", "lag_z", "quad",
            "WX2", "hdi_diff", "hdi_pos", "nbrhd_prosperous", "type")]
m5[, c("z", "lag_z", "WX2", "hdi_diff")] <- round(m5[, c("z", "lag_z", "WX2", "hdi_diff")], 4)
write.csv(m5, file.path(here, "M5_typology.csv"), row.names = FALSE)
log("wrote M5_typology.csv, M5_crosstab.csv")

## sanity check: the 5 M3 FDR-significant LISA units
rule("M5  SANITY CHECK - the 5 FDR-significant LISA units")
lz <- read.csv(file.path(here, "M3_lisaY_all_schemes.csv"), stringsAsFactors = FALSE)
sig5 <- lz[lz$scheme == "QueenAug" & lz$variable == "Y" &
             (lz$sig_fdr == "TRUE" | lz$sig_fdr == TRUE), c("bps", "name", "quad")]
names(sig5)[3] <- "lisa_class"
sig5 <- merge(sig5, m5[, c("bps", "type", "quad", "hdi_diff")], by = "bps")
logcap(sig5[order(sig5$type), ])
coherent <- all(sig5$type[sig5$lisa_class == "Low-Low"] == "Prosperous core") &&
            all(sig5$type[sig5$lisa_class == "High-High"] %in%
                  c("Lagging region", "Agglomeration shadow"))
log("coherent (4 LL -> Prosperous core, HH -> poverty-side)? ", coherent)

## map
rule("M5  MAP - fig_M5_typology.png")
TPAL <- c("Prosperous core" = "#2166ac", "Sheltered fringe" = "#92c5de",
          "Agglomeration shadow" = "#b2182b", "Lagging region" = "#ef8a62")
TPAL <- TPAL[names(TPAL) %in% unique(type)]
png(file.path(here, "fig_M5_typology.png"), width = 2200, height = 1100, res = 175)
par(mar = c(0.4, 0.4, 3.2, 0.4))
plot(st_geometry(java), col = TPAL[type], border = "grey55", lwd = 0.3)
title(main = "M5  Typology of poverty position  (augmented Queen, all 119 units)",
      adj = 0, font.main = 1, cex.main = 1.25)
mtext("Moran quadrant on Y (signs only)  x  HDI differential  X2 - W.X2",
      side = 3, line = 0.3, adj = 0, cex = 0.8, col = "grey30")
legend("bottomleft", bty = "n", cex = 0.95, fill = TPAL, border = "grey55",
       legend = paste0(names(TPAL), "  (n=", as.integer(table(type)[names(TPAL)]), ")"))
dev.off(); log("wrote fig_M5_typology.png")

###############################################################################
# M6.  SPILLOVER-AWARE TARGETING
###############################################################################
rule("M6  TARGETING - net regional benefit of a +1 HDI unit in regency i")
imp <- read.csv(file.path(here, "M4b_impacts_with_CI.csv"), stringsAsFactors = FALSE)
b_dir <- abs(imp$semi_elasticity[imp$variable == "X2" & imp$effect == "DIRECT"])
th    <-     imp$semi_elasticity[imp$variable == "X2" & imp$effect == "INDIRECT"]
th_lo <-     imp$ci_lo_pct[imp$variable == "X2" & imp$effect == "INDIRECT"] / 100
th_hi <-     imp$ci_hi_pct[imp$variable == "X2" & imp$effect == "INDIRECT"] / 100
log("SDEM homogeneous impacts on log(Y): direct b = -", round(b_dir, 5),
    " ; theta (W.X2) = +", round(th, 5), " ; 95% CI [", round(th_lo, 5), ", ",
    round(th_hi, 5), "] (p = 0.053)")
log("impacts homogeneous (no rho) -> nothing to rank by effect size; net_i varies")
log("only through the neighbourhood (w_ji = 1/k_j, j's row weight on i)")

Y <- d$Y
own_gain <- b_dir * Y
nbr_loss <- function(theta) theta * as.vector(t(Wm) %*% Y)   # sum_j theta*Y_j*w_ji
net_of   <- function(theta) own_gain - nbr_loss(theta)
net <- net_of(th)
R0  <- rank(-Y,   ties.method = "min")   # 1 = highest poverty rate
R1  <- rank(-net, ties.method = "min")   # 1 = highest net regional benefit
rank_change <- R0 - R1                   # >0 promoted by R1 ; <0 demoted

m6 <- data.frame(bps = d$bps, name = d$name, prov = d$prov, type = d$type,
  Y = round(Y, 2), k = kdeg, own_gain = round(own_gain, 4),
  nbr_loss = round(nbr_loss(th), 4), net = round(net, 4),
  R0_rank_by_Y = R0, R1_rank_by_net = R1, rank_change_R0_minus_R1 = rank_change)
m6 <- m6[order(m6$R1_rank_by_net), ]
write.csv(m6, file.path(here, "M6_rankings.csv"), row.names = FALSE)

sp <- suppressWarnings(cor(R0, R1, method = "spearman"))
log("Spearman(R0, R1) = ", round(sp, 4),
    if (sp > 0.98) " -> near-identical; the spillover correction changes almost nothing (that is the finding)"
    else " -> the spillover correction re-orders the ranking (the divergence is the contribution)")
log("top 10 PROMOTED by the spillover-aware ranking (higher priority under R1):")
logcap(head(m6[order(-m6$rank_change_R0_minus_R1),
    c("name", "prov", "type", "Y", "k", "net", "R0_rank_by_Y", "R1_rank_by_net",
      "rank_change_R0_minus_R1")], 10))
log("top 10 DEMOTED (lower priority under R1):")
logcap(head(m6[order(m6$rank_change_R0_minus_R1),
    c("name", "prov", "type", "Y", "k", "net", "R0_rank_by_Y", "R1_rank_by_net",
      "rank_change_R0_minus_R1")], 10))
log("direction: DEMOTES high-poverty units with many/poor neighbours (positive HDI",
    " spillover makes helping them 'cost' the neighbourhood); PROMOTES weakly-",
    "connected enclave cities (k=1-2)")

## mandatory sensitivity on theta
rule("M6  SENSITIVITY - theta = 0, point, 95% CI bounds")
scen <- list(`theta=0 (no spillover)` = 0, `theta=CI low` = th_lo,
             `theta=point (0.0214)` = th, `theta=CI high` = th_hi)
sens <- do.call(rbind, lapply(names(scen), function(nm) {
  Rk <- rank(-net_of(scen[[nm]]), ties.method = "min")
  data.frame(scenario = nm, theta = round(scen[[nm]], 5),
    spearman_vs_R0 = round(cor(R0, Rk, method = "spearman"), 4),
    spearman_vs_R1pt = round(cor(R1, Rk, method = "spearman"), 4),
    max_abs_rank_move_vs_R0 = max(abs(R0 - Rk)),
    mean_abs_rank_move_vs_R0 = round(mean(abs(R0 - Rk)), 2),
    n_moving_gt10_vs_R0 = sum(abs(R0 - Rk) > 10), row.names = NULL)
}))
logcap(sens)
write.csv(sens, file.path(here, "M6_sensitivity.csv"), row.names = FALSE)
log("theta = 0 makes net_i = b*Y_i -> R1(theta=0) == R0 exactly (Spearman 1)")
log("top-10 priority list under each theta scenario:")
for (nm in names(scen)) {
  Rk <- rank(-net_of(scen[[nm]]), ties.method = "min")   # original unit order
  log("  [", nm, "]  ", paste(d$name[order(Rk)][1:10], collapse = " > "))
}
max_swing <- max(sens$max_abs_rank_move_vs_R0)
log("stability across the theta CI: max single-unit move vs R0 = ", max_swing,
    " places; units moving > 10 places ", min(sens$n_moving_gt10_vs_R0), "-",
    max(sens$n_moving_gt10_vs_R0), " of 119")
log(if (max_swing <= 10) ">>> ranking STABLE across the theta 95% CI - a league table is defensible"
    else ">>> ranking NOT stable across the theta CI - report the FRAMEWORK and DIRECTION, not a league table")

## rate vs headcount
rule("M6  RATE vs HEADCOUNT")
hc_cols <- grep("penduduk|population|jiwa|headcount|miskin_ribu|pop_",
                names(d), ignore.case = TRUE, value = TRUE)
if (length(hc_cols)) {
  log("found candidate headcount/population column(s): ", paste(hc_cols, collapse = ", "))
} else {
  log("NO headcount/population column in code/ -> M6 run in RATES ONLY; headcount")
  log("ranking PENDING. A rate ranking treats Kota Depok's 2.31% (large pop) and")
  log("Sampang's 20.61% (small pop) as the same policy quantity. If wanted, the")
  log("input is BPS 'Jumlah Penduduk Miskin (Ribu Jiwa) Menurut Kabupaten/Kota'.")
  log("Nothing downloaded in this run.")
}

rule("M5 + M6 COMPLETE - STOP")
log("M6 framing (ranking the paper stands behind vs illustrative demo) = Anita's call")
close(log_con)
cat("\nDONE - see ", OUT_LOG, "\n")
