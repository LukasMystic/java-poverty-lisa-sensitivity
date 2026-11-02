###############################################################################
# Paper 3 - spatial_M4b.R : RE-EXAMINATION of the spillover-model specification
#   seed 20260826  |  reads java_merged.gpkg (from spatial_M1_M3.R)
#
# Re-opens the M4 choice (SDM). The decision rule is written to
# M4b_decision_rule.txt BEFORE any new model is fitted and applied even if it
# lands on worse p-values / AIC. Adds SLX; full set OLS/SLX/SAR/SEM/SDM/SDEM.
# Diagnoses the cross-scheme instability of the indirect effects.
#
# Outputs (code/): M4b_decision_rule.txt, M4b_model_comparison.csv,
#   M4b_impacts_with_CI.csv, M4b_instability_diagnosis.csv, M4b_results.txt
# STOPS after this.
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
if (!file.exists("java_merged.gpkg")) stop("run spatial_M1_M3.R first")

OUT_LOG  <- file.path(here, "M4b_results.txt")
RULE_TXT <- file.path(here, "M4b_decision_rule.txt")
if (file.exists(OUT_LOG)) invisible(file.remove(OUT_LOG))
log_con <- file(OUT_LOG, open = "a")
log <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, log_con); flush(log_con) }
logcap <- function(x) { o <- capture.output(print(x)); cat(o, sep = "\n"); cat("\n")
  writeLines(o, log_con); flush(log_con) }
rule <- function(t = "") log("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78))
`%||%` <- function(a, b) if (is.null(a)) b else a

rule("PAPER 3  -  M4b  SPECIFICATION RE-EXAMINATION")
log("run ", as.character(Sys.time()), " | seed ", SEED,
    " | spdep ", as.character(packageVersion("spdep")),
    " | spatialreg ", as.character(packageVersion("spatialreg")))
log("response = log(Y); primary W = augmented Queen")

## ---- STEP 0: decision rule written to file BEFORE any model is fitted -----
rule_lines <- c(
"M4b DECISION RULE  -  written before fitting SLX or recomputing any estimate",
"===========================================================================",
"",
"FRAMING. We are choosing the specification the diagnostics and the data",
"structure support. We are NOT trying to make the indirect effects significant;",
"under any specification we expect them to stay weak, and that is acceptable.",
"If this rule and the p-values / AIC point different ways, the RULE wins and the",
"disagreement is stated explicitly.",
"",
"CANDIDATE SET: OLS, SLX, SAR, SEM, SDM, SDEM.",
"  SLX  : lm(logY ~ X + WX)                 - local spillovers, OLS, no rho/lambda",
"  SAR  : global lag of y (rho)",
"  SEM  : spatial error (lambda)",
"  SDM  : SAR + WX  (global multiplier + local-X)",
"  SDEM : SEM + WX  (spatial error + local-X, NO global multiplier)",
"",
"APPLY IN ORDER:",
"",
"R1. Spatial dependence present at all?",
"    If neither LMerr nor LMlag rejects at 0.05 -> OLS. Otherwise continue.",
"",
"R2. Which channel, per the ROBUST LM tests (each conditions on the other)?",
"    - robust LMerr rejects & robust LMlag does NOT -> ERROR / LOCAL channel.",
"      Global-lag models (SAR, SDM) are NOT supported by the score tests,",
"      regardless of how rho behaves once forced into the model. Restrict the",
"      primary candidates to {SEM, SLX, SDEM}.",
"    - robust LMlag rejects & robust LMerr does NOT -> GLOBAL-LAG channel.",
"      Primary candidates {SAR, SDM}.",
"    - both robust tests reject -> keep {SDM, SDEM}; let R4/R5 separate them.",
"    - neither robust test rejects -> fall back to non-robust LM + AIC; note it.",
"",
"R3. Local vs global, cross-checked against the M3 correlogram.",
"    If Moran's I on Y has decayed to non-significance by lag 3-4 (short-range),",
"    that is the signature of LOCAL spillovers and counts AGAINST a global",
"    multiplier (rho). Short-range correlogram + R2 = ERROR/LOCAL  =>  exclude",
"    SAR and SDM from the primary set.",
"",
"R4. Is plain SEM adequate?",
"    Hausman.test on SEM: if it REJECTS (SEM betas inconsistent with OLS), SEM is",
"    inadequate -> the model needs local-X (theta) terms -> choose among",
"    {SLX, SDEM}. If Hausman does NOT reject AND no theta LR test rejects -> SEM.",
"",
"R5. SLX vs SDEM  (SLX is nested in SDEM under lambda = 0):",
"    - Fit SLX by OLS; test its residuals with moran.mc (9999) under queen_aug.",
"    - SLX residuals NOT significantly autocorrelated (p >= 0.05): SLX absorbs the",
"      dependence alone -> CHOOSE SLX (Elhorst: SLX is the natural, most",
"      parsimonious starting point when spillovers are local).",
"    - SLX residuals STILL significantly autocorrelated: SLX insufficient ->",
"      CHOOSE SDEM, provided LR(SDEM vs SLX) rejects lambda = 0. If that LR does",
"      not reject, keep SLX.",
"",
"R6. Residual absorption is NECESSARY. The chosen model's residuals must show no",
"    significant spatial autocorrelation (moran.mc, 9999, p >= 0.05). If the",
"    R1-R5 choice fails this, escalate to the smallest superset that passes and",
"    flag that a higher-order / combined spec may be needed.",
"",
"R7. AIC is a TIEBREAK ONLY - to separate two models R1-R6 leave equally",
"    admissible. It never overturns R1-R6.",
"",
"PRE-REGISTERED EXPECTATION. R2 + R3 point to the error/local family; R4 (Hausman",
"already rejected SEM in M4, p = 0.0023) forces theta terms; so the choice is",
"SLX vs SDEM, decided by R5 (whether SLX's own residuals are clean). SDM is",
"expected to be EXCLUDED by R2/R3 despite its significant rho and lowest AIC.",
"",
"WHAT WOULD CHANGE THE OUTCOME:",
"  - robust LMlag significant & robust LMerr not  -> R2 admits SAR/SDM, drops",
"    the error family.",
"  - correlogram on Y still significant at lag 5-6 -> R3 no longer argues against",
"    a global multiplier.",
"  - Hausman does NOT reject SEM and no theta LR rejects -> R4 selects plain SEM.",
"  - SLX residuals come back clean -> R5 selects SLX over SDEM even though SDEM",
"    will have the better AIC.",
"",
paste0("Written: ", as.character(Sys.time())))
writeLines(rule_lines, RULE_TXT)
log("decision rule -> ", basename(RULE_TXT), " (", length(rule_lines),
    " lines) BEFORE any model below")

## ---- STEP 1: data + weight schemes -------------------------------------
rule("1.  DATA AND WEIGHT SCHEMES")
java <- st_read(file.path(here, "java_merged.gpkg"), quiet = TRUE)
java <- java[order(java$bps), ]
jdf  <- st_drop_geometry(java)
jdf$logY <- base::log(jdf$Y)                    # log() is shadowed by the helper
FRM  <- logY ~ X1 + X2 + X3

aeqd <- "+proj=aeqd +lat_0=-7.5 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
java_p <- st_transform(java, aeqd)
coords <- st_coordinates(st_point_on_surface(st_geometry(java_p)))
qnb <- poly2nb(java_p, queen = TRUE, row.names = as.character(java$bps))
add_sym_link <- function(nb, i, j) {
  put <- function(v, x) sort(unique(as.integer(c(v[v != 0L], x))))
  nb[[i]] <- put(nb[[i]], j); nb[[j]] <- put(nb[[j]], i); nb
}
ix <- function(b) which(java$bps == b)
AUG <- rbind(c(3526L, 3578L), c(3101L, 3175L))
qanb <- qnb
for (r in seq_len(nrow(AUG))) qanb <- add_sym_link(qanb, ix(AUG[r, 1]), ix(AUG[r, 2]))
qanb <- sym.attr.nb(qanb); attr(qanb, "region.id") <- attr(qnb, "region.id")
knn4_nb <- knn2nb(knearneigh(coords, k = 4), row.names = as.character(java$bps), sym = FALSE)

## 4th scheme: inverse-distance within the smallest band that connects the graph
nn1  <- knn2nb(knearneigh(coords, k = 1))
dmax <- max(unlist(nbdists(nn1, coords)))
db   <- dnearneigh(coords, 0, dmax, row.names = as.character(java$bps))
while (n.comp.nb(db)$nc > 1) { dmax <- dmax * 1.15
  db <- dnearneigh(coords, 0, dmax, row.names = as.character(java$bps)) }
idw_glist <- lapply(nbdists(db, coords), function(x) 1 / (x / 1000))   # inverse km
lw_idw    <- nb2listw(db, glist = idw_glist, style = "W", zero.policy = TRUE)
LW <- list(queen_aug = nb2listw(qanb,    style = "W", zero.policy = TRUE),
           queen     = nb2listw(qnb,     style = "W", zero.policy = TRUE),
           knn4      = nb2listw(knn4_nb, style = "W", zero.policy = TRUE),
           idw       = lw_idw)
W <- LW$queen_aug
log("queen_aug ", sum(card(qanb)), " links / ", n.comp.nb(qanb)$nc, " comp / ",
    sum(card(qanb) == 0), " isolate | queen bare ", sum(card(qnb)), " links / ",
    n.comp.nb(qnb)$nc, " subgraphs | IDW band 0-", round(dmax / 1000, 1),
    " km mean ", round(mean(card(db)), 1), " nb")

morI <- function(res, listw = W, ns = 9999) {
  set.seed(SEED)
  m <- moran.mc(res, listw, nsim = ns, zero.policy = TRUE, alternative = "greater")
  c(I = unname(m$statistic), p = m$p.value)
}

## ---- STEP 2: spec tests + six candidates -----------------------------
rule("2.  SPECIFICATION TESTS AND SIX CANDIDATE MODELS")
m_ols <- lm(FRM, data = jdf)
rst <- lm.RStests(m_ols, W, test = "all", zero.policy = TRUE)
spec <- do.call(rbind, lapply(names(rst), function(nm) data.frame(
  test = nm,
  alias = c(RSerr = "LMerr", RSlag = "LMlag", adjRSerr = "robustLMerr",
            adjRSlag = "robustLMlag", SARMA = "SARMA")[nm],
  statistic = round(unname(rst[[nm]]$statistic), 4),
  p_value = signif(unname(rst[[nm]]$p.value), 5), row.names = NULL)))
spec$verdict <- ifelse(spec$p_value < 0.05, "reject H0", "n.s.")
logcap(spec)
p_arserr <- spec$p_value[spec$test == "adjRSerr"]
p_arslag <- spec$p_value[spec$test == "adjRSlag"]

fit <- list()
fit$OLS  <- m_ols
fit$SLX  <- lmSLX(FRM, data = jdf, listw = W, Durbin = TRUE, zero.policy = TRUE)
fit$SAR  <- lagsarlm(FRM,  data = jdf, listw = W, Durbin = FALSE, zero.policy = TRUE)
fit$SEM  <- errorsarlm(FRM, data = jdf, listw = W, Durbin = FALSE, zero.policy = TRUE)
fit$SDM  <- lagsarlm(FRM,  data = jdf, listw = W, Durbin = TRUE,  zero.policy = TRUE)
fit$SDEM <- errorsarlm(FRM, data = jdf, listw = W, Durbin = TRUE,  zero.policy = TRUE)
spar <- function(m) {
  if (!is.null(m$rho))
    return(list(type = "rho", val = unname(m$rho), p = 2 * pnorm(-abs(m$rho / m$rho.se))))
  if (!is.null(m$lambda))
    return(list(type = "lambda", val = unname(m$lambda), p = 2 * pnorm(-abs(m$lambda / m$lambda.se))))
  list(type = "-", val = NA_real_, p = NA_real_)
}
cmp <- do.call(rbind, lapply(names(fit), function(nm) {
  m <- fit[[nm]]; sp <- spar(m); rI <- morI(residuals(m))
  data.frame(model = nm, k = attr(logLik(m), "df"),
             logLik = round(as.numeric(logLik(m)), 3), AIC = round(AIC(m), 2),
             spatial_par = sp$type,
             par_est = if (is.na(sp$val)) NA else round(sp$val, 4),
             par_p   = if (is.na(sp$p))   NA else signif(sp$p, 4),
             resid_MoranI = round(unname(rI["I"]), 4),
             resid_Moran_p = signif(unname(rI["p"]), 4),
             row.names = NULL, stringsAsFactors = FALSE)
}))
log("six-candidate comparison (residual Moran's I under queen_aug, 9999 sims):")
logcap(cmp)
write.csv(cmp, file.path(here, "M4b_model_comparison.csv"), row.names = FALSE)
for (nm in c("SLX", "SEM", "SDEM")) {
  log("[", nm, "] AIC ", round(AIC(fit[[nm]]), 2))
  co <- if (inherits(fit[[nm]], "lm") && !inherits(fit[[nm]], "Sarlm"))
    summary(fit[[nm]])$coefficients else summary(fit[[nm]])$Coef
  logcap(round(co, 5))
}

## ---- STEP 3: apply the rule -----------------------------------------
rule("3.  APPLYING THE DECISION RULE")
lr <- function(full, restr, label) {
  ll_f <- as.numeric(logLik(full)); ll_r <- as.numeric(logLik(restr))
  df <- attr(logLik(full), "df") - attr(logLik(restr), "df")
  stat <- 2 * (ll_f - ll_r); p <- pchisq(stat, df, lower.tail = FALSE)
  data.frame(comparison = label, chisq = round(stat, 4), df = df,
             p_value = signif(p, 5), row.names = NULL)
}
lrt <- rbind(
  lr(fit$SLX,  fit$OLS, "SLX vs OLS     (H0: theta = 0)"),
  lr(fit$SDEM, fit$SEM, "SDEM vs SEM    (H0: theta = 0)"),
  lr(fit$SDEM, fit$SLX, "SDEM vs SLX    (H0: lambda = 0)"),
  lr(fit$SDM,  fit$SLX, "SDM vs SLX     (H0: rho = 0)"),
  lr(fit$SDM,  fit$SAR, "SDM vs SAR     (H0: theta = 0)"),
  lr(fit$SDM,  fit$SEM, "SDM vs SEM     (H0: common factor)"))
logcap(lrt)
haus <- Hausman.test(fit$SEM)
log("Hausman on SEM: chisq ", round(unname(haus$statistic), 3), " df ",
    unname(haus$parameter), " p ", signif(unname(haus$p.value), 5),
    if (haus$p.value < 0.05) " -> REJECT (SEM inadequate)" else " -> not rejected")

slx_rI <- cmp[cmp$model == "SLX", c("resid_MoranI", "resid_Moran_p")]
lr_sdem_slx_p <- lrt$p_value[lrt$comparison == "SDEM vs SLX    (H0: lambda = 0)"]
r2_error <- p_arserr < 0.05 && p_arslag >= 0.05
r2_lag   <- p_arslag < 0.05 && p_arserr >= 0.05
log("R2 robust LMerr p ", signif(p_arserr, 3), " / robust LMlag p ", signif(p_arslag, 3),
    if (r2_error) " -> error/local channel; SAR & SDM excluded"
    else if (r2_lag) " -> global-lag channel" else " -> ambiguous")
log("R3 M3 correlogram on Y short-range (I gone by lag 4) -> against a global multiplier")
log("R4 Hausman on SEM p ", signif(unname(haus$p.value), 3),
    if (haus$p.value < 0.05) " -> SEM inadequate, need theta -> {SLX, SDEM}" else " -> SEM ok")
log("R5 SLX residual Moran's I ", slx_rI$resid_MoranI, " p ", slx_rI$resid_Moran_p,
    if (slx_rI$resid_Moran_p >= 0.05) " -> SLX clean -> SLX"
    else paste0(" -> SLX not clean -> SDEM iff LR(SDEM vs SLX) rejects (p ", signif(lr_sdem_slx_p, 3), ")"))
chosen_name <- (
  if (r2_error && haus$p.value < 0.05 && slx_rI$resid_Moran_p >= 0.05) "SLX"
  else if (r2_error && haus$p.value < 0.05 && slx_rI$resid_Moran_p < 0.05 && lr_sdem_slx_p < 0.05) "SDEM"
  else if (r2_error && haus$p.value < 0.05) "SLX"
  else if (r2_lag) "SDM"
  else "SEM"
)
chosen <- fit[[chosen_name]]
r6_rI  <- cmp[cmp$model == chosen_name, c("resid_MoranI", "resid_Moran_p")]
log("R6 chosen (", chosen_name, ") residual Moran's I ", r6_rI$resid_MoranI, " p ",
    r6_rI$resid_Moran_p, if (r6_rI$resid_Moran_p >= 0.05) " -> absorbed" else " -> NOT absorbed")
if (r6_rI$resid_Moran_p < 0.05 && chosen_name == "SLX" &&
    cmp$resid_Moran_p[cmp$model == "SDEM"] >= 0.05) {
  log("   -> escalating SLX -> SDEM (SLX fails R6, SDEM passes)")
  chosen_name <- "SDEM"; chosen <- fit$SDEM
}
aic_best <- cmp$model[cmp$model != "OLS"][which.min(cmp$AIC[cmp$model != "OLS"])]
log("R7 lowest-AIC spatial model ", aic_best, " (AIC ",
    min(cmp$AIC[cmp$model != "OLS"]), ") - tiebreak only")
log(">>> M4b CHOSEN MODEL: ", chosen_name,
    if (chosen_name != aic_best) paste0("  (differs from lowest-AIC ", aic_best,
      ", delta ", round(cmp$AIC[cmp$model == chosen_name] - min(cmp$AIC[cmp$model != "OLS"]), 2),
      "; rule wins, AIC is tiebreak only)") else "")

## ---- STEP 4: instability diagnosis --------------------------------
rule("4.  CROSS-SCHEME INSTABILITY OF THE INDIRECT EFFECTS")
fit_spec <- function(nm, data, listw) switch(nm,
  SLX  = lmSLX(FRM, data = data, listw = listw, Durbin = TRUE, zero.policy = TRUE),
  SDEM = errorsarlm(FRM, data = data, listw = listw, Durbin = TRUE, zero.policy = TRUE),
  SDM  = lagsarlm(FRM, data = data, listw = listw, Durbin = TRUE, zero.policy = TRUE),
  SEM  = errorsarlm(FRM, data = data, listw = listw, Durbin = FALSE, zero.policy = TRUE),
  SAR  = lagsarlm(FRM, data = data, listw = listw, Durbin = FALSE, zero.policy = TRUE))

## impacts + normal-approx 95% CI (LagImpact and WXimpact)
imp_ci <- function(model, listw, scenario, R = 2000) {
  tr <- trW(as(as_dgRMatrix_listw(listw), "CsparseMatrix"), type = "mult")
  set.seed(SEED)
  s <- summary(impacts(model, tr = tr, R = R), zstats = TRUE, short = TRUE)
  est <- if (!is.null(s$res)) cbind(s$res$direct, s$res$indirect, s$res$total) else s$mat
  se  <- s$semat
  vn  <- sub(" .*$", "", rownames(se))
  z95 <- qnorm(0.975)
  do.call(rbind, lapply(seq_len(nrow(se)), function(i) do.call(rbind, lapply(1:3, function(j) {
    e <- est[i, j]; s_ <- se[i, j]
    data.frame(scenario = scenario, variable = vn[i],
               effect = c("DIRECT", "INDIRECT", "TOTAL")[j],
               semi_elasticity = round(e, 6), pct_per_unit = round(100 * e, 4),
               se = round(s_, 6),
               ci_lo_pct = round(100 * (e - z95 * s_), 4),
               ci_hi_pct = round(100 * (e + z95 * s_), 4),
               p_value = signif(s$pzmat[i, j], 5), row.names = NULL)
  }))))
}

## (4a) bare Queen n=119 (disconnected) vs isolate dropped (n=118)
q_full_imp <- imp_ci(fit_spec(chosen_name, jdf, LW$queen), LW$queen, "queen_bare_n119")
keep  <- java$bps != 3101L
jdf_s <- jdf[keep, ]
qnb_s <- poly2nb(java_p[keep, ], queen = TRUE, row.names = as.character(java$bps[keep]))
lw_qs <- nb2listw(qnb_s, style = "W", zero.policy = TRUE)
log("bare Queen minus Kepulauan Seribu: n ", sum(keep), " / ", n.comp.nb(qnb_s)$nc,
    " subgraphs (Madura still detached) / ", sum(card(qnb_s) == 0), " isolates")
q_drop_imp <- imp_ci(fit_spec(chosen_name, jdf_s, lw_qs), lw_qs, "queen_bare_n118_noKepSeribu")

scen <- rbind(
  imp_ci(chosen, LW$queen_aug, "queen_aug_CHOSEN"),
  q_full_imp, q_drop_imp,
  imp_ci(fit_spec(chosen_name, jdf, LW$knn4), LW$knn4, "knn4"),
  imp_ci(fit_spec(chosen_name, jdf, LW$idw),  LW$idw,  "idw_distband"))
write.csv(scen, file.path(here, "M4b_instability_diagnosis.csv"), row.names = FALSE)

log("INDIRECT (spillover) effect across weight scenarios [", chosen_name, "]:")
ind <- scen[scen$effect == "INDIRECT",
            c("scenario", "variable", "pct_per_unit", "ci_lo_pct", "ci_hi_pct", "p_value")]
logcap(ind)
x1a <- ind$pct_per_unit[ind$scenario == "queen_aug_CHOSEN"           & ind$variable == "X1"]
x1f <- ind$pct_per_unit[ind$scenario == "queen_bare_n119"            & ind$variable == "X1"]
x1d <- ind$pct_per_unit[ind$scenario == "queen_bare_n118_noKepSeribu" & ind$variable == "X1"]
log("(4a) X1 indirect: queen_aug ", round(x1a, 3), " | bare Queen n119 ", round(x1f, 3),
    " | bare Queen n118 (no Kep. Seribu) ", round(x1d, 3), " -> ",
    if (!is.na(x1d) && abs(x1d - x1a) < abs(x1f - x1a))
      "dropping the isolate moves it toward queen_aug (disconnection was inflating it)"
    else "the difference is not the isolate")
ch <- scen[scen$scenario == "queen_aug_CHOSEN", ]
ch$ci_width_pct <- ch$ci_hi_pct - ch$ci_lo_pct
log("(4b) CI widths under the chosen model (queen_aug):")
logcap(ch[, c("variable", "effect", "pct_per_unit", "ci_lo_pct", "ci_hi_pct",
              "ci_width_pct", "p_value")])
log("mean INDIRECT 95% CI width ",
    round(mean(ch$ci_width_pct[ch$effect == "INDIRECT"]), 2),
    " pp -> low power at n=119, not sign-flip fragility per se")

## ---- STEP 5: impacts under the winner (with CIs) -----------------
rule(paste0("5.  IMPACTS UNDER THE CHOSEN MODEL (", chosen_name, ")"))
win <- scen[scen$scenario == "queen_aug_CHOSEN",
            c("variable", "effect", "semi_elasticity", "pct_per_unit",
              "se", "ci_lo_pct", "ci_hi_pct", "p_value")]
win$model <- chosen_name
write.csv(win, file.path(here, "M4b_impacts_with_CI.csv"), row.names = FALSE)
logcap(win)
log("semi-elasticities: percent change in the poverty rate per one-unit change.")
if (chosen_name %in% c("SLX", "SDEM"))
  log(chosen_name, " has no rho: INDIRECT = theta exactly, TOTAL = beta + theta; the spillover is a LOCAL quantity")

rule("M4b COMPLETE - STOP")
log("chosen specification: ", chosen_name, " (was SDM in M4)")
close(log_con)
cat("\nDONE - see ", OUT_LOG, "\n")
