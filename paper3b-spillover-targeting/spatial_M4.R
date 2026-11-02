###############################################################################
# Paper 3 - spatial_M4.R : SPILLOVER MODEL   |   seed 20260826
# Locked (Anita, 2026-08-27): response log(Y); scope Java 119; primary W =
# augmented Queen. Reads java_merged.gpkg (from spatial_M1_M3.R); rebuilds W.
#
# 1 spec tests (lm.RStests)  2 candidates OLS/SAR/SEM/SDM/SDEM (ML)
# 3 model choice (LR + Hausman, rule-led)  4 theta (W.X) test on the chosen model
# 5 impacts direct/indirect/total (+p)  6 residual moran.mc  7 robustness (queen, kNN4)
#
# Outputs (code/): M4_specification_tests.csv, M4_model_comparison.csv,
#   M4_coefficients.csv, M4_impacts.csv, M4_robustness_impacts.csv, M4_results.txt
# STOPS after M4.  (Superseded by spatial_M4b.R, which re-decides SDM -> SDEM.)
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
MERGED  <- file.path(here, "java_merged.gpkg")
OUT_LOG <- file.path(here, "M4_results.txt")
if (!file.exists(MERGED)) stop("run spatial_M1_M3.R first: ", MERGED, " not found")

if (file.exists(OUT_LOG)) invisible(file.remove(OUT_LOG))
log_con <- file(OUT_LOG, open = "a")
log <- function(...) { msg <- paste0(...); cat(msg, "\n")
  writeLines(msg, log_con); flush(log_con) }
logcap <- function(x) { out <- capture.output(print(x))
  cat(out, sep = "\n"); cat("\n"); writeLines(out, log_con); flush(log_con) }
rule <- function(t = "") log("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78))
`%||%` <- function(a, b) if (is.null(a)) b else a

rule("PAPER 3  -  M4  SPILLOVER MODEL")
log("run ", as.character(Sys.time()), " | seed ", SEED,
    " | R ", R.version$major, ".", R.version$minor,
    " | spdep ", as.character(packageVersion("spdep")),
    " | spatialreg ", as.character(packageVersion("spatialreg")))
log("response = log(Y); scope Java 119; primary W = augmented Queen")

## ---- load + rebuild weights (as in spatial_M1_M3.R) ----------------------
java <- st_read(MERGED, quiet = TRUE)
java <- java[order(java$bps), ]
stopifnot(nrow(java) == 119, !anyNA(java$Y), all(java$Y > 0))
jdf <- st_drop_geometry(java)
jdf$logY <- base::log(jdf$Y)          # log() is shadowed by the helper
FRM <- logY ~ X1 + X2 + X3

aeqd <- "+proj=aeqd +lat_0=-7.5 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
java_p <- st_transform(java, aeqd)
coords <- st_coordinates(st_point_on_surface(st_geometry(java_p)))
queen_nb <- poly2nb(java_p, queen = TRUE, row.names = as.character(java$bps))
add_sym_link <- function(nb, i, j) {
  put <- function(v, x) sort(unique(as.integer(c(v[v != 0L], x))))
  nb[[i]] <- put(nb[[i]], j); nb[[j]] <- put(nb[[j]], i); nb
}
idx <- function(b) { k <- which(java$bps == b); stopifnot(length(k) == 1); k }
AUG_LINKS <- rbind(c(3526L, 3578L), c(3101L, 3175L))
queen_aug_nb <- queen_nb
for (r in seq_len(nrow(AUG_LINKS)))
  queen_aug_nb <- add_sym_link(queen_aug_nb, idx(AUG_LINKS[r, 1]), idx(AUG_LINKS[r, 2]))
queen_aug_nb <- sym.attr.nb(queen_aug_nb)
attr(queen_aug_nb, "region.id") <- attr(queen_nb, "region.id")
knn4_nb <- knn2nb(knearneigh(coords, k = 4), row.names = as.character(java$bps), sym = FALSE)
LW <- list(queen_aug = nb2listw(queen_aug_nb, style = "W", zero.policy = TRUE),
           queen     = nb2listw(queen_nb,     style = "W", zero.policy = TRUE),
           knn4      = nb2listw(knn4_nb,      style = "W", zero.policy = TRUE))
W <- LW$queen_aug
log("weights: queen_aug ", sum(card(queen_aug_nb)), " links / ",
    n.comp.nb(queen_aug_nb)$nc, " comp / ", sum(card(queen_aug_nb) == 0), " isolate")

m_ols <- lm(FRM, data = jdf)
log("OLS log(Y)~X1+X2+X3: adj R2 ", round(summary(m_ols)$adj.r.squared, 4),
    " logLik ", round(as.numeric(logLik(m_ols)), 3), " AIC ", round(AIC(m_ols), 2),
    " k ", attr(logLik(m_ols), "df"))

###############################################################################
# 1.  SPECIFICATION TESTS  (lm.RStests = old lm.LMtests; queen_aug)
###############################################################################
rule("1.  SPECIFICATION TESTS (RSerr=LMerr, RSlag=LMlag, adjRS*=robust, SARMA)")
rst <- lm.RStests(m_ols, W, test = "all", zero.policy = TRUE)
spec <- do.call(rbind, lapply(names(rst), function(nm) data.frame(
  test = nm,
  alias = c(RSerr = "LMerr", RSlag = "LMlag", adjRSerr = "robustLMerr",
            adjRSlag = "robustLMlag", SARMA = "SARMA")[nm],
  statistic = round(unname(rst[[nm]]$statistic), 4),
  df = unname(rst[[nm]]$parameter),
  p_value = signif(unname(rst[[nm]]$p.value), 5), stringsAsFactors = FALSE)))
spec$verdict <- ifelse(spec$p_value < 0.05, "reject H0", "n.s.")
rownames(spec) <- NULL
logcap(spec)
write.csv(spec, file.path(here, "M4_specification_tests.csv"), row.names = FALSE)
p_arserr <- spec$p_value[spec$test == "adjRSerr"]
p_arslag <- spec$p_value[spec$test == "adjRSlag"]
log("robust LMerr p=", signif(p_arserr, 3), " robust LMlag p=", signif(p_arslag, 3),
    if (p_arserr < 0.05 && p_arslag >= 0.05) " -> error/local channel (tree: SEM)"
    else if (p_arslag < 0.05 && p_arserr >= 0.05) " -> lag channel (tree: SAR)"
    else " -> ambiguous")
log("LM tree does not test lag-of-X; Durbin family still estimated (see M3b correlogram)")

###############################################################################
# 2.  CANDIDATE MODELS  (ML, queen_aug)
###############################################################################
rule("2.  CANDIDATE MODELS (ML, log(Y) ~ X1+X2+X3, queen_aug)")
fit <- list()
fit$OLS  <- m_ols
fit$SAR  <- lagsarlm(FRM,  data = jdf, listw = W, Durbin = FALSE, method = "eigen", zero.policy = TRUE)
fit$SEM  <- errorsarlm(FRM, data = jdf, listw = W, Durbin = FALSE, method = "eigen", zero.policy = TRUE)
fit$SDM  <- lagsarlm(FRM,  data = jdf, listw = W, Durbin = TRUE,  method = "eigen", zero.policy = TRUE)
fit$SDEM <- errorsarlm(FRM, data = jdf, listw = W, Durbin = TRUE,  method = "eigen", zero.policy = TRUE)
get_spar <- function(m) {
  if (inherits(m, "lm")) return(c(par = NA_real_, type = NA_character_))
  if (!is.null(m$rho))    return(c(par = unname(m$rho),    type = "rho"))
  if (!is.null(m$lambda)) return(c(par = unname(m$lambda), type = "lambda"))
  c(par = NA_real_, type = NA_character_)
}
cmp <- do.call(rbind, lapply(names(fit), function(nm) {
  m <- fit[[nm]]; sp <- get_spar(m)
  data.frame(model = nm, k = attr(logLik(m), "df"),
             logLik = round(as.numeric(logLik(m)), 3), AIC = round(AIC(m), 2),
             spatial_par = ifelse(is.na(sp["type"]), "-", sp["type"]),
             estimate = ifelse(is.na(as.numeric(sp["par"])), NA, round(as.numeric(sp["par"]), 4)),
             row.names = NULL, stringsAsFactors = FALSE)
}))
logcap(cmp)
write.csv(cmp, file.path(here, "M4_model_comparison.csv"), row.names = FALSE)

for (nm in names(fit)) {
  m <- fit[[nm]]
  if (inherits(m, "lm")) { log("[", nm, "] AIC ", round(AIC(m), 2))
    logcap(round(summary(m)$coefficients, 5))
  } else {
    sm <- summary(m)
    sp <- if (!is.null(m$rho)) paste0(" rho=", round(m$rho, 4)) else paste0(" lambda=", round(m$lambda, 4))
    log("[", nm, "] AIC ", round(AIC(m), 2), sp, " (Wald p ", signif(sm$Wald1$p.value, 3), ")")
    logcap(round(sm$Coef, 5))
  }
}

###############################################################################
# 3.  MODEL CHOICE  (nested LR tests + Hausman; rule-led, not AIC alone)
###############################################################################
rule("3.  MODEL CHOICE")
lr_row <- function(label, full, restr) {
  t <- LR.Sarlm(full, restr)
  data.frame(comparison = label, chisq = round(unname(t$statistic), 4),
             df = unname(t$parameter), p_value = signif(unname(t$p.value), 5), row.names = NULL)
}
lrt <- rbind(
  lr_row("SDM vs SAR  (H0: theta = 0)",                       fit$SDM,  fit$SAR),
  lr_row("SDM vs SEM  (H0: common factor, theta=-rho*beta)",  fit$SDM,  fit$SEM),
  lr_row("SDEM vs SEM (H0: theta = 0)",                       fit$SDEM, fit$SEM),
  lr_row("SDM vs OLS",                                        fit$SDM,  fit$OLS),
  lr_row("SDEM vs OLS",                                       fit$SDEM, fit$OLS))
logcap(lrt)
haus <- tryCatch(Hausman.test(fit$SEM), error = function(e) NULL)
if (!is.null(haus))
  log("Hausman on SEM: chisq ", round(unname(haus$statistic), 4), " df ",
      unname(haus$parameter), " p ", signif(unname(haus$p.value), 5),
      if (haus$p.value < 0.05) " -> REJECT (SEM inadequate)" else " -> not rejected")

## Rule: (b) common-factor LR rejected -> keep Durbin; (c) theta LR rejected ->
## keep W.X terms; (d) SDM vs SDEM on rho significance + AIC within 2; (e) SDEM
## wins ties given the short-range correlogram.
cf_reject   <- lrt$p_value[lrt$comparison == "SDM vs SEM  (H0: common factor, theta=-rho*beta)"] < 0.05
theta_sdm   <- lrt$p_value[lrt$comparison == "SDM vs SAR  (H0: theta = 0)"] < 0.05
theta_sdem  <- lrt$p_value[lrt$comparison == "SDEM vs SEM (H0: theta = 0)"] < 0.05
rho_sdm_p   <- summary(fit$SDM)$Wald1$p.value
rho_sdm_sig <- !is.na(rho_sdm_p) && rho_sdm_p < 0.05
aic_pref    <- cmp$model[cmp$model != "OLS"][which.min(cmp$AIC[cmp$model != "OLS"])]
log("common-factor rejected? ", cf_reject, " | theta!=0 SDMvsSAR ", theta_sdm,
    " SDEMvsSEM ", theta_sdem, " | rho(SDM) sig ", rho_sdm_sig,
    " (p ", signif(rho_sdm_p, 3), ") | lowest-AIC ", aic_pref)
chosen_name <- if (!cf_reject && !theta_sdm && !theta_sdem) "SEM" else
  if (rho_sdm_sig && cmp$AIC[cmp$model == "SDM"] <= cmp$AIC[cmp$model == "SDEM"] + 2) "SDM" else "SDEM"
log(">>> CHOSEN MODEL: ", chosen_name,
    if (chosen_name == aic_pref) "  (tests and AIC agree)"
    else paste0("  (differs from lowest-AIC ", aic_pref, "; rule is test-led)"))
log("SDM vs SDEM AIC gap ", round(abs(cmp$AIC[cmp$model == "SDM"] - cmp$AIC[cmp$model == "SDEM"]), 2),
    " - SDEM stays a defensible alternative")
chosen <- fit[[chosen_name]]

###############################################################################
# 4.  THETA TEST  (W.X coefficients on the chosen model)
###############################################################################
rule(paste0("4.  THETA TEST (", chosen_name, ")"))
Coef <- summary(chosen)$Coef
rn <- rownames(Coef)
lag_rows <- grep("^lag\\.", rn)
if (!length(lag_rows)) {
  log("chosen model has no lagged-X terms - theta test N/A")
} else {
  theta_tab <- data.frame(term = sub("^lag\\.", "W.", rn[lag_rows]),
    estimate = round(Coef[lag_rows, 1], 5), std_error = round(Coef[lag_rows, 2], 5),
    z_value = round(Coef[lag_rows, 3], 3), p_value = signif(Coef[lag_rows, 4], 5), row.names = NULL)
  logcap(theta_tab)
  log("W.X1 = neighbours' TPT ; W.X2 = neighbours' HDI ; W.X3 = neighbours' TPAK")
  sig_theta <- theta_tab$term[theta_tab$p_value < 0.05]
  sp_nm  <- if (!is.null(chosen$rho)) "rho" else "lambda"
  sp_val <- if (!is.null(chosen$rho)) chosen$rho else chosen$lambda
  sp_p   <- summary(chosen)$Wald1$p.value
  log("significant theta: ", if (length(sig_theta)) paste(sig_theta, collapse = ", ") else "NONE",
      " | ", sp_nm, " = ", round(sp_val, 4), " (Wald p ", signif(sp_p, 3), ")")
  drv_sig  <- any(c("W.X1", "W.X3") %in% sig_theta)
  sp_modest <- !is.na(sp_p) && (sp_p >= 0.05 || abs(sp_val) < 0.20)
  log("scale-mismatch hypothesis (W.X1/W.X3 sig + modest ", sp_nm, "): ",
      if (length(sig_theta) == 0) "NOT supported (no theta significant)"
      else if (drv_sig && sp_modest) "SUPPORTED"
      else if (drv_sig && !sp_modest) "PARTIAL (driver lags sig but spatial par strong)"
      else paste0("NOT supported (only ", paste(sig_theta, collapse = ","),
                  " sig, ", sp_nm, " strong)"))
}
full_coef <- data.frame(term = rn, estimate = round(Coef[, 1], 6),
  std_error = round(Coef[, 2], 6), z_value = round(Coef[, 3], 4),
  p_value = signif(Coef[, 4], 5), row.names = NULL)
sp_add <- if (!is.null(chosen$rho))
  data.frame(term = "rho", estimate = round(chosen$rho, 6), std_error = round(chosen$rho.se, 6),
             z_value = round(chosen$rho / chosen$rho.se, 4),
             p_value = signif(summary(chosen)$Wald1$p.value, 5)) else
  data.frame(term = "lambda", estimate = round(chosen$lambda, 6), std_error = round(chosen$lambda.se, 6),
             z_value = round(chosen$lambda / chosen$lambda.se, 4),
             p_value = signif(summary(chosen)$Wald1$p.value, 5))
full_coef <- rbind(full_coef, sp_add)
full_coef$model <- chosen_name
write.csv(full_coef, file.path(here, "M4_coefficients.csv"), row.names = FALSE)
log("wrote M4_coefficients.csv")

###############################################################################
# 5.  IMPACT DECOMPOSITION  (semi-elasticities; INDIRECT = the spillover)
###############################################################################
rule("5.  IMPACT DECOMPOSITION (trW type='mult', impacts R=1000)")
make_impacts <- function(model, listw, seed = SEED, R = 1000) {
  if (is.null(model$rho) && !any(grepl("^lag\\.", rownames(summary(model)$Coef)))) {
    Cf <- summary(model)$Coef
    b  <- Cf[!rownames(Cf) %in% "(Intercept)", , drop = FALSE]
    vs <- rownames(b)
    return(data.frame(variable = rep(vs, 3),
      effect = rep(c("DIRECT", "INDIRECT", "TOTAL"), each = length(vs)),
      semi_elasticity = round(c(b[, 1], rep(0, length(vs)), b[, 1]), 6),
      pct_per_unit = round(100 * c(b[, 1], rep(0, length(vs)), b[, 1]), 4),
      p_value = signif(c(b[, 4], rep(NA, length(vs)), b[, 4]), 5), row.names = NULL))
  }
  tr <- trW(as(as_dgRMatrix_listw(listw), "CsparseMatrix"), type = "mult")
  set.seed(seed)
  s <- summary(impacts(model, tr = tr, R = R), zstats = TRUE, short = TRUE)
  pe <- if (!is.null(s$res)) s$res else s$impacts   # LagImpact vs WXimpact
  vars <- names(pe$direct)
  data.frame(variable = rep(vars, 3),
    effect = rep(c("DIRECT", "INDIRECT", "TOTAL"), each = length(vars)),
    semi_elasticity = round(c(pe$direct, pe$indirect, pe$total), 6),
    pct_per_unit = round(100 * c(pe$direct, pe$indirect, pe$total), 4),
    p_value = signif(c(s$pzmat[, "Direct"], s$pzmat[, "Indirect"], s$pzmat[, "Total"]), 5),
    row.names = NULL)
}
imp_tab <- make_impacts(chosen, W)
imp_tab$scheme <- "queen_aug"
logcap(imp_tab[, c("variable", "effect", "semi_elasticity", "pct_per_unit", "p_value")])
write.csv(imp_tab, file.path(here, "M4_impacts.csv"), row.names = FALSE)
sig_dir <- imp_tab$variable[imp_tab$effect == "DIRECT"   & imp_tab$p_value < 0.05]
sig_ind <- imp_tab$variable[imp_tab$effect == "INDIRECT" & imp_tab$p_value < 0.05]
sig_tot <- imp_tab$variable[imp_tab$effect == "TOTAL"    & imp_tab$p_value < 0.05]
log("significant DIRECT: ",   if (length(sig_dir)) paste(sig_dir, collapse = ", ") else "NONE",
    " | INDIRECT: ",          if (length(sig_ind)) paste(sig_ind, collapse = ", ") else "NONE",
    " | TOTAL: ",             if (length(sig_tot)) paste(sig_tot, collapse = ", ") else "NONE")
if (chosen_name == "SDM" && !is.null(fit$SDEM)) {
  sdem_imp <- make_impacts(fit$SDEM, W); sdem_imp$scheme <- "queen_aug_SDEM"
  log("SDEM impacts under queen_aug (close AIC rival):")
  logcap(sdem_imp[, c("variable", "effect", "semi_elasticity", "pct_per_unit", "p_value")])
}

###############################################################################
# 6.  RESIDUAL SPATIAL DEPENDENCE  (moran.mc, 9999 sims, queen_aug)
###############################################################################
rule("6.  RESIDUAL SPATIAL DEPENDENCE")
OLS_RESID_I <- 0.266
set.seed(SEED)
mc <- moran.mc(residuals(chosen), W, nsim = 9999, zero.policy = TRUE, alternative = "greater")
set.seed(SEED)
mc_ols <- moran.mc(residuals(m_ols), W, nsim = 9999, zero.policy = TRUE, alternative = "greater")
log(chosen_name, " residuals: Moran's I ", round(unname(mc$statistic), 4),
    " perm p ", signif(mc$p.value, 4), " | OLS resid I here ",
    round(unname(mc_ols$statistic), 4), " (baseline ", OLS_RESID_I, ")")
log(if (mc$p.value < 0.05)
  "=> residual dependence NOT fully absorbed (finding: higher-order/combined spec or omitted covariate)"
  else "=> residual spatial dependence absorbed (p >= 0.05)")

###############################################################################
# 7.  ROBUSTNESS  (chosen spec under bare Queen and kNN4)
###############################################################################
rule(paste0("7.  ROBUSTNESS (", chosen_name, " under bare Queen, kNN4)"))
refit_same <- function(listw) {
  durbin <- chosen_name %in% c("SDM", "SDEM")
  if (chosen_name %in% c("SAR", "SDM"))
    lagsarlm(FRM, data = jdf, listw = listw, Durbin = durbin, method = "eigen", zero.policy = TRUE)
  else errorsarlm(FRM, data = jdf, listw = listw, Durbin = durbin, method = "eigen", zero.policy = TRUE)
}
rob <- list(queen_aug = imp_tab)
for (sc in c("queen", "knn4")) {
  m_sc <- refit_same(LW[[sc]])
  it <- make_impacts(m_sc, LW[[sc]]); it$scheme <- sc
  rob[[sc]] <- it
  log("[", sc, "] ", if (!is.null(m_sc$rho)) paste0("rho=", round(m_sc$rho, 4))
      else paste0("lambda=", round(m_sc$lambda, 4)), " AIC ", round(AIC(m_sc), 2))
  logcap(it[, c("variable", "effect", "semi_elasticity", "pct_per_unit", "p_value")])
}
rob_all <- do.call(rbind, rob)
write.csv(rob_all, file.path(here, "M4_robustness_impacts.csv"), row.names = FALSE)

wide <- reshape(rob_all[, c("variable", "effect", "scheme", "semi_elasticity")],
                idvar = c("variable", "effect"), timevar = "scheme", direction = "wide")
names(wide) <- sub("semi_elasticity\\.", "", names(wide))
pw <- reshape(rob_all[, c("variable", "effect", "scheme", "p_value")],
              idvar = c("variable", "effect"), timevar = "scheme", direction = "wide")
names(pw) <- sub("p_value\\.", "p_", names(pw))
side <- merge(wide, pw, by = c("variable", "effect"))
side <- side[order(side$effect, side$variable), ]
logcap(side)
flip <- subset(side, sign(queen_aug) != sign(queen) | sign(queen_aug) != sign(knn4))
sigchg <- subset(side, (p_queen_aug < 0.05) != (p_queen < 0.05) |
                       (p_queen_aug < 0.05) != (p_knn4 < 0.05))
log("SIGN FLIPS across schemes: ", if (nrow(flip)) "" else "none")
if (nrow(flip)) logcap(flip[, c("variable", "effect", "queen_aug", "queen", "knn4")])
log("SIGNIFICANCE CHANGES across schemes: ", if (nrow(sigchg)) "" else "none")
if (nrow(sigchg)) logcap(sigchg[, c("variable", "effect", "p_queen_aug", "p_queen", "p_knn4")])
stable_dir <- sig_dir[sig_dir %in% side$variable[side$effect == "DIRECT" &
  side$p_queen_aug < .05 & side$p_queen < .05 & side$p_knn4 < .05]]
log("stable in sign AND significance across all three schemes: ",
    if (length(stable_dir)) paste0("DIRECT ", paste(stable_dir, collapse = ", ")) else "NONE")

rule("M4 COMPLETE - STOP")
log("chosen ", chosen_name, " | resid Moran's I ", round(unname(mc$statistic), 4),
    " p ", signif(mc$p.value, 3))
log("only significant theta = W.X2; no INDIRECT effect significant; only DIRECT HDI robust")
log("=> M6 'rank by indirect effect' needs rethinking (see spatial_M4b.R)")
close(log_con)
cat("\nDONE - see ", OUT_LOG, "\n")
