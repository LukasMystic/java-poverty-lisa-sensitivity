###############################################################################
# Paper 3 - Spatial Spillovers and Priority Targeting of Poverty: Java, Indonesia
# spatial_M1_M3.R   |   seed 20260826   |   run from the code/ directory
#
# Raw files -> merged 119-unit dataset -> validation -> M1 weights -> M2 global
# autocorrelation -> M3 LISA/Gi* -> M3b correlogram -> OLS baseline. Stops before M4.
#
# Outputs (code/): java_merged.gpkg/.csv, M1_M3_results.txt, M2_global.csv,
#   M3_lisaY_all_schemes.csv, M3_lisa_predictors_primary.csv, M3_correlogram.csv,
#   OLS_baseline.csv, M3_LISA_map_augmentedQueen.png
###############################################################################

suppressPackageStartupMessages({
  library(sf); library(spdep); library(readxl); library(writexl)
})
sf::sf_use_s2(FALSE)

SEED <- 20260826L
set.seed(SEED)

## portable working dir: dir of this script (Rscript), else current dir
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else normalizePath(getwd())
})
setwd(here)

## help Sys.which() find a modern ogr2ogr for the GDB fallback (field-type-14 bug)
for (p in c("/opt/homebrew/bin", "/usr/local/bin", "/opt/local/bin"))
  if (dir.exists(p)) Sys.setenv(PATH = paste(p, Sys.getenv("PATH"), sep = ":"))

GDB     <- file.path(here, "idn_admin_boundaries.gdb")
POV_CSV <- file.path(here,
  "Percentage of Population Living Below the National Poverty Line by Regency_Municipality, 2025.csv")
IPM_DIR <- file.path(here, "IPM")
LAB_DIR <- file.path(here, "TPT+TPAK")
XLSX    <- file.path(here, "dataset.xlsx")
OUT_GPKG <- file.path(here, "java_merged.gpkg")
OUT_CSV  <- file.path(here, "java_merged.csv")
OUT_LOG  <- file.path(here, "M1_M3_results.txt")

if (!file.exists(XLSX)) stop("Run from code/: dataset.xlsx not found in ", here)

if (file.exists(OUT_LOG)) invisible(file.remove(OUT_LOG))
log_con <- file(OUT_LOG, open = "a")
log <- function(...) { msg <- paste0(...); cat(msg, "\n")
  writeLines(msg, log_con); flush(log_con) }
logcap <- function(x) { out <- capture.output(print(x))
  cat(out, sep = "\n"); cat("\n"); writeLines(out, log_con); flush(log_con) }
rule <- function(t = "") log("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78))

rule("PAPER 3  -  M1 to M3")
log("run ", as.character(Sys.time()), " | seed ", SEED, " | ", R.version.string)
log("sf/GEOS/GDAL: ", paste(sf::sf_extSoftVersion()[c("sf", "GEOS", "GDAL")], collapse = " / "))

## ---------------------------------------------------------------------------
## 0.  Build the merged dataset
## ---------------------------------------------------------------------------
## GDAL <=3.8 OpenFileGDB aborts on ESRI field type 14 in idn_admin2 (0 rows);
## fall back to a modern ogr2ogr (brew GDAL >=3.9) transcoding the layer to GPKG.
## Release package: when the full 277 MB geodatabase is absent, read the bundled
## Java-only extract of the same layer (idn_admin2_java.gpkg: all 122 Java
## features, all fields, made with ogr2ogr -where "adm1_pcode IN (ID31..ID36)").
## The filtering below is unchanged, so both inputs give the same 119 units.
JAVA_EXTRACT <- file.path(here, "idn_admin2_java.gpkg")
read_admin2 <- function(gdb) {
  if (!dir.exists(gdb) && file.exists(JAVA_EXTRACT)) {
    x <- sf::st_read(JAVA_EXTRACT, layer = "idn_admin2", quiet = TRUE)
    log("idn_admin2: read bundled Java extract ", basename(JAVA_EXTRACT), ", ", nrow(x), " features")
    return(x)
  }
  direct <- tryCatch(
    suppressWarnings(sf::st_read(gdb, layer = "idn_admin2", quiet = TRUE)),
    error = function(e) NULL)
  if (!is.null(direct) && nrow(direct) > 100) {
    log("idn_admin2: read directly from GDB, ", nrow(direct), " features"); return(direct)
  }
  cand <- c(Sys.which("ogr2ogr"), "/opt/homebrew/bin/ogr2ogr",
            "/opt/homebrew/opt/gdal/bin/ogr2ogr", "/usr/local/bin/ogr2ogr")
  ogr <- cand[cand != "" & file.exists(cand)][1]
  if (is.na(ogr)) stop("No ogr2ogr to transcode the GDB (need GDAL >= 3.9, e.g. brew install gdal).")
  gpkg <- file.path(tempdir(), "idn_admin2_transcoded.gpkg")
  if (file.exists(gpkg)) file.remove(gpkg)
  st <- system2(ogr, c("-f", "GPKG", shQuote(gpkg), shQuote(gdb), "idn_admin2"),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(gpkg)) stop("ogr2ogr transcode failed:\n", paste(st, collapse = "\n"))
  x <- sf::st_read(gpkg, quiet = TRUE)
  log("idn_admin2: transcoded via ", ogr, ", ", nrow(x), " features"); x
}
a2_raw <- read_admin2(GDB)

## Java provinces: DKI 31, Jabar 32, Jateng 33, DIY 34, Jatim 35, Banten 36
JAVA_ADM1 <- c("ID31", "ID32", "ID33", "ID34", "ID35", "ID36")
PROV_NAME <- c(ID31 = "DKI Jakarta", ID32 = "Jawa Barat", ID33 = "Jawa Tengah",
               ID34 = "DI Yogyakarta", ID35 = "Jawa Timur", ID36 = "Banten")
nm <- names(a2_raw)
adm2_pcode_col <- nm[grep("^adm2_pcode$",     nm, ignore.case = TRUE)][1]
adm1_pcode_col <- nm[grep("^adm1_pcode$",     nm, ignore.case = TRUE)][1]
adm2_en_col    <- nm[grep("^adm2_(en|name)$", nm, ignore.case = TRUE)][1]
adm1_en_col    <- nm[grep("^adm1_(en|name)$", nm, ignore.case = TRUE)][1]
stopifnot(!is.na(adm2_pcode_col), !is.na(adm1_pcode_col), !is.na(adm2_en_col))

geo <- a2_raw
geo$adm1_pcode <- as.character(geo[[adm1_pcode_col]])
geo$adm2_pcode <- as.character(geo[[adm2_pcode_col]])
geo$name_gdb   <- as.character(geo[[adm2_en_col]])
geo$prov_gdb   <- if (!is.na(adm1_en_col)) as.character(geo[[adm1_en_col]]) else NA_character_
geo <- geo[geo$adm1_pcode %in% JAVA_ADM1, ]
geo$bps  <- suppressWarnings(as.integer(sub("^ID", "", geo$adm2_pcode)))
geo$prov <- PROV_NAME[geo$adm1_pcode]

## drop pseudo-polygons (lakes/forest/reservoirs: BPS code ends 88 or 99)
pseudo <- which(geo$bps %% 100L %in% c(88L, 99L))
log("pseudo-polygons dropped (code ends 88/99): ", length(pseudo))
if (length(pseudo)) logcap(data.frame(bps = geo$bps[pseudo], name = geo$name_gdb[pseudo]))
geo <- geo[setdiff(seq_len(nrow(geo)), pseudo), ]
bad_code <- which(is.na(geo$bps) | geo$bps < 3100L | geo$bps >= 3700L)
if (length(bad_code)) { log("WARNING dropped odd BPS codes:")
  logcap(data.frame(bps = geo$bps[bad_code], name = geo$name_gdb[bad_code]))
  geo <- geo[setdiff(seq_len(nrow(geo)), bad_code), ] }
geo <- geo[order(geo$bps), ]
log("Java admin-2 units after filtering: ", nrow(geo), " (expect 119)")
logcap(as.data.frame(table(geo$prov)))
if (any(duplicated(geo$bps))) stop("duplicate BPS codes in geometry")

## poverty CSV -> Y (national file: ALL-CAPS province rows interleave kab/kota)
pov_raw <- read.csv(POV_CSV, skip = 3, header = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8-BOM", col.names = c("name", "value"))
pov_raw$name  <- trimws(pov_raw$name)
pov_raw <- pov_raw[pov_raw$name != "", ]
pov_raw$value <- suppressWarnings(as.numeric(pov_raw$value))
is_prov_hdr <- !grepl("[a-z]", pov_raw$name)      # province header = no lowercase
pov_raw$prov_hdr <- NA_character_
cur <- NA_character_
for (i in seq_len(nrow(pov_raw))) {
  if (is_prov_hdr[i]) cur <- pov_raw$name[i]; pov_raw$prov_hdr[i] <- cur
}
JAVA_POV_HDR <- c("DKI JAKARTA" = "DKI Jakarta", "JAWA BARAT" = "Jawa Barat",
                  "JAWA TENGAH" = "Jawa Tengah", "D I YOGYAKARTA" = "DI Yogyakarta",
                  "JAWA TIMUR" = "Jawa Timur", "BANTEN" = "Banten")
pov <- pov_raw[!is_prov_hdr & pov_raw$prov_hdr %in% names(JAVA_POV_HDR), ]
pov <- data.frame(prov = JAVA_POV_HDR[pov$prov_hdr], name_csv = pov$name,
                  Y = pov$value, stringsAsFactors = FALSE)
log("poverty CSV: ", nrow(pov), " Java kab/kota rows (expect 119)")
if (anyNA(pov$Y)) stop("NA poverty values after parse")

## IPM -> X2 ; TPT+TPAK (August) -> X1, X3   (each provincial CSV ends with a
## province-total row; drop it. February columns are '...' at kab/kota level.)
prov_from_file <- function(f) {
  f <- basename(f)
  if (grepl("Jawa Timur",  f)) return("Jawa Timur")
  if (grepl("Jawa Barat",  f)) return("Jawa Barat")
  if (grepl("Jawa Tengah", f)) return("Jawa Tengah")
  if (grepl("Yogyakarta",  f)) return("DI Yogyakarta")
  if (grepl("DKI|Jakarta", f)) return("DKI Jakarta")
  if (grepl("Banten",      f)) return("Banten")
  stop("cannot map file to province: ", f)
}
PROV_TOTAL_LABELS <- c("Jawa Timur", "Jawa Barat", "Jawa Tengah", "DI Yogyakarta",
  "D.I. Yogyakarta", "Daerah Istimewa Yogyakarta", "DKI Jakarta", "Banten")
drop_total <- function(d) {
  d <- d[-nrow(d), , drop = FALSE]
  d[!(d$name %in% PROV_TOTAL_LABELS), , drop = FALSE]
}
read_ipm <- function(f) {
  d <- read.csv(f, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(d) <- c("name", "X2")
  d$name <- trimws(d$name); d$X2 <- suppressWarnings(as.numeric(d$X2))
  d <- drop_total(d); d$prov <- prov_from_file(f); d[, c("prov", "name", "X2")]
}
read_lab <- function(f) {
  d <- read.csv(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
                fileEncoding = "UTF-8-BOM")
  names(d)[1] <- "name"
  tpt_ag  <- grep("Pengangguran.*Agustus|TPT.*Agustus|Agustus.*Pengangguran", names(d))
  tpak_ag <- grep("Partisipasi.*Agustus|TPAK.*Agustus|Agustus.*Partisipasi", names(d))
  stopifnot(length(tpt_ag) == 1, length(tpak_ag) == 1)
  d$name <- trimws(d$name)
  d$X1 <- suppressWarnings(as.numeric(d[[tpt_ag]]))     # TPT August
  d$X3 <- suppressWarnings(as.numeric(d[[tpak_ag]]))    # TPAK August
  d <- drop_total(d); d$prov <- prov_from_file(f); d[, c("prov", "name", "X1", "X3")]
}
ipm <- do.call(rbind, lapply(list.files(IPM_DIR, "\\.csv$", full.names = TRUE), read_ipm))
lab <- do.call(rbind, lapply(list.files(LAB_DIR, "\\.csv$", full.names = TRUE), read_lab))
log("IPM rows ", nrow(ipm), " | TPT+TPAK rows ", nrow(lab), " (expect 119 each)")
if (anyNA(ipm$X2)) stop("NA in IPM/X2")
if (anyNA(lab$X1) || anyNA(lab$X3)) stop("NA in TPT/TPAK")
if (any(ipm$X2 < 55 | ipm$X2 > 95)) log("WARNING HDI outside 55-95 -> wrong BPS table?")
log("HDI range across 119 units: ", round(min(ipm$X2), 2), " - ", round(max(ipm$X2), 2))

## join the three CSV sources, then join to geometry, on normalised names
norm_name <- function(x) {
  x <- tolower(trimws(x))
  x <- iconv(x, to = "ASCII//TRANSLIT"); x[is.na(x)] <- ""
  x <- gsub("[.'`-]", " ", x)
  x <- gsub("\\bkab(upaten)?\\b", "", x)      # keep 'kota' - distinguishes city
  x <- gsub("\\badm\\.?\\b", "", x)
  gsub("\\s+", "", x)                          # drop spaces: gunung kidul = gunungkidul
}
merge_csv <- Reduce(function(a, b) merge(a, b, by = c("prov", "name"), all = TRUE),
  list(transform(pov, name = name_csv)[, c("prov", "name", "Y")],
       setNames(ipm, c("prov", "name", "X2")),
       setNames(lab, c("prov", "name", "X1", "X3"))))
bad <- merge_csv[!complete.cases(merge_csv[, c("Y", "X1", "X2", "X3")]), ]
if (nrow(bad)) { logcap(bad); stop("CSV sources disagree on kab/kota spelling") }
if (nrow(merge_csv) != 119) stop("CSV merge is not 119 rows")

geo$key <- paste(geo$prov, norm_name(geo$name_gdb))
merge_csv$key <- paste(merge_csv$prov, norm_name(merge_csv$name))
merged <- merge(as.data.frame(geo)[, c("bps", "prov", "name_gdb", "adm2_pcode", "key")],
                merge_csv[, c("prov", "name", "key", "Y", "X1", "X2", "X3")],
                by = c("prov", "key"), all = FALSE)
if (nrow(merged) != 119) {
  logcap(geo[!geo$key %in% merged$key, c("prov", "name_gdb", "bps")])
  logcap(merge_csv[!merge_csv$key %in% merged$key, c("prov", "name", "Y")])
  stop("geometry x CSV join did not produce 119 units")
}
if (anyNA(merged[, c("Y", "X1", "X2", "X3")])) stop("NA in merged variables")
if (any(duplicated(merged$bps))) stop("duplicate BPS after join")

merged <- merged[order(merged$bps), ]
geo_ord <- geo[match(merged$bps, geo$bps), ]
java <- st_sf(bps = merged$bps, adm2_pcode = merged$adm2_pcode, name = merged$name_gdb,
              prov = merged$prov, Y = merged$Y, X1 = merged$X1, X2 = merged$X2,
              X3 = merged$X3, geometry = st_geometry(geo_ord))
st_crs(java) <- st_crs(a2_raw)
log("merged dataset: ", nrow(java), " units, all with geometry")
logcap(summary(st_drop_geometry(java)[, c("Y", "X1", "X2", "X3")]))
log("Y: mean ", round(mean(java$Y), 2), " median ", round(median(java$Y), 2),
    " sd ", round(sd(java$Y), 2), " min ", round(min(java$Y), 2),
    " (", java$name[which.min(java$Y)], ") max ", round(max(java$Y), 2),
    " (", java$name[which.max(java$Y)], ")")

## ---------------------------------------------------------------------------
## V.  Validation - 38 East Java rows must reproduce dataset.xlsx exactly
## ---------------------------------------------------------------------------
rule("VALIDATION vs dataset.xlsx (38 East Java rows)")
anchor <- as.data.frame(read_excel(XLSX))
ej <- java[java$prov == "Jawa Timur", ]
ej <- ej[order(ej$bps), ]
stopifnot(nrow(ej) == 38, nrow(anchor) == 38)   # dataset.xlsx is in BPS order
cmp <- data.frame(bps = ej$bps, name = ej$name,
  Y_new = ej$Y, Y_ref = anchor$Y, X1_new = ej$X1, X1_ref = anchor$X1,
  X2_new = ej$X2, X2_ref = anchor$X2, X3_new = ej$X3, X3_ref = anchor$X3)
tol <- 0.005
mm <- list(Y  = which(abs(cmp$Y_new  - cmp$Y_ref)  > tol),
           X1 = which(abs(cmp$X1_new - cmp$X1_ref) > tol),
           X2 = which(abs(cmp$X2_new - cmp$X2_ref) > tol),
           X3 = which(abs(cmp$X3_new - cmp$X3_ref) > tol))
n_mm <- sapply(mm, length)
log("mismatched cells (|new-ref| > ", tol, "): ",
    paste(names(n_mm), n_mm, sep = "=", collapse = " "))
if (any(n_mm > 0)) {
  for (v in names(mm)) if (length(mm[[v]]))
    logcap(cmp[mm[[v]], c("bps", "name", grep(v, names(cmp), value = TRUE))])
  close(log_con); stop("VALIDATION FAILED")
}
log("VALIDATION PASSED: 0 mismatched cells on Y, X1, X2, X3 across all 38 units")

## centroid cross-check vs dataset.xlsx u,v (multipart regencies expected off)
cent_ll <- suppressWarnings(st_coordinates(st_centroid(st_geometry(ej))))
uv <- data.frame(bps = ej$bps, name = ej$name, u_ref = anchor$u, v_ref = anchor$v,
  u_cent = round(cent_ll[, 1], 3), v_cent = round(cent_ll[, 2], 3),
  d_deg = round(sqrt((cent_ll[, 1] - anchor$u)^2 + (cent_ll[, 2] - anchor$v)^2), 3))
log("centroid vs u,v - largest 6 deviations (deg):")
logcap(head(uv[order(-uv$d_deg), ], 6))

if (file.exists(OUT_GPKG)) invisible(file.remove(OUT_GPKG))
st_write(java, OUT_GPKG, quiet = TRUE)
write.csv(st_drop_geometry(java), OUT_CSV, row.names = FALSE)
log("wrote ", basename(OUT_GPKG), " and ", basename(OUT_CSV))

###############################################################################
# M1.  SPATIAL WEIGHTS
###############################################################################
rule("M1  SPATIAL WEIGHTS")
aeqd <- "+proj=aeqd +lat_0=-7.5 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
java_p <- st_transform(java, aeqd)
coords <- st_coordinates(st_point_on_surface(st_geometry(java_p)))
row.names(coords) <- as.character(java$bps)

## bare Queen (robustness variant)
queen_nb <- poly2nb(java_p, queen = TRUE, row.names = as.character(java$bps))
queen_lw <- nb2listw(queen_nb, style = "W", zero.policy = TRUE)
log("\n-- bare Queen (style W) : ROBUSTNESS VARIANT --")
logcap(queen_nb)
card_q <- card(queen_nb)
logcap(as.data.frame(table(k = card_q)))
comp_q <- n.comp.nb(queen_nb)
log("mean links ", round(mean(card_q), 3), " | symmetric ", is.symmetric.nb(queen_nb),
    " | components ", comp_q$nc, " (sizes ",
    paste(sort(table(comp_q$comp), decreasing = TRUE), collapse = ","), ")")
iso_q <- which(card_q == 0)
if (length(iso_q)) { log("isolated (zero.policy=TRUE):")
  logcap(data.frame(bps = java$bps[iso_q], name = java$name[iso_q])) }
main_comp <- as.integer(names(which.max(table(comp_q$comp))))
off_main  <- which(comp_q$comp != main_comp)
if (length(off_main)) {
  log("not connected to the mainland component (", length(off_main), " units):")
  logcap(data.frame(bps = java$bps[off_main], name = java$name[off_main],
                    Y = java$Y[off_main], component = comp_q$comp[off_main],
                    k = card_q[off_main]))
}
log("most-connected:")
logcap(head(data.frame(bps = java$bps, name = java$name, k = card_q)[order(-card_q), ], 5))

## augmented Queen (PRIMARY): + Bangkalan(3526)-Kota Surabaya(3578) [Suramadu,
## 2009] and Kepulauan Seribu(3101)-Kota Jakarta Utara(3175) [admin/ferry link]
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
attr(queen_aug_nb, "ncomp") <- n.comp.nb(queen_aug_nb)
queen_aug_lw <- nb2listw(queen_aug_nb, style = "W", zero.policy = TRUE)
log("\n-- augmented Queen (style W) : PRIMARY --")
logcap(queen_aug_nb)
card_a <- card(queen_aug_nb); comp_a <- n.comp.nb(queen_aug_nb)
logcap(as.data.frame(table(k = card_a)))
log("mean links ", round(mean(card_a), 3), " (bare ", round(mean(card_q), 3),
    ") | symmetric ", is.symmetric.nb(queen_aug_nb),
    " | isolates ", sum(card_a == 0), " | components ", comp_a$nc)
if (!(comp_a$nc == 1 && sum(card_a == 0) == 0))
  stop("augmented Queen not fully connected - check AUG_LINKS")
log("=> fully connected, no isolates. Endpoints changing degree:")
chg <- which(card_a != card_q)
logcap(data.frame(bps = java$bps[chg], name = java$name[chg],
                  k_bare = card_q[chg], k_aug = card_a[chg]))

## kNN k = 4,5,6 (robustness only)
knn_nb <- list(); knn_lw <- list()
for (k in c(4, 5, 6)) {
  kn <- knn2nb(knearneigh(coords, k = k), row.names = as.character(java$bps), sym = FALSE)
  knn_nb[[as.character(k)]] <- kn
  knn_lw[[as.character(k)]] <- nb2listw(kn, style = "W", zero.policy = TRUE)
  log("kNN k=", k, ": symmetrised mean links ",
      round(mean(card(make.sym.nb(kn))), 3))
}
for (b in c(3525L, 3529L)) {            # Gresik (point on Bawean) / Sumenep (Queen degree 1); see checks/paper3a_seed_probe.R
  i <- which(java$bps == b)
  for (k in c("4", "5", "6")) {
    nb_i <- knn_nb[[k]][[i]]
    log("kNN k=", k, " ", java$name[i], " -> ",
        paste0(java$name[nb_i], "(", java$bps[nb_i], ")", collapse = ", "))
  }
}

W  <- list(QueenAug = queen_aug_lw, Queen = queen_lw,
           kNN4 = knn_lw[["4"]], kNN5 = knn_lw[["5"]], kNN6 = knn_lw[["6"]])
NB <- list(QueenAug = queen_aug_nb, Queen = queen_nb,
           kNN4 = knn_nb[["4"]], kNN5 = knn_nb[["5"]], kNN6 = knn_nb[["6"]])
PRIMARY <- "QueenAug"

###############################################################################
# M2.  GLOBAL AUTOCORRELATION
###############################################################################
NSIM_M2 <- 9999L
NSIM_M3 <- 99999L
rule(paste0("M2  GLOBAL AUTOCORRELATION (Moran's I & Geary's C, ", NSIM_M2, " perm)"))
set.seed(SEED)   # M2 results depend on this seed AND the loop order below
vars <- c("Y", "X1", "X2", "X3")
m2 <- data.frame()
for (v in vars) for (s in names(W)) {
  mi <- moran.mc(java[[v]], W[[s]], nsim = NSIM_M2, zero.policy = TRUE, alternative = "greater")
  gc <- geary.mc(java[[v]], W[[s]], nsim = NSIM_M2, zero.policy = TRUE, alternative = "greater")
  m2 <- rbind(m2, data.frame(variable = v, scheme = s,
    Moran_I = round(unname(mi$statistic), 4), Moran_p = signif(mi$p.value, 4),
    Geary_C = round(unname(gc$statistic), 4), Geary_p = signif(gc$p.value, 4),
    Moran_rank = unname(mi$parameter), Geary_rank = unname(gc$parameter)))
}
rownames(m2) <- NULL
## geary.mc alt='greater' is spdep's rank convention for positive autocorrelation
logcap(m2)
write.csv(m2, file.path(here, "M2_global.csv"), row.names = FALSE)
log("wrote M2_global.csv")
logcap(m2[m2$scheme %in% c(PRIMARY, "Queen"), ])

###############################################################################
# M3.  LOCAL STRUCTURE (LISA + Gi*, FDR q = 0.05)
###############################################################################
rule(paste0("M3  LOCAL STRUCTURE (localmoran_perm + localG_perm, ", NSIM_M3,
            " perm, FDR-BH q=0.05, iseed=", SEED, ")"))
Q <- 0.05
## LISA p = 'Pr(z != E(Ii)) Sim' (two-sided cond. perm); class = quadr$mean;
## FDR = Benjamini-Hochberg (Caldas de Castro & Singer 2006). NSIM_M3 raised
## from 9,999 because the FDR count on Y was seed-unstable (4-6 across the six
## probe seeds; recorded by checks/paper3a_seed_probe.R).
lisa_one <- function(v, s) {
  x <- java[[v]]; lw <- W[[s]]
  lm <- localmoran_perm(x, lw, nsim = NSIM_M3, zero.policy = TRUE,
                        alternative = "two.sided", iseed = SEED)
  pcol <- grep("^Pr\\(z .* E\\(Ii\\)\\) Sim$", colnames(lm), value = TRUE)
  if (!length(pcol)) stop("LISA sim p-value column not found")
  p_raw <- lm[, pcol[1]]
  p_fdr <- rep(NA_real_, length(p_raw))
  ok <- !is.na(p_raw)                       # isolates -> NA; BH over defined only
  p_fdr[ok] <- p.adjust(p_raw[ok], method = "BH")
  quad <- attr(lm, "quadr")$mean
  lw_self <- nb2listw(include.self(NB[[s]]), style = "W", zero.policy = TRUE)
  g <- localG_perm(x, lw_self, nsim = NSIM_M3, zero.policy = TRUE, iseed = SEED)
  gi_int <- attr(g, "internals")
  gpcol <- if (!is.null(gi_int)) grep("^Pr\\(z .*\\) Sim$", colnames(gi_int), value = TRUE) else character(0)
  gp_raw <- if (length(gpcol)) gi_int[, gpcol[1]] else 2 * pnorm(-abs(as.numeric(g)))
  gp_fdr <- rep(NA_real_, length(gp_raw))
  okg <- !is.na(gp_raw)
  gp_fdr[okg] <- p.adjust(gp_raw[okg], method = "BH")
  data.frame(bps = java$bps, name = java$name, prov = java$prov, variable = v, scheme = s,
    x = x, Ii = lm[, "Ii"], quad = as.character(quad), p_raw = p_raw, p_fdr = p_fdr,
    sig_raw = !is.na(p_raw) & p_raw < Q, sig_fdr = !is.na(p_fdr) & p_fdr < Q,
    Gi_z = as.numeric(g), Gp_raw = gp_raw, Gp_fdr = gp_fdr,
    G_sig_raw = !is.na(gp_raw) & gp_raw < Q, G_sig_fdr = !is.na(gp_fdr) & gp_fdr < Q,
    defined = ok, stringsAsFactors = FALSE)
}
lisaY <- do.call(rbind, lapply(names(W), function(s) lisa_one("Y", s)))

## seed-stability check on the headline FDR count (kept as a reproducibility record)
probe_seeds <- c(1L, 20260826L, 42L, 2024L, 99999L, 7L)
stab <- sapply(probe_seeds, function(sd) {
  lmv <- localmoran_perm(java$Y, W[[PRIMARY]], nsim = NSIM_M3, zero.policy = TRUE,
                         alternative = "two.sided", iseed = sd)
  pc <- grep("^Pr\\(z .* E\\(Ii\\)\\) Sim$", colnames(lmv), value = TRUE)[1]
  pr <- lmv[, pc]; okp <- !is.na(pr)
  pf <- rep(NA_real_, length(pr)); pf[okp] <- p.adjust(pr[okp], "BH")
  c(uncorrected = sum(pr < Q, na.rm = TRUE), FDR_q05 = sum(pf < Q, na.rm = TRUE))
})
colnames(stab) <- paste0("seed_", probe_seeds)
log("FDR-count seed stability (LISA on Y, ", PRIMARY, "):")
logcap(t(stab))

log("LISA on Y - significant-unit counts by scheme:")
tab <- do.call(rbind, lapply(names(W), function(s) {
  d <- lisaY[lisaY$scheme == s, ]
  data.frame(scheme = s, n_defined = sum(d$defined), undefined = sum(!d$defined),
             uncorrected = sum(d$sig_raw), FDR_q05 = sum(d$sig_fdr))
}))
logcap(tab)

log("LISA on Y - class breakdown (FDR-significant):")
for (s in names(W)) {
  d <- lisaY[lisaY$scheme == s & lisaY$sig_fdr, ]
  log("  ", s, ": ", nrow(d))
  if (nrow(d)) logcap(as.data.frame(table(class = d$quad)))
}
log("LISA on Y - class breakdown (uncorrected):")
for (s in names(W)) {
  d <- lisaY[lisaY$scheme == s & lisaY$sig_raw, ]
  log("  ", s, ": ", nrow(d))
  if (nrow(d)) logcap(as.data.frame(table(class = d$quad)))
}

log(PRIMARY, " - every FDR-significant LISA unit on Y:")
dq <- lisaY[lisaY$scheme == PRIMARY & lisaY$sig_fdr,
            c("bps", "name", "prov", "x", "Ii", "quad", "p_raw", "p_fdr")]
logcap(dq[order(dq$quad, -dq$Ii), ])
log("bare Queen FDR units:")
dqb <- lisaY[lisaY$scheme == "Queen" & lisaY$sig_fdr,
             c("bps", "name", "prov", "x", "Ii", "quad", "p_raw", "p_fdr")]
logcap(dqb[order(dqb$quad, -dqb$Ii), ])

log("Gi* (localG_perm, include.self) on Y - hot/cold counts:")
gtab <- do.call(rbind, lapply(names(W), function(s) {
  d <- lisaY[lisaY$scheme == s, ]
  data.frame(scheme = s,
    hot_raw = sum(d$G_sig_raw & d$Gi_z > 0), cold_raw = sum(d$G_sig_raw & d$Gi_z < 0),
    hot_fdr = sum(d$G_sig_fdr & d$Gi_z > 0), cold_fdr = sum(d$G_sig_fdr & d$Gi_z < 0))
}))
logcap(gtab)

log("LISA on predictors (", PRIMARY, ") - counts:")
lisaP <- do.call(rbind, lapply(c("X1", "X2", "X3"), function(v) lisa_one(v, PRIMARY)))
ptab <- do.call(rbind, lapply(c("X1", "X2", "X3"), function(v) {
  d <- lisaP[lisaP$variable == v, ]
  data.frame(variable = v, uncorrected = sum(d$sig_raw), FDR_q05 = sum(d$sig_fdr))
}))
logcap(ptab)

rule("MADURA SPOTLIGHT - Bangkalan (3526), Sampang (3527), Sumenep (3529)")
for (b in c(3526L, 3527L, 3529L)) {
  i <- which(java$bps == b)
  sm <- lisaY[lisaY$bps == b, ]
  log(java$name[i], " (", b, ") Y=", java$Y[i], " (Java mean ", round(mean(java$Y), 2), ")")
  logcap(sm[, c("scheme", "Ii", "quad", "p_raw", "p_fdr", "sig_raw", "sig_fdr",
                "Gi_z", "Gp_raw", "Gp_fdr", "G_sig_fdr")])
  for (sc in c("QueenAug", "Queen")) {
    nb_i <- NB[[sc]][[i]]; has <- length(nb_i) && nb_i[1] != 0
    log("  ", sc, " nb: ", if (has) paste0(java$name[nb_i], "(", java$bps[nb_i],
        ", Y=", java$Y[nb_i], ")", collapse = ", ") else "(isolated)")
  }
}
write.csv(lisaY, file.path(here, "M3_lisaY_all_schemes.csv"), row.names = FALSE)
write.csv(lisaP, file.path(here, "M3_lisa_predictors_primary.csv"), row.names = FALSE)
log("wrote M3_lisaY_all_schemes.csv, M3_lisa_predictors_primary.csv")

###############################################################################
# M3b.  SPATIAL CORRELOGRAM
###############################################################################
rule("M3b  SPATIAL CORRELOGRAM (augmented Queen nb, Moran's I, lags 1-6)")
## sp.correlogram(method='I', randomisation=TRUE); sig = two-sided z = (I-E)/sd(I)
corr_tab <- data.frame()
for (v in c("Y", "X1", "X2", "X3")) {
  sc <- sp.correlogram(queen_aug_nb, java[[v]], order = 6, method = "I",
                       style = "W", randomisation = TRUE, zero.policy = TRUE)
  res <- sc$res
  z <- (res[, 1] - res[, 2]) / sqrt(res[, 3]); p <- 2 * pnorm(-abs(z))
  corr_tab <- rbind(corr_tab, data.frame(variable = v, lag = 1:6,
    MoranI = round(res[, 1], 4), E = round(res[, 2], 4), sd = round(sqrt(res[, 3]), 4),
    z = round(z, 3), p = signif(p, 3),
    sig = ifelse(p < .01, "**", ifelse(p < .05, "*", "")), stringsAsFactors = FALSE))
}
rownames(corr_tab) <- NULL
logcap(corr_tab)
write.csv(corr_tab, file.path(here, "M3_correlogram.csv"), row.names = FALSE)
log("wrote M3_correlogram.csv")

###############################################################################
# OLS.  BASELINE + RESIDUAL DIAGNOSTICS (response Y, then log(Y))
###############################################################################
rule("OLS BASELINE + RESIDUAL DIAGNOSTICS")
jdf <- st_drop_geometry(java)
jdf$logY <- base::log(jdf$Y)        # log() is shadowed by the logging helper
preds <- c("X1", "X2", "X3")
vif <- sapply(preds, function(p) {  # VIF depends only on the predictor set
  f <- as.formula(paste(p, "~", paste(setdiff(preds, p), collapse = " + ")))
  1 / (1 - summary(lm(f, data = jdf))$r.squared)
})
log("VIF (X1,X2,X3; same for Y and log(Y) models):")
logcap(data.frame(predictor = preds, VIF = round(vif, 3)))

m_Y    <- lm(Y    ~ X1 + X2 + X3, data = jdf)   # top level so lm.morantest re-evaluates
m_logY <- lm(logY ~ X1 + X2 + X3, data = jdf)
ols_diag <- function(resp_label, m) {
  sm <- summary(m)
  log("-- response = ", resp_label, " --")
  logcap(round(sm$coefficients, 5))
  log("adj R2 ", round(sm$adj.r.squared, 4), " | resid SE ", round(sm$sigma, 4),
      " on ", m$df.residual, " df")
  sh <- shapiro.test(residuals(m))
  log("Shapiro-Wilk resid: W ", round(unname(sh$statistic), 4), " p ", signif(sh$p.value, 4))
  mt  <- lm.morantest(m, queen_aug_lw, zero.policy = TRUE, alternative = "greater")
  mtp <- moran.mc(residuals(m), queen_aug_lw, nsim = NSIM_M2, zero.policy = TRUE,
                  alternative = "greater")
  log("resid Moran's I (aug Queen): lm.morantest I ", round(unname(mt$estimate[1]), 4),
      " z ", round(unname(mt$statistic), 3), " p ", signif(mt$p.value, 4),
      " | moran.mc I ", round(unname(mtp$statistic), 4), " p ", signif(mtp$p.value, 4),
      " | paper 1 East Java ref 0.1764")
  data.frame(response = resp_label, term = rownames(sm$coefficients),
    estimate = round(sm$coefficients[, 1], 5), std_error = round(sm$coefficients[, 2], 5),
    t_value = round(sm$coefficients[, 3], 3), p_value = signif(sm$coefficients[, 4], 4),
    adj_r2 = round(sm$adj.r.squared, 4), shapiro_W = round(unname(sh$statistic), 4),
    shapiro_p = signif(sh$p.value, 4), resid_moranI = round(unname(mt$estimate[1]), 4),
    resid_moran_p_normal = signif(mt$p.value, 4), resid_moran_p_perm = signif(mtp$p.value, 4),
    row.names = NULL, stringsAsFactors = FALSE)
}
set.seed(SEED)
ols_out <- rbind(ols_diag("Y", m_Y), ols_diag("log(Y)", m_logY))
write.csv(ols_out, file.path(here, "OLS_baseline.csv"), row.names = FALSE)
log("wrote OLS_baseline.csv")

###############################################################################
# WEEK-3 DECISION GATE
###############################################################################
rule("WEEK-3 DECISION GATE")
hl <- do.call(rbind, lapply(names(W), function(s) {
  d <- lisaY[lisaY$scheme == s, ]
  cl <- function(sub, k) sum(sub & d$quad == k)
  data.frame(scheme = s,
    unc_total = sum(d$sig_raw), unc_HH = cl(d$sig_raw, "High-High"),
    unc_LL = cl(d$sig_raw, "Low-Low"), unc_HL = cl(d$sig_raw, "High-Low"),
    unc_LH = cl(d$sig_raw, "Low-High"),
    fdr_total = sum(d$sig_fdr), fdr_HH = cl(d$sig_fdr, "High-High"),
    fdr_LL = cl(d$sig_fdr, "Low-Low"), fdr_HL = cl(d$sig_fdr, "High-Low"),
    fdr_LH = cl(d$sig_fdr, "Low-High"))
}))
logcap(hl)
q_primary <- sum(lisaY$scheme == PRIMARY & lisaY$sig_fdr)
prim_units <- lisaY[lisaY$scheme == PRIMARY & lisaY$sig_fdr, ]
log("PRIMARY (", PRIMARY, ") FDR LISA clusters on Y: ", q_primary,
    " (", sum(prim_units$quad == "High-High"), " HH, ",
    sum(prim_units$quad == "Low-Low"), " LL, ",
    sum(prim_units$quad == "High-Low"), " HL, ",
    sum(prim_units$quad == "Low-High"), " LH)")
log("  HH: ", paste(prim_units$name[prim_units$quad == "High-High"], collapse = ", "))
log("  LL: ", paste(prim_units$name[prim_units$quad == "Low-Low"], collapse = ", "))
for (s in c("Queen", "kNN4", "kNN5", "kNN6"))
  log("  ", s, ": ", sum(lisaY$scheme == s & lisaY$sig_fdr), " FDR / ",
      sum(lisaY$scheme == s & lisaY$sig_raw), " uncorrected")
for (b in c(3526L, 3527L, 3529L)) {
  r <- lisaY[lisaY$scheme == PRIMARY & lisaY$bps == b, ]
  log("  ", sprintf("%-10s", r$name), " FDR-sig ", r$sig_fdr, " (", r$quad,
      ", perm p ", sprintf("%.4f", r$p_raw), ", FDR p ", sprintf("%.3f", r$p_fdr), ")")
}
log("GATE: ", q_primary, if (q_primary < 10) " < ~10 (below threshold)" else " >= ~10")

## LISA cluster map on PRIMARY weights
rule("LISA CLUSTER MAP (augmented Queen, FDR q=0.05)")
MAP_PNG <- file.path(here, "M3_LISA_map_augmentedQueen.png")
dmap <- lisaY[lisaY$scheme == PRIMARY, ]
dmap <- dmap[match(java$bps, dmap$bps), ]
cls <- ifelse(dmap$sig_fdr, as.character(dmap$quad), "Not significant")
cls <- factor(cls, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant"))
pal <- c("High-High" = "#b2182b", "Low-Low" = "#2166ac", "High-Low" = "#ef8a62",
         "Low-High" = "#67a9cf", "Not significant" = "grey88")
gll <- st_geometry(java)
cn  <- suppressWarnings(st_coordinates(st_point_on_surface(gll)))
png(MAP_PNG, width = 1800, height = 1200, res = 160)
op <- par(mar = c(1, 1, 3.4, 1))
plot(gll, col = pal[as.character(cls)],
     border = ifelse(cls == "Not significant", "grey55", "grey15"),
     lwd = ifelse(cls == "Not significant", 0.3, 0.9),
     main = paste0("LISA clusters - poverty rate (Y), Java 2025\n",
                   "augmented Queen W, conditional permutation (", NSIM_M3, "), FDR q = 0.05"))
for (r in seq_len(nrow(AUG_LINKS))) {
  a <- idx(AUG_LINKS[r, 1]); b <- idx(AUG_LINKS[r, 2])
  lines(rbind(cn[a, ], cn[b, ]), col = "black", lwd = 1.8, lty = 2)
}
lab <- match(c(3528L, 3101L), java$bps)
text(cn[lab, 1], cn[lab, 2] + 0.15, labels = java$name[lab], cex = 0.7, font = 2, pos = 3, offset = 0.2)
legend("bottomleft", bty = "n", cex = 0.9, title = "LISA class (FDR q = 0.05)",
       legend = c(paste0(levels(cls),
                  c("  (Pamekasan)", "  (Jabodetabek metro core)", "", "", ""))),
       fill = unname(pal), border = "grey35")
legend("bottom", bty = "n", cex = 0.8, lty = 2, lwd = 1.8, col = "black",
       legend = "documented augmented link (Suramadu; Kep. Seribu - Jkt Utara)")
par(op); dev.off()
log("wrote ", basename(MAP_PNG))

###############################################################################
# RELEASE ADDITION: two checks quoted in the paper. They run after every other
# computation, so they cannot change any result above.
###############################################################################
## (a) permutation resolution: the six-seed probe above, repeated at 9,999
rule("PERMUTATION RESOLUTION - FDR count on Y, QueenAug, six seeds")
probe_counts <- function(nsim) sapply(probe_seeds, function(sd) {
  lmv <- localmoran_perm(java$Y, W[[PRIMARY]], nsim = nsim, zero.policy = TRUE,
                         alternative = "two.sided", iseed = sd)
  pc <- grep("^Pr\\(z .* E\\(Ii\\)\\) Sim$", colnames(lmv), value = TRUE)[1]
  pr <- lmv[, pc]; okp <- !is.na(pr)
  pf <- rep(NA_real_, length(pr)); pf[okp] <- p.adjust(pr[okp], "BH")
  c(uncorrected = sum(pr < Q, na.rm = TRUE), FDR_q05 = sum(pf < Q, na.rm = TRUE))
})
stab9 <- probe_counts(9999L)
colnames(stab9) <- paste0("seed_", probe_seeds)
log("nsim = 9,999:")
logcap(t(stab9))
log("FDR count range: 9,999 perm ", min(stab9["FDR_q05", ]), "-", max(stab9["FDR_q05", ]),
    " | 99,999 perm ", min(stab["FDR_q05", ]), "-", max(stab["FDR_q05", ]))

## (b) bridge decomposition: Moran's I of Y with each added link on its own
rule("ADDED LINKS - global Moran's I of Y")
with_links <- function(rows) {
  nb <- queen_nb
  for (r in rows) nb <- add_sym_link(nb, idx(AUG_LINKS[r, 1]), idx(AUG_LINKS[r, 2]))
  nb <- sym.attr.nb(nb); attr(nb, "region.id") <- attr(queen_nb, "region.id"); nb
}
set.seed(SEED)
mi_links <- sapply(list(Queen = integer(0), Queen_Suramadu = 1L, Queen_Seribu = 2L,
                        QueenAug = 1:2), function(rows) {
  lw <- nb2listw(with_links(rows), style = "W", zero.policy = TRUE)
  unname(moran.mc(java$Y, lw, nsim = 99, zero.policy = TRUE)$statistic)
})
logcap(round(mi_links, 4))
log("drop Queen -> QueenAug ", round(mi_links["Queen"] - mi_links["QueenAug"], 4),
    " | Suramadu alone ", round(mi_links["Queen"] - mi_links["Queen_Suramadu"], 4),
    " | QueenAug undirected edges ", sum(card(queen_aug_nb)) / 2)

close(log_con)
cat("\nDONE. Log: ", OUT_LOG, "\n")
