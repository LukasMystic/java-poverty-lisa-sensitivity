# How Neighbour Definition and Permutation Resolution Shape Poverty Cluster Detection in Java, Indonesia

Code, data and results for the paper of the same title by Anita Rahayu and
Stanley Pratama Teguh (Bina Nusantara University), accepted for the 1st
International Conference on Computing and Digital Technology (ICoDiT 2026),
proceedings in Springer Communications in Computer and Information Science (CCIS).

The paper holds the data, the statistic and the multiple-testing correction fixed
and varies only the spatial weights and the permutation count. Under five
weighting schemes the number of FDR-significant local Moran clusters of the
2025 poverty rate across Java's 119 regencies and cities ranges from 4 to 16. At
9,999 permutations the corrected count depends on the random seed (4 to 6 across
six seeds); at 99,999 permutations it is 5 under every seed.

## Reproduce the results

```sh
Rscript spatial_M1_M3.R
```

The script builds the 119-unit dataset from the raw files, checks it against a
validation anchor, constructs the spatial weights, and computes every statistic
behind Tables 1 to 3. It takes about three minutes on an Apple M-series laptop.
All outputs are committed, and a run overwrites them in place, so `git diff`
afterwards shows whether anything changed. On the authors' machine a clean run
reproduces every result CSV byte for byte.

| Output | Content |
|---|---|
| `java_merged.gpkg`, `java_merged.csv` | the 119-unit dataset with geometry |
| `M1_M3_results.txt` | full run log: data checks, weights, all statistics, seed probes |
| `M2_global.csv` | global Moran's I and Geary's C, every variable and scheme |
| `M3_lisaY_all_schemes.csv` | local Moran and Gi* for the poverty rate, every unit and scheme |
| `M3_lisa_predictors_primary.csv` | local Moran for the three covariates under Queen+ |
| `M3_correlogram.csv` | Moran's I by neighbour lag order |
| `OLS_baseline.csv` | baseline regression, not used in the paper |
| `M3_LISA_map_augmentedQueen.png` | quick-look cluster map |

## Where each result in the paper comes from

| In the paper | Source |
|---|---|
| Table 1, global Moran's I | `M2_global.csv` |
| Table 2, significant local clusters by scheme | `M3_lisaY_all_schemes.csv` (columns `sig_raw`, `sig_fdr`, `G_sig_fdr`); counts also in `M1_M3_results.txt` |
| Table 3, seed stability at 99,999 permutations | `M1_M3_results.txt`, block "FDR-count seed stability" |
| Corrected count of 4 to 6 at 9,999 permutations | `M1_M3_results.txt`, block "PERMUTATION RESOLUTION" |
| Moran's I change from the two added links (0.484 to 0.406, Suramadu 0.069) | `M1_M3_results.txt`, block "ADDED LINKS" |
| Gresik and Sumenep kNN neighbours | `M1_M3_results.txt`, section M1 |
| Fig. 2 data | `M3_lisaY_all_schemes.csv`, scheme `QueenAug` |
| Fig. 3 data | `M3_correlogram.csv` |

## Spatial weights

`weights/` holds the five matrices as spdep GAL files (`<scheme>.gal`) and as
edge lists with the row-standardised weight of each directed link
(`<scheme>_edges.csv`). `weights/weights_summary.csv` gives the link counts and
the global Moran's I of the poverty rate under each one, which matches
`M2_global.csv`. They were exported from the same construction as the script.

## Settings

| Setting | Value |
|---|---|
| Units | 119 regencies and cities of Java (BPS codes 3101 to 3674) |
| Primary weights (Queen+) | Queen contiguity plus two links: Bangkalan (3526) to Kota Surabaya (3578), the Suramadu Bridge; Kepulauan Seribu (3101) to Kota Jakarta Utara (3175). 263 undirected links (526 directed), mean degree 4.42, one component |
| Comparison weights | bare Queen; kNN with k = 4, 5, 6 from `st_point_on_surface()` in an azimuthal equidistant projection (lat 0 = -7.5, lon 0 = 110), not symmetrised |
| Standardisation | row-standardised (`style = "W"`) for every scheme |
| Global statistics | `moran.mc`, `geary.mc`, 9,999 permutations, one-sided (`alternative = "greater"`), `set.seed(20260826)` |
| Local statistics | `localmoran_perm` and `localG_perm` (Gi*, `include.self`), two-sided conditional permutation, 99,999 permutations, `iseed = 20260826` |
| Multiple testing | Benjamini-Hochberg (`p.adjust(method = "BH")`), q = 0.05, separately per scheme and statistic over defined units (118 for local Moran under bare Queen, where Kepulauan Seribu is an isolate) |
| Stability probe | local Moran on Y under Queen+, seeds 1, 20260826, 42, 2024, 99999, 7, at 99,999 and at 9,999 permutations |
| Correlogram | `sp.correlogram`, Queen+, lags 1 to 6, analytical variance under randomisation |

## Requirements

R 4.6.1 with sf 1.1.2, spdep 1.4.2 and readxl 1.5.0 (GDAL 3.8.5, GEOS 3.13.0,
PROJ 9.5.1, as bundled with the sf binary). The environment of the reference
run is in `sessionInfo.txt`. Other versions should give the same statistics,
but permutation results are only guaranteed identical with the same spdep
version and seeds.

## Data

| File | Content | Source and licence |
|---|---|---|
| `Percentage of Population Living Below the National Poverty Line by Regency_Municipality, 2025.csv` | poverty rate, Y | BPS-Statistics Indonesia, [bps.go.id](https://www.bps.go.id) |
| `IPM/*.csv` | Human Development Index, X2, one file per province | BPS-Statistics Indonesia |
| `TPT+TPAK/*.csv` | open unemployment rate X1 and labour force participation X3, August Sakernas round | BPS-Statistics Indonesia |
| `dataset.xlsx` | an earlier extract of the 38 East Java units, used only as a validation gate: the rebuilt data must match it exactly | derived from BPS statistics |
| `idn_admin2_java.gpkg` | the 122 Java features of layer `idn_admin2` from the Indonesia COD-AB geodatabase, version 01 (valid on 2020-04-01), all fields, unmodified | OCHA / HDX, [Indonesia - Subnational Administrative Boundaries](https://data.humdata.org/dataset/cod-ab-idn), source BPS, licensed [CC BY-IGO 3.0](http://creativecommons.org/licenses/by/3.0/igo/legalcode) |

The BPS tables are redistributed with attribution to BPS-Statistics Indonesia.

**Boundary file.** The full geodatabase is 277 MB, too large for this
repository, so a Java-only extract is included. It was made with

```sh
ogr2ogr -f GPKG idn_admin2_java.gpkg idn_admin_boundaries.gdb idn_admin2 \
  -nln idn_admin2 -where "adm1_pcode IN ('ID31','ID32','ID33','ID34','ID35','ID36')"
```

The script reads this extract when `idn_admin_boundaries.gdb` is absent. To
start from the original file instead, download `idn_admin_boundaries.gdb.zip`
from the HDX page above and unzip it here; reading it needs a GDAL 3.9 or newer
`ogr2ogr` on `PATH`, because GDAL 3.8.5 cannot read one of its field types. The
script detects this and converts the layer automatically. Both routes give the
same 119 units and identical results. Three non-administrative polygons (codes
ending in 88 or 99: Hutan, Waduk Cirata, Wadung Kedungombo) are dropped by the
script and reported in its log.

## Licence

Code: MIT (see `LICENSE`). Data files keep the terms of their sources listed above.

## Citation

See `CITATION.cff`. Please cite the paper and the archived version of this
repository.
