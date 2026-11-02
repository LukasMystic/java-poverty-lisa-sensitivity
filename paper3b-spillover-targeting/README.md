# Agglomeration Shadows and the Limits of Spillover Aware Poverty Targeting in Java, Indonesia

Code and results for the second paper of this project, by Anita Rahayu and
Stanley Pratama Teguh (Bina Nusantara University), accepted for the 1st
International Conference on Computing and Digital Technology (ICoDiT 2026).

It shares its dataset and its spatial weights with the paper in the repository
root (neighbour definition and permutation resolution in LISA cluster
detection). Both are built from the same 119-unit Java dataset and the same
augmented Queen contiguity graph; this folder adds the spatial econometric
model, the impact decomposition and the targeting exercise.

The paper estimates an SDEM on the 2025 poverty rate, decomposes its impacts,
repeats the estimation across eight weighting schemes, builds a priority
ranking from the estimated spillover, and tests how far that ranking is driven
by the geometry of the neighbour graph rather than by measured poverty.

## Reproduce the results

```sh
Rscript ../spatial_M1_M3.R    # at the repository root: builds java_merged.gpkg
Rscript run_3b.R              # here
```

`run_3b.R` copies the two inputs it needs from the root and runs the four
scripts in order. The copies are git-ignored. Everything uses seed 20260826 and
every output is committed, so a re-run overwrites the files in place and
`git diff` shows whether anything moved.

| Script | What it does |
|---|---|
| `spatial_M4.R` | first pass: LM and robust LM tests, the six candidate models, impacts. Ends on SDM by AIC |
| `spatial_M4b.R` | re-examination. Writes `M4b_decision_rule.txt` **before** fitting anything, then applies it and lands on SDEM despite its 0.5-point AIC disadvantage |
| `spatial_M4c.R` | re-estimates the neighbour effect across eight well-posed graphs, reports the ill-posed ones separately, rebuilds the two figures |
| `spatial_M5_M6.R` | the typology, and the spillover-aware ranking with its sensitivity analysis |

The typology outputs (`M5_typology.csv`, `M5_crosstab.csv`, `fig_M5_typology.png`) are
produced by the last script but are not used in either paper; they are included so a run of
`run_3b.R` can be compared against committed files in full.

`M4b_decision_rule.txt` is the selection rule as it was written, before the
final comparison. It is included because the paper's argument about specification
choice rests on the rule having been fixed in advance.

## Where each result in the paper comes from

| In the paper | Source |
|---|---|
| Table 1, candidate model comparison | `M4b_model_comparison.csv` |
| Table 2, impact decomposition | `M4b_impacts_with_CI.csv` |
| Table 3, neighbour effects and 95% CIs across eight graphs | `M4c_theta_robustness.csv` |
| The two excluded ill-posed graphs | `M4c_illposed_graphs.csv`, and `M4c_results.txt` |
| Table 4, ranking sensitivity across the confidence interval | `M6_sensitivity.csv` |
| Table 5, mean rank change by neighbour count | `M6_rankings.csv`, columns `k` and `rank_change_R0_minus_R1` |
| Rank correlation of −0.753 between neighbour count and rank change | `M6_rankings.csv`, same two columns |
| Specification tests (LM, robust LM) | `M4_specification_tests.csv` |
| LR tests and the Hausman test | `M4b_results.txt` |
| Fig. 1, model comparison | `fig_M4b_model_comparison.png` |
| Fig. 2, impact decomposition with 95% CIs | `fig_M4b_impacts.png` |

`M6_rankings.csv` carries one row per regency: `k` is the number of
administrative neighbours, `own_gain` and `nbr_loss` are the two components of
the model-implied score, both scaled by poverty rates, and
`rank_change_R0_minus_R1` is the move between the poverty ranking and the
spillover-aware one.

## Settings

Seed 20260826. Response is log poverty rate. Scope is the 119 regencies and
cities of Java. The primary weights matrix is the augmented Queen graph built
by the root script, row-standardised. Package versions and the R version are in
`../sessionInfo.txt`.

## Licence

Code: MIT, as in the repository root. Data files keep the terms of their
sources, listed in the root README.
