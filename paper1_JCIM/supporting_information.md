# Supporting Information

## For

Graph Symmetry and Kekulé Structure Localization Correlate with Carcinogenic Potency of Polycyclic Aromatic Hydrocarbons

Zhiwei Liu

### Contents

- Table S1. PubChem confirmations and graph-verification notes for corrected or high-risk structures (extended commentary on five cases requiring explicit annotation).
- Table S2. Complete 27-PAH descriptor table used in the main analysis. The `CID` column lists the canonical PubChem identifier for every entry; all 27 structures were programmatically retrieved from PubChem and the resulting molecular graphs were independently verified before being passed to the K and |Aut(G)| computations.
- Table S3. Sensitivity analyses for floor assignment and benzo[c]phenanthrene proxy assignment.
- Table S4. Leave-one-out predictions and residuals for every molecule. The three methylated PAHs (3-methylcholanthrene, 5-methylchrysene, and 7,12-dimethylbenz[a]anthracene) — the systematic failure mode discussed in the main text — appear with LOO residuals of +1.30, +1.78, and +2.41 log₁₀ units respectively, all in the direction of underprediction; this is the basis for the decomposition RMSE_unsubstituted = 1.15 vs RMSE_methylated = 1.88.
- Table S5. Dibenzopyrene isomer comparison, including the qualitative comparator dibenzo[e,l]pyrene.
- Figure S1. ROC curve for the binary classification task at threshold K/|Aut(G)| ≥ 5.0 across the full 27-PAH dataset (full-set AUC = 0.95; EPA Priority subset n=14 AUC = 0.95). Exported as `figS1_ROC.pdf` and `figS1_ROC.png`.

### Selection of the binary threshold K/|Aut(G)| ≥ 5.0

The threshold was selected at the natural gap in the descriptor distribution (Table S2): dibenz[a,h]anthracene (K/|Aut| = 5.00, PEF = 5.0) is the lowest-K/|Aut| compound classified as a regulatory carcinogen in the primary hierarchy and serves as the anchor point. Shifting the cutoff by ±0.5 around 5.0 trades a single boundary case in each direction without improving the overall confusion matrix; the value 5.0 was therefore retained as the simplest defensible cut. The two false positives (benzo[ghi]perylene, benzo[e]pyrene) and two false negatives (benzo[k]fluoranthene, dibenz[a,h]anthracene) are explicitly retained and discussed in the main text rather than tuned away.

### Key SI Notes

- All 27 structures are cross-verified against PubChem entries (Table S2, `CID` column). The five cases below required explicit annotation because they involved nomenclature ambiguity, prior literature errors, or boundary failure modes; the remaining 22 structures matched their canonical PubChem entries on the first retrieval.
- 5-Methylchrysene is confirmed against PubChem CID 19427.
- 7,12-Dimethylbenz[a]anthracene (DMBA) is confirmed against PubChem CID 6001.
- Benzo[k]fluoranthene remains a real boundary false negative after `pynauty` confirmation of `|Aut| = 2`.
- Dibenzo[a,e]pyrene is corrected to PubChem CID 9126 with `K = 17` and `|Aut| = 1`.
- Dibenzo[e,l]pyrene maps to PubChem CID 9122 and gives `K = 20`, `|Aut| = 4`, `K/|Aut| = 5.0`.

### Exported Files

- [Table S1](SI_table_S1_pubchem_confirmations.csv)
- [Table S2](SI_table_S2_full_descriptors.csv)
- [Table S3](SI_table_S3_sensitivity.csv)
- [Table S4](SI_table_S4_LOO_predictions.csv)
- [Table S5](SI_table_S5_dibenzopyrenes.csv)
- [Figure S1 PDF](SI_figS1_ROC.pdf)
- [Figure S1 PNG](SI_figS1_ROC_preview.png)

### LOO Summary

- LOO Spearman rho = +0.708 (95% bootstrap CI: +0.480 to +0.834)
- LOO RMSE = 1.25 log10 units
- RMSE unsubstituted (n=24) = 1.15
- RMSE methylated (n=3) = 1.88
- LOO classification accuracy = 19/27 (70.4%)

### Binary Summary

- Threshold `K/|Aut| >= 5.0`
- TP = 11, TN = 12, FP = 2, FN = 2 (sum = 27)
- AUC = 0.95 (full dataset); AUC = 0.95 (EPA Priority PAHs only, n=14)

### Bootstrap 95% CI

- Primary ρ = +0.745 (95% CI: +0.546 to +0.855, n=10000 resamples, seed=42)
- EPA-only (n=14): ρ = +0.759, p = 0.0016, AUC = 0.95

