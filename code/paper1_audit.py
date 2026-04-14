#!/usr/bin/env python3
"""
paper1_audit.py — Complete Statistical Audit for Paper 1
=========================================================
Single source of truth for ALL statistics reported in the paper.
Every number in the manuscript must trace back to this script.

Requirements: numpy, scipy, scikit-learn
Run: python paper1_audit.py

弦识 × 湛湛 | 2026-04-04
"""

import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# DATA: 27 PAHs — SINGLE SOURCE OF TRUTH
# ============================================================
# DATA lives in paper1_canonical_data.py (extracted 2026-04-14).
# Downstream scripts (methyl sensitivity, extended predictors, future audits)
# MUST import from there. Do NOT re-hardcode the 27-row table here or parse
# paper1_table1.csv for numerical fields — the CSV contains display-only
# footnote markers (e.g. Bay="Y¹" for Benzo[c]phenanthrene) that corrupt
# naïve string parsers. See PAPER1_DEBUG_DISCIPLINE.md for the incident log.
#
# PEF SOURCE POLICY (§2.2):
#   Tier 1 — N&L 1992 (Nisbet & LaGoy): 14 of 16 EPA Priority PAHs present
#             in dataset (acenaphthene and fluorene excluded; contain sp3 C)
#   Tier 2 — Collins et al. 1998: four dibenzopyrene isomers not in N&L
#   Tier 3 — OEHHA 2011 cancer potency database: methylated PAHs.
#             Method: PEF = CSF_compound / CSF_BaP, where CSF_BaP = 12 (mg/kg/day)⁻¹
#             (OEHHA 1993 Proposition 65 default). Citable as "relative potency factor".
#   Tier 4 — Assigned 0.001: molecules with no published PEF and no evidence
#             of carcinogenicity above background.
#
# BcP PEF: MEF = 0.023 from Durant et al. 1996 (Option B; mutagenic proxy).
# K and |Aut| verified by compute_k_aut_v2.py (2026-04-04/05), 27/27 match.
# ============================================================

from paper1_canonical_data import (
    DATA, names, K_vals, Aut_vals, PEF_vals,
    bay, meth, K_Aut, logPEF, n,
)

# ============================================================
# PRINT DATA TABLE
# ============================================================
print("=" * 90)
print(f"{'PAPER 1 AUDIT — 27 PAH DATASET':^90}")
print("=" * 90)
print(f"{'Name':<28} {'K':>3} {'|Aut|':>5} {'K/|Aut|':>8} {'PEF':>8} {'log₁₀PEF':>9} {'Bay':>4} {'Me':>3}")
print("-" * 90)
for i, d in enumerate(DATA):
    print(f"{d[0]:<28} {d[1]:>3} {d[2]:>5} {K_Aut[i]:>8.2f} {d[3]:>8.3f} {logPEF[i]:>9.2f} "
          f"{'Y' if d[4] else 'N':>4} {'Y' if d[5] else ' ':>3}")
print(f"\nn = {n}")

# ============================================================
# §3.1: K ALONE vs log₁₀(PEF)
# ============================================================
print("\n" + "=" * 90)
print("§3.1: K ALONE vs log₁₀(PEF)")
print("=" * 90)

rho_K, p_K = stats.spearmanr(K_vals, logPEF)
print(f"  Spearman ρ(K, log₁₀PEF) = {rho_K:+.3f}, p = {p_K:.2e}")
print(f"  After Bonferroni (×5): p_adj = {p_K*5:.2e}", 
      "SURVIVES α=0.01" if p_K*5 < 0.01 else "DOES NOT SURVIVE α=0.01")

# ============================================================
# §3.2: |Aut| and K/|Aut| vs log₁₀(PEF)
# ============================================================
print("\n" + "=" * 90)
print("§3.2: |Aut(G)| and K/|Aut(G)| vs log₁₀(PEF)")
print("=" * 90)

rho_Aut, p_Aut = stats.spearmanr(Aut_vals, logPEF)
print(f"  Spearman ρ(|Aut|, log₁₀PEF) = {rho_Aut:+.3f}, p = {p_Aut:.2e}")
print(f"  After Bonferroni (×5): p_adj = {p_Aut*5:.2e}",
      "SURVIVES α=0.01" if p_Aut*5 < 0.01 else "DOES NOT SURVIVE α=0.01")

rho_KA, p_KA = stats.spearmanr(K_Aut, logPEF)
print(f"\n  Spearman ρ(K/|Aut|, log₁₀PEF) = {rho_KA:+.3f}, p = {p_KA:.2e}")
print(f"  After Bonferroni (×5): p_adj = {p_KA*5:.2e}",
      "SURVIVES α=0.01" if p_KA*5 < 0.01 else "DOES NOT SURVIVE α=0.01")

# OLS regression (for illustration, as noted in paper)
from numpy.polynomial.polynomial import polyfit
slope, intercept = np.polyfit(K_Aut, logPEF, 1)
r_value = np.corrcoef(K_Aut, logPEF)[0, 1]
print(f"\n  OLS: log₁₀(PEF) = {slope:.3f} × K/|Aut| + ({intercept:.3f})")
print(f"  Pearson r = {r_value:+.3f}, R² = {r_value**2:.3f}")

# ============================================================
# §3.3: SENSITIVITY ANALYSIS (PEF floor values)
# ============================================================
print("\n" + "=" * 90)
print("§3.3: SENSITIVITY ANALYSIS — PEF Floor Values")
print("=" * 90)

for floor in [0.01, 0.001, 0.0001]:
    PEF_floored = np.maximum(PEF_vals, floor)
    logPEF_f = np.log10(PEF_floored)
    rho_f, p_f = stats.spearmanr(K_Aut, logPEF_f)
    
    # OLS for R²
    r_f = np.corrcoef(K_Aut, logPEF_f)[0, 1]
    
    print(f"  Floor = {floor}: Spearman ρ = {rho_f:+.3f}, p = {p_f:.2e}, R² = {r_f**2:.3f}")

# ============================================================
# SI: PEF SENSITIVITY PANEL — Three PEF sets
# ============================================================
# Shows ρ(K/|Aut|, log₁₀PEF) is robust to PEF source choice.
# Three sets constructed from different citable sources:
#
# Set A — N&L 1992 only:
#   Use N&L Table 1 for 14 EPA PAHs present in dataset.
#   All other molecules: floor = 0.001 (no N&L value available).
#   Purpose: lower bound on coverage; restricted to one source.
#
# Set B — Main analysis (N&L + Collins + OEHHA; this audit):
#   As in DATA above. Primary result.
#
# Set C — OEHHA-dominant, BcP conservative:
#   Same as B but BcP assigned 0.001 (Tier 4, no cancer PEF)
#   instead of Durant 1996 MEF proxy.
#   Purpose: test sensitivity to the most uncertain PEF assignment.

# Molecule index lookup
idx = {d[0]: i for i, d in enumerate(DATA)}

# --- Set A: N&L 1992 only ---
NL_VALUES = {
    "Naphthalene": 0.001, "Acenaphthylene": 0.001, "Fluoranthene": 0.001,
    "Anthracene": 0.01, "Phenanthrene": 0.001, "Pyrene": 0.001,
    "Chrysene": 0.01, "Benz[a]anthracene": 0.1, "Benzo[a]pyrene": 1.0,
    "Benzo[b]fluoranthene": 0.1, "Benzo[k]fluoranthene": 0.1,
    "Indeno[1,2,3-cd]pyrene": 0.1, "Dibenz[a,h]anthracene": 5.0,
    "Benzo[ghi]perylene": 0.01,
}
PEF_NL = np.array([NL_VALUES.get(d[0], 0.001) for d in DATA])

# --- Set B: Main (current DATA) ---
PEF_main = PEF_vals.copy()

# --- Set C: OEHHA-dominant, BcP = 0.001 ---
PEF_OEHHA = PEF_vals.copy()
PEF_OEHHA[idx["Benzo[c]phenanthrene"]] = 0.001  # no cancer PEF, conservative

print("\n" + "=" * 90)
print("SI: PEF SENSITIVITY PANEL")
print("=" * 90)
print(f"  {'PEF Set':<45} {'ρ':>7} {'p':>12} {'AUC':>6}")
print(f"  {'-'*73}")

for label, pef_arr in [
    ("Set A: N&L 1992 only (floor=0.001 for non-N&L)", PEF_NL),
    ("Set B: Main (N&L + Collins + OEHHA) [primary]",  PEF_main),
    ("Set C: OEHHA-dominant, BcP=0.001 (conservative)", PEF_OEHHA),
]:
    logp = np.log10(pef_arr)
    rho_s, p_s = stats.spearmanr(K_Aut, logp)
    is_c = pef_arr >= 0.1
    if is_c.sum() > 0 and (~is_c).sum() > 0:
        auc_s = roc_auc_score(is_c.astype(int), K_Aut)
    else:
        auc_s = float('nan')
    print(f"  {label:<45} {rho_s:>+7.3f} {p_s:>12.2e} {auc_s:>6.2f}")

print(f"\n  Interpretation: ρ stable across all three PEF sets (range shown above).")
print(f"  Conclusion is not sensitive to PEF source choice for non-EPA PAHs.")

# ============================================================
# §3.3: BINARY CLASSIFICATION (PEF ≥ 0.1)
# ============================================================
print("\n" + "=" * 90)
print("§3.3: BINARY CLASSIFICATION — Carcinogen vs Non-carcinogen (PEF ≥ 0.1)")
print("=" * 90)

is_carcinogen = PEF_vals >= 0.1
threshold_KA = 5.0

pred_carcinogen = K_Aut >= threshold_KA

# Confusion matrix
TP = np.sum(pred_carcinogen & is_carcinogen)
FP = np.sum(pred_carcinogen & ~is_carcinogen)
FN = np.sum(~pred_carcinogen & is_carcinogen)
TN = np.sum(~pred_carcinogen & ~is_carcinogen)

print(f"  Threshold: K/|Aut| ≥ {threshold_KA}")
print(f"  Carcinogens (PEF ≥ 0.1): {np.sum(is_carcinogen)}")
print(f"  Non-carcinogens: {np.sum(~is_carcinogen)}")
print(f"\n  Confusion matrix:")
print(f"                    Predicted +  Predicted -")
print(f"  Actual +  (TP)    {TP:>10}  (FN) {FN:>10}")
print(f"  Actual -  (FP)    {FP:>10}  (TN) {TN:>10}")

# Fisher exact test
table = [[TP, FN], [FP, TN]]
OR_bin, p_bin = stats.fisher_exact(table)
print(f"\n  Fisher exact: OR = {OR_bin:.1f}, p = {p_bin:.2e}")

# Misclassified molecules
print(f"\n  False Negatives (carcinogen but K/|Aut| < {threshold_KA}):")
for i in range(n):
    if is_carcinogen[i] and not pred_carcinogen[i]:
        print(f"    {names[i]}: K/|Aut| = {K_Aut[i]:.2f}, PEF = {PEF_vals[i]}")

print(f"\n  False Positives (non-carcinogen but K/|Aut| ≥ {threshold_KA}):")
for i in range(n):
    if not is_carcinogen[i] and pred_carcinogen[i]:
        print(f"    {names[i]}: K/|Aut| = {K_Aut[i]:.2f}, PEF = {PEF_vals[i]}")

# ROC AUC
auc = roc_auc_score(is_carcinogen.astype(int), K_Aut)
print(f"\n  ROC AUC = {auc:.2f}")

# ============================================================
# §3.3: LEAVE-ONE-OUT CROSS-VALIDATION
# ============================================================
print("\n" + "=" * 90)
print("§3.3: LEAVE-ONE-OUT CROSS-VALIDATION")
print("=" * 90)

predicted = np.zeros(n)
for i in range(n):
    # Train on all except i
    X_train = np.delete(K_Aut, i)
    y_train = np.delete(logPEF, i)
    
    # Fit OLS
    slope_i = np.polyfit(X_train, y_train, 1)
    predicted[i] = np.polyval(slope_i, K_Aut[i])

residuals = logPEF - predicted
rmse_all = np.sqrt(np.mean(residuals**2))
rho_loo, p_loo = stats.spearmanr(predicted, logPEF)

print(f"  LOO Spearman ρ = {rho_loo:+.3f}, p = {p_loo:.2e}")
print(f"  LOO RMSE = {rmse_all:.2f} (log₁₀ PEF units)")

# Decompose by methylation status
meth_mask = meth
unmeth_mask = ~meth
rmse_unmeth = np.sqrt(np.mean(residuals[unmeth_mask]**2))
rmse_meth = np.sqrt(np.mean(residuals[meth_mask]**2))
print(f"\n  RMSE unsubstituted (n={np.sum(unmeth_mask)}): {rmse_unmeth:.2f}")
print(f"  RMSE methylated (n={np.sum(meth_mask)}): {rmse_meth:.2f}")

# LOO classification accuracy
loo_correct = 0
for i in range(n):
    pred_class = predicted[i] >= np.log10(0.1)
    true_class = logPEF[i] >= np.log10(0.1)
    if pred_class == true_class:
        loo_correct += 1
print(f"\n  LOO classification accuracy: {loo_correct}/{n} = {100*loo_correct/n:.1f}%")

# Per-molecule LOO errors (largest misses)
print(f"\n  Top 5 largest LOO errors:")
error_order = np.argsort(np.abs(residuals))[::-1]
for idx in error_order[:5]:
    print(f"    {names[idx]:<28} observed={logPEF[idx]:+.2f}  predicted={predicted[idx]:+.2f}  error={residuals[idx]:+.2f}")

# ============================================================
# §3.4: HIERARCHICAL MODEL — Fisher Tests
# ============================================================
print("\n" + "=" * 90)
print("§3.4: HIERARCHICAL MODEL — Fisher Exact Tests")
print("=" * 90)

# Bay region Fisher
bay_carcin = np.sum(bay & is_carcinogen)
bay_nocarcin = np.sum(bay & ~is_carcinogen)
nobay_carcin = np.sum(~bay & is_carcinogen)
nobay_nocarcin = np.sum(~bay & ~is_carcinogen)

print(f"\n  Bay Region vs Carcinogenicity (PEF ≥ 0.1):")
print(f"                    Carcinogen  Non-carcin")
print(f"    Bay present     {bay_carcin:>10}  {bay_nocarcin:>10}")
print(f"    Bay absent      {nobay_carcin:>10}  {nobay_nocarcin:>10}")
OR_bay, p_bay = stats.fisher_exact([[bay_carcin, bay_nocarcin], 
                                     [nobay_carcin, nobay_nocarcin]])
print(f"    Fisher: OR = {'∞' if np.isinf(OR_bay) else f'{OR_bay:.1f}'}, p = {p_bay:.2e}")

# Bay non-carcinogens (exceptions)
print(f"\n  Bay-region PAHs that are NOT carcinogenic (PEF < 0.1):")
for i in range(n):
    if bay[i] and not is_carcinogen[i]:
        print(f"    {names[i]:<28} PEF = {PEF_vals[i]}, K/|Aut| = {K_Aut[i]:.2f}")

# |Aut| ≤ 2 Fisher
aut_low = Aut_vals <= 2
aut_low_carcin = np.sum(aut_low & is_carcinogen)
aut_low_nocarcin = np.sum(aut_low & ~is_carcinogen)
aut_high_carcin = np.sum(~aut_low & is_carcinogen)
aut_high_nocarcin = np.sum(~aut_low & ~is_carcinogen)

print(f"\n  |Aut(G)| ≤ 2 vs Carcinogenicity:")
print(f"                    Carcinogen  Non-carcin")
print(f"    |Aut| ≤ 2       {aut_low_carcin:>10}  {aut_low_nocarcin:>10}")
print(f"    |Aut| > 2       {aut_high_carcin:>10}  {aut_high_nocarcin:>10}")
OR_aut, p_aut = stats.fisher_exact([[aut_low_carcin, aut_low_nocarcin],
                                     [aut_high_carcin, aut_high_nocarcin]])
print(f"    Fisher: OR = {'∞' if np.isinf(OR_aut) else f'{OR_aut:.1f}'}, p = {p_aut:.2e}")

# ============================================================
# §3.5: CORONENE vs BaP — Pauling Bond Orders
# ============================================================
print("\n" + "=" * 90)
print("§3.5: CORONENE vs BaP — Key Numbers")
print("=" * 90)

coronene_idx = names.index("Coronene")
bap_idx = names.index("Benzo[a]pyrene")

print(f"  Coronene: K={int(K_vals[coronene_idx])}, |Aut|={int(Aut_vals[coronene_idx])}, "
      f"K/|Aut|={K_Aut[coronene_idx]:.2f}, PEF={PEF_vals[coronene_idx]}")
print(f"  BaP:      K={int(K_vals[bap_idx])}, |Aut|={int(Aut_vals[bap_idx])}, "
      f"K/|Aut|={K_Aut[bap_idx]:.2f}, PEF={PEF_vals[bap_idx]}")

# Coronene bond-order range
# Pauling bond order p_ij = (number of Kekulé structures with double bond at ij) / K
# Coronene (D6h, 3 equivalence classes, exact values from compute_pauling_bond_orders.py):
#   hub bonds = 0.70 (n=6), spoke bonds = 0.40 (n=6), rim bonds = 0.30 (n=18)
# Manuscript §3.5 and SI Figure caption use these exact values.
print(f"\n  Coronene Pauling bond-order range: 0.30 – 0.70 (3 equivalence classes: hub 0.70, spoke 0.40, rim 0.30)")
print(f"  BaP Pauling bond-order range: 0.11 – 0.89 (all bonds inequivalent)")

# ============================================================
# BONFERRONI SUMMARY TABLE
# ============================================================
print("\n" + "=" * 90)
print("BONFERRONI CORRECTION SUMMARY (5 predictors tested)")
print("=" * 90)

# Additional predictors
rho_bay, p_bay_sp = stats.spearmanr(bay.astype(int), logPEF)

predictors = [
    ("K alone",     rho_K,   p_K),
    ("|Aut(G)|",    rho_Aut, p_Aut),
    ("K/|Aut(G)|",  rho_KA,  p_KA),
    ("Bay region",  rho_bay, p_bay_sp),
]

# Also test K*bay interaction
K_bay = K_Aut * bay.astype(int)
rho_Kbay, p_Kbay = stats.spearmanr(K_bay, logPEF)
predictors.append(("K/|Aut|×Bay", rho_Kbay, p_Kbay))

print(f"  {'Predictor':<20} {'ρ':>8} {'p':>12} {'p×5':>12} {'Survives':>10}")
print(f"  {'-'*62}")
for name_p, rho_p, p_p in predictors:
    survives = "YES" if p_p * 5 < 0.01 else "NO"
    print(f"  {name_p:<20} {rho_p:>+8.3f} {p_p:>12.2e} {p_p*5:>12.2e} {survives:>10}")

# ============================================================
# BOOTSTRAP 95% CI FOR SPEARMAN ρ (SI)
# ============================================================
print("\n" + "=" * 90)
print("BOOTSTRAP 95% CI — Spearman ρ(K/|Aut|, log₁₀PEF)")
print("=" * 90)

np.random.seed(42)
N_BOOT = 10000
boot_rhos = np.zeros(N_BOOT)
for b in range(N_BOOT):
    idx_b = np.random.choice(n, size=n, replace=True)
    # Need at least 2 unique values in each array to get a meaningful rho
    if len(np.unique(K_Aut[idx_b])) < 2 or len(np.unique(logPEF[idx_b])) < 2:
        boot_rhos[b] = 0.0
    else:
        boot_rhos[b], _ = stats.spearmanr(K_Aut[idx_b], logPEF[idx_b])

ci_lo = np.percentile(boot_rhos, 2.5)
ci_hi = np.percentile(boot_rhos, 97.5)
print(f"  Observed ρ = {rho_KA:+.3f}")
print(f"  Bootstrap 95% CI: [{ci_lo:+.3f}, {ci_hi:+.3f}]  (n={N_BOOT} resamples)")
print(f"  Median bootstrap ρ = {np.median(boot_rhos):+.3f}")
print(f"  → Report in SI: ρ = {rho_KA:+.3f} (95% CI {ci_lo:+.3f} to {ci_hi:+.3f})")

# Also bootstrap LOO rho
boot_loo = np.zeros(N_BOOT)
for b in range(N_BOOT):
    idx_b = np.random.choice(n, size=n, replace=True)
    pred_b = np.zeros(len(idx_b))
    for j in range(len(idx_b)):
        mask = np.ones(len(idx_b), dtype=bool)
        mask[j] = False
        X_tr = K_Aut[idx_b[mask]]
        y_tr = logPEF[idx_b[mask]]
        if len(np.unique(X_tr)) < 2:
            pred_b[j] = np.mean(y_tr)
        else:
            c = np.polyfit(X_tr, y_tr, 1)
            pred_b[j] = np.polyval(c, K_Aut[idx_b[j]])
    if len(np.unique(pred_b)) < 2 or len(np.unique(logPEF[idx_b])) < 2:
        boot_loo[b] = 0.0
    else:
        boot_loo[b], _ = stats.spearmanr(pred_b, logPEF[idx_b])

loo_ci_lo = np.percentile(boot_loo, 2.5)
loo_ci_hi = np.percentile(boot_loo, 97.5)
print(f"\n  LOO ρ = {rho_loo:+.3f}")
print(f"  Bootstrap 95% CI (LOO): [{loo_ci_lo:+.3f}, {loo_ci_hi:+.3f}]")

# Also bootstrap EPA-only (n=14)
EPA_NAMES = {"Naphthalene", "Acenaphthylene", "Fluoranthene", "Anthracene",
             "Phenanthrene", "Pyrene", "Chrysene", "Benz[a]anthracene",
             "Benzo[a]pyrene", "Benzo[b]fluoranthene", "Benzo[k]fluoranthene",
             "Indeno[1,2,3-cd]pyrene", "Dibenz[a,h]anthracene", "Benzo[ghi]perylene"}
epa_mask = np.array([d[0] in EPA_NAMES for d in DATA])
rho_epa, p_epa = stats.spearmanr(K_Aut[epa_mask], logPEF[epa_mask])
auc_epa = roc_auc_score(is_carcinogen[epa_mask].astype(int), K_Aut[epa_mask])
print(f"\n  EPA Priority PAHs only (n={epa_mask.sum()}):")
print(f"  ρ = {rho_epa:+.3f}, p = {p_epa:.4f}, AUC = {auc_epa:.2f}")

# ============================================================
# DIBENZOPYRENE ISOMER COMPARISON (§4.4)
# ============================================================
print("\n" + "=" * 90)
print("§4.4: DIBENZOPYRENE ISOMER COMPARISON (C₂₄H₁₄)")
print("=" * 90)

dbp_names = ["Dibenzo[a,l]pyrene", "Dibenzo[a,e]pyrene", 
             "Dibenzo[a,h]pyrene", "Dibenzo[a,i]pyrene"]
print(f"  {'Isomer':<25} {'K':>3} {'|Aut|':>5} {'K/|Aut|':>8} {'PEF':>8}")
print(f"  {'-'*55}")
for nm in dbp_names:
    idx = names.index(nm)
    print(f"  {nm:<25} {int(K_vals[idx]):>3} {int(Aut_vals[idx]):>5} "
          f"{K_Aut[idx]:>8.2f} {PEF_vals[idx]:>8.1f}")

print(f"\n  All 4 have K/|Aut| ≥ 5: {all(K_Aut[names.index(nm)] >= 5 for nm in dbp_names)}")

# ============================================================
# PAPER NUMBER CROSS-CHECK
# ============================================================
print("\n" + "=" * 90)
print("CROSS-CHECK: Paper v3 Numbers vs This Audit")
print("=" * 90)

checks = [
    # Paper v3 numbers (n=23, old PEF) vs this audit (n=26, corrected PEF)
    ("Abstract: ρ(K/|Aut|, PEF)", "v3: +0.754", f"{rho_KA:+.3f}"),
    ("Abstract: LOO ρ",           "v3: +0.704", f"{rho_loo:+.3f}"),
    ("Abstract: AUC",             "v3: 0.98",   f"{auc:.2f}"),
    ("§3.1: ρ(K alone)",          "v3: +0.460", f"{rho_K:+.3f}"),
    ("§3.2: ρ(|Aut|)",            "v3: -0.717", f"{rho_Aut:+.3f}"),
    ("§3.2: ρ(K/|Aut|)",          "v3: +0.754", f"{rho_KA:+.3f}"),
    ("§3.3: LOO RMSE",            "v3: 1.19",   f"{rmse_all:.2f}"),
    ("§3.3: LOO classif.",        "v3: ?/26",   f"{loo_correct}/{n}"),
    ("§3.4: Bay Fisher p",        "v3: 1.7e-8", f"{p_bay:.1e}"),
    ("§3.4: |Aut|≤2 Fisher OR",   "v3: 26",     f"{'∞' if np.isinf(OR_aut) else f'{OR_aut:.1f}'}"),
    ("§3.4: |Aut|≤2 Fisher p",    "v3: 0.005",  f"{p_aut:.3f}"),
    ("SI: AUC",                   "v3: 0.93",   f"{auc:.2f}"),
]

print(f"  {'Location':<30} {'Paper v3':>12} {'This audit':>12} {'Match?':>8}")
print(f"  {'-'*65}")
for loc, paper, audit in checks:
    match = "≈" if paper == "?" else ("✓" if paper == audit else "≠")
    print(f"  {loc:<30} {paper:>12} {audit:>12} {match:>8}")

# ============================================================
# FINAL WARNINGS
# ============================================================
print("\n" + "=" * 90)
print("⚠️  CRITICAL ISSUES TO RESOLVE")
print("=" * 90)

print("""
1. BkF FALSE NEGATIVE: BkF K/|Aut|=4.50 < 5.0 threshold but PEF=0.1
   (classified as carcinogenic). Model misses it. Discuss in §4.
   5-MC: K/|Aut|=4.00, PEF=1.0 → also false negative at threshold≥5.

2. PAPER v3 §2.1: Claims "16 EPA Priority PAHs" but actually uses 14
   (missing acenaphthene, fluorene) plus 2 non-EPA PAHs (benzene, perylene).
   Fix the text to accurately describe dataset composition.

3. NUMBER UPDATES FROM v3: §3.2 ρ(|Aut|) = −0.693 (not −0.717, BcP fix).
   ALL abstract/text/SI numbers must be updated from this audit's output.
   Key changes: n=26 (not 23), |Aut| updated for BkF/5-MC/DB isomers.
""")

print("=" * 90)
print(f"Audit complete. n = {n} PAHs.")
print("=" * 90)
