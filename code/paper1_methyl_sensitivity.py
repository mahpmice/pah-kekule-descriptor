#!/usr/bin/env python3
"""
paper1_methyl_sensitivity.py — Methylated-PAH exclusion sensitivity for Paper 1
===============================================================================
Response to the methylated-PAH limitation discussed in §4.3 of the manuscript.

§2.3 of the manuscript excludes exocyclic methyl substituents from the π-graph.
A3.6 requires an honest accounting: if the three methylated PAHs in the dataset
(3-Methylcholanthrene, 5-Methylchrysene, and DMBA) are dropped entirely,
does K/|Aut(G)| still predict log10 PEF?

Three outcomes are informative regardless of direction:

  (i)  ρ on N=24 >> +0.745 → methylated PAHs are "noise"; paper's scope is
       rigorously "unsubstituted PAHs" and should say so.
  (ii) ρ on N=24 ≈ +0.745 → methyl exclusion is a faithful abstraction.
  (iii) ρ on N=24 << +0.745 → methylated PAHs are *contributors* to the
        correlation, which tempers the mechanistic claim.

This script computes (ρ, raw p, Bonferroni p × 17, AUC against PEF ≥ 1.0) for
the five pre-registered predictors on the N=24 subset, and writes results +
discussion hooks to `paper1_JCIM/extended_predictors_N17/`.

Data: `data/paper1_table1.csv` (canonical v4, CID-cleaned v6.1 2026-04-14).

Requirements: numpy, scipy, scikit-learn
Run from repo root:
    python3 code/paper1_methyl_sensitivity.py

Zhiwei Liu | 2026-04-14 | methyl-excluded subset sensitivity
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_CSV = REPO_ROOT / "data" / "paper1_table1.csv"
OUT_DIR = REPO_ROOT / "paper1_JCIM" / "extended_predictors_N17"
OUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_CSV = OUT_DIR / "paper1_methyl_sensitivity.csv"
LOG_TXT = OUT_DIR / "paper1_methyl_sensitivity_log.txt"

# The three methylated PAHs excluded per §2.3 discussion
METHYL_NAMES = {"3-Methylcholanthrene", "5-Methylchrysene", "DMBA"}

# Bonferroni factor — keep consistent with the N=17 audit so reviewers
# reading §3.2 don't see a different correction on sensitivity tables.
BONFERRONI_FACTOR = 17

# Five pre-registered predictors (order matches manuscript §3.2)
#
# IMPORTANT: the CSV's K_over_Aut column is rounded to 2 decimals for human
# readability (e.g. Benzene = 0.17 instead of 2/12 = 0.16666...). Using the
# rounded column for Spearman ρ collapses ties that shouldn't exist and drifts
# ρ(K/|Aut|×Bay) from +0.834 to +0.812. Always recompute K/|Aut| from the
# integer-valued K and Aut columns for statistics.
def _bay(row: dict) -> float:
    # CSV uses "Y" or "Y¹" (with a footnote superscript for Benzo[c]phenanthrene's
    # fjord-region classification, which behaves as a bay region for the
    # potency-correlation purposes of this paper). Strict "== Y" would drop
    # BcP's bay=1 and drift ρ(K/|Aut|×Bay) from +0.834 to +0.812. Match a
    # leading "Y" (any footnote suffix allowed) to stay consistent with the
    # hard-coded DATA in paper1_audit.py.
    return 1.0 if row["Bay"].strip().startswith("Y") else 0.0


def _k_over_aut(row: dict) -> float:
    return float(row["K"]) / float(row["Aut"])


PREDICTORS = [
    ("K", lambda row: float(row["K"])),
    ("|Aut|", lambda row: float(row["Aut"])),
    ("K/|Aut|", _k_over_aut),
    ("Bay", _bay),
    ("K/|Aut| x Bay", lambda row: _k_over_aut(row) * _bay(row)),
]


def _parse_log10_pef(cell: str) -> float:
    """CSV uses Unicode minus sign '−' for negatives; coerce to float."""
    s = cell.strip().replace("−", "-").replace("\u2212", "-")
    # handle leading "+"
    if s.startswith("+"):
        s = s[1:]
    return float(s)


def load_table() -> list[dict]:
    rows: list[dict] = []
    with DATA_CSV.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            if not r.get("Name"):
                continue
            rows.append(r)
    if len(rows) != 27:
        sys.exit(f"[FATAL] paper1_table1.csv must contain 27 PAHs; got {len(rows)}")
    return rows


def _compute_predictor_stats(values: np.ndarray, pef: np.ndarray) -> dict:
    rho, raw_p = stats.spearmanr(values, pef)
    bonf_p = min(1.0, raw_p * BONFERRONI_FACTOR)
    # AUC at threshold PEF >= 1.0 (log10 PEF >= 0)
    pef_binary = (pef >= 0.0).astype(int)
    if pef_binary.sum() == 0 or pef_binary.sum() == len(pef_binary):
        auc = float("nan")
    else:
        auc = roc_auc_score(pef_binary, values)
    return {
        "rho": float(rho),
        "raw_p": float(raw_p),
        "bonferroni_p": float(bonf_p),
        "auc_pef_ge_1": float(auc),
    }


def main() -> None:
    rows = load_table()
    names_all = [r["Name"] for r in rows]
    pef_all = np.array([_parse_log10_pef(r["log10PEF"]) for r in rows])

    # Full N=27 and N=24 (methyl-excluded)
    keep_mask = np.array([n not in METHYL_NAMES for n in names_all])
    n_full = len(rows)
    n_sub = int(keep_mask.sum())
    assert n_sub == 24, f"expected 24 non-methyl PAHs, got {n_sub}"

    # Sanity: verify the three excluded names are actually in the table
    excluded_names = [n for n in names_all if n in METHYL_NAMES]
    if set(excluded_names) != METHYL_NAMES:
        sys.exit(f"[FATAL] expected {METHYL_NAMES}, found {excluded_names}")

    print(f"N_full = {n_full}; N_methyl_excluded = {n_sub}")
    print(f"Excluded: {sorted(excluded_names)}")
    print()

    # Collect per-predictor statistics on both N=27 and N=24
    result_rows = []
    log_lines = [
        "paper1_methyl_sensitivity — N=24 (methyl-excluded) sensitivity audit",
        "=====================================================================",
        f"Excluded: {', '.join(sorted(excluded_names))}",
        f"Bonferroni factor: x{BONFERRONI_FACTOR} (consistent with N=17 audit)",
        "",
        f"{'predictor':<20s}{'N=27 ρ':>10s}{'N=27 Bonf p':>14s}{'N=27 AUC':>10s}"
        f"{'N=24 ρ':>10s}{'N=24 Bonf p':>14s}{'N=24 AUC':>10s}{'Δρ':>10s}",
        "-" * 98,
    ]

    for name, getter in PREDICTORS:
        vals_all = np.array([getter(r) for r in rows])
        vals_sub = vals_all[keep_mask]
        pef_sub = pef_all[keep_mask]

        full_stats = _compute_predictor_stats(vals_all, pef_all)
        sub_stats = _compute_predictor_stats(vals_sub, pef_sub)
        delta = sub_stats["rho"] - full_stats["rho"]

        result_rows.append({
            "predictor": name,
            "n_full": n_full,
            "rho_full": full_stats["rho"],
            "bonf_p_full": full_stats["bonferroni_p"],
            "auc_full": full_stats["auc_pef_ge_1"],
            "n_sub": n_sub,
            "rho_sub": sub_stats["rho"],
            "bonf_p_sub": sub_stats["bonferroni_p"],
            "auc_sub": sub_stats["auc_pef_ge_1"],
            "delta_rho": delta,
        })
        log_lines.append(
            f"{name:<20s}{full_stats['rho']:>+10.3f}{full_stats['bonferroni_p']:>14.3e}"
            f"{full_stats['auc_pef_ge_1']:>10.3f}"
            f"{sub_stats['rho']:>+10.3f}{sub_stats['bonferroni_p']:>14.3e}"
            f"{sub_stats['auc_pef_ge_1']:>10.3f}{delta:>+10.3f}"
        )

    # Write CSV
    with RESULTS_CSV.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(result_rows[0].keys()))
        w.writeheader()
        for row in result_rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})

    # Interpretation block
    k_over_aut = next(r for r in result_rows if r["predictor"] == "K/|Aut|")
    rho_full = k_over_aut["rho_full"]
    rho_sub = k_over_aut["rho_sub"]
    delta = k_over_aut["delta_rho"]

    log_lines += [
        "",
        "Interpretation (K/|Aut| as reference predictor)",
        "-------------------------------------------------",
        f"N=27 Spearman ρ = {rho_full:+.3f}",
        f"N=24 Spearman ρ = {rho_sub:+.3f}",
        f"Δρ (N=24 − N=27) = {delta:+.3f}",
        "",
    ]

    if delta > 0.05:
        verdict = (
            "Δρ > +0.05: removing methylated PAHs strengthens the correlation. "
            "The manuscript's scope is honestly restricted to unsubstituted PAHs; "
            "§4.3 should state this as a scope clarification rather than a limitation."
        )
    elif delta < -0.05:
        verdict = (
            "Δρ < −0.05: the methylated PAHs contribute positively to the correlation; "
            "they are NOT dilutive noise. §4.3 must acknowledge that the aggregate ρ "
            "depends in part on the inclusion of DMBA/3-MC/5-MC, and §2.3 cannot "
            "frame methyl exclusion as loss-free."
        )
    else:
        verdict = (
            "|Δρ| ≤ 0.05: the K/|Aut| correlation is robust to methyl PAH inclusion/"
            "exclusion. §4.3 can cite this sensitivity as evidence that the π-graph "
            "abstraction is not carried by the three methylated outliers."
        )

    log_lines.append(verdict)
    log_lines.append("")
    log_lines.append(
        "Reviewer-safe reporting template (§3.3 insertion):"
    )
    log_lines.append(
        f'"Excluding the three methylated PAHs (DMBA, 3-methylcholanthrene, '
        f'5-methylchrysene), the N=24 subset yields Spearman ρ = {rho_sub:+.3f} '
        f'for K/|Aut(G)| (Bonferroni p × 17 = {k_over_aut["bonf_p_sub"]:.3e}, '
        f'AUC at PEF ≥ 1 = {k_over_aut["auc_sub"]:.2f}). '
        f'Δρ = {delta:+.3f} relative to the full 27-PAH dataset."'
    )

    log_text = "\n".join(log_lines) + "\n"
    LOG_TXT.write_text(log_text, encoding="utf-8")

    # Also echo to stdout
    print(log_text)
    print(f"[written] {RESULTS_CSV}")
    print(f"[written] {LOG_TXT}")


if __name__ == "__main__":
    main()
