# Reproducibility — Paper 1 (JCIM)

**Manuscript.** Liu Z. (2026). *Graph-algebraic invariants predict PAH carcinogenic potency without machine learning.* Journal of Chemical Information and Modeling (submitted).

**Repository.** `https://github.com/mahpmice/pah-kekule-descriptor`

**Patent notice.** Methods used in this work are covered by U.S. Provisional Patent Application No. 64/037,961 (filed 2026-04-13). The code in this repository is released under the LICENSE at repository root; the release does not grant any patent rights beyond the applicable license terms.

---

## 1. System requirements

| Component | Tested version | Minimum |
|-----------|----------------|---------|
| Python    | 3.11.9         | 3.9     |
| networkx  | 3.2.1          | 3.0     |
| numpy     | 1.26.4         | 1.24    |
| scipy     | 1.11.4         | 1.10    |
| scikit-learn | 1.4.0       | 1.3     |
| matplotlib | 3.8.2         | 3.7     |
| reportlab  | 4.0.9         | 4.0     |

Platforms tested: macOS 14 (Arial available) and Linux sandbox (DejaVu fallback built into `build_paper1_pdf.py`). No GPU required. Total wall time for a full rebuild <3 minutes on a 2023 laptop.

Network: `paper1_extended_predictors.py` requires outbound HTTPS to `https://pubchem.ncbi.nlm.nih.gov` (for XLogP3 and MolecularWeight on 27 CIDs; ~1 s per call). The other four scripts are offline.

---

## 2. One-command reproduction

From the repository root:

```bash
python3 -m pip install -r requirements.txt
python3 code/paper1_audit.py                 # Table 1 statistics, Bonferroni × 5, bootstrap CI, LOO
python3 code/paper1_extended_predictors.py   # N=17 extended predictors, Bonferroni × 17, BH q<0.05
python3 code/paper1_methyl_sensitivity.py    # N=24 methyl-excluded sensitivity (§3.3)
python3 code/generate_paper1_si.py           # SI Table S1–S5 (CSVs)
python3 code/generate_paper1_figures.py      # Figures 1–4 + Figure S1 (PDF + PNG)
python3 code/build_paper1_pdf.py             # manuscript PDF (14 pages)
```

Each script prints a log on stdout; written artifacts land in `data/`, `figures/`, and `paper1_JCIM/extended_predictors_N17/`.

---

## 3. Expected numerical outputs

The following values from `paper1_audit.py` must match the manuscript to the stated precision. Any deviation is a reproducibility failure — please file an issue.

| Statistic                              | Expected value | Location in manuscript |
|----------------------------------------|----------------|------------------------|
| Spearman ρ(K alone)                    | +0.420         | §3.1                   |
| Spearman ρ(\|Aut\|)                    | −0.688         | §3.2 Table 2           |
| Spearman ρ(K/\|Aut\|)                  | +0.745         | §3.2 Table 2           |
| Bootstrap 95% CI for ρ(K/\|Aut\|)     | [+0.546, +0.855] | §3.2                 |
| Leave-one-out Spearman ρ               | +0.708         | §3.3                   |
| LOO 95% CI                             | [+0.480, +0.834] | §3.3                 |
| Spearman ρ(K/\|Aut\| × bay)            | +0.834         | §3.2 Table 2           |
| Spearman ρ(bay region)                 | +0.802         | §3.2 Table 2           |
| EPA Priority (n=14) ρ                  | +0.759         | §3.5                   |
| EPA Priority AUC (threshold PEF≥1)     | 0.95           | §3.5                   |

From `paper1_extended_predictors.py` (N=17 audit, Bonferroni × 17 + BH q<0.05):

| Statistic                              | Expected value | Location |
|----------------------------------------|----------------|----------|
| K/\|Aut\| × bay, Bonferroni p (×17)    | 1.14 × 10⁻⁶    | §3.2 extended table |
| K/\|Aut\|, Bonferroni p (×17)          | 1.42 × 10⁻⁴    | §3.2 extended table |
| K alone, Bonferroni p (×17)            | 0.495          | §3.1 ("K alone does not survive") |
| Predictors surviving Bonferroni × 17   | 9 / 17         | §3.2                 |
| Predictors surviving BH q<0.05         | 16 / 17        | §3.2                 |

From `paper1_methyl_sensitivity.py` (N=24 methyl-excluded subset):

| Statistic                               | Expected value | Location |
|-----------------------------------------|----------------|----------|
| N=24 Spearman ρ(K/\|Aut\|)              | +0.775         | §3.3     |
| N=24 Spearman ρ(K/\|Aut\| × bay)        | +0.868         | §3.3     |
| Δρ for K/\|Aut\| (N=24 − N=27)         | +0.031         | §3.3     |
| Δρ for K/\|Aut\| × bay (N=24 − N=27)   | +0.034         | §3.3     |
| AUC for K/\|Aut\| at PEF ≥ 1 (N=24)     | 0.87           | §3.3     |

---

## 4. Directory layout

```
.
├── README.md
├── REPRODUCIBILITY.md           (this file)
├── LICENSE
├── requirements.txt
├── code/
│   ├── compute_k_aut_v2.py          Source of truth for K and |Aut| on 27 PAH π-graphs
│   ├── paper1_audit.py              Primary statistical audit (Table 1, Bonferroni × 5, bootstrap CI, LOO)
│   ├── paper1_extended_predictors.py  N=17 extended-predictor audit (Bonferroni × 17, BH q<0.05)
│   ├── paper1_methyl_sensitivity.py   N=24 methyl-excluded sensitivity (response to §2.3)
│   ├── generate_paper1_si.py        Supporting Information tables S1–S5
│   ├── generate_paper1_figures.py   Figures 1–4 + Figure S1 (PDF + PNG)
│   └── build_paper1_pdf.py          Manuscript PDF build
├── data/
│   ├── paper1_table1.csv            27 PAH × (name, CID, SMILES, K, |Aut|, K/|Aut|, bay, log10PEF)
│   └── paper1_si_table_s1..s5.csv   SI data tables
├── figures/
│   ├── fig1_main_results_v4.(pdf|png)
│   ├── fig2_sensitivity_loo_v4.(pdf|png)
│   ├── fig3_hierarchy_v4.(pdf|png)
│   ├── fig4_coronene_bap_v4.(pdf|png)
│   └── figS1_ROC_v4.(pdf|png)
└── paper1_JCIM/
    ├── manuscript_Liu2026_PAH_Kekule.pdf   14-page main manuscript (includes Notes on patent)
    ├── manuscript_Liu2026_PAH_Kekule.md
    ├── supporting_information.pdf / .md
    ├── provenance_v6.md                    Versioned audit trail of every number
    └── extended_predictors_N17/
        ├── paper1_extended_results.csv         17 × (predictor, ρ, raw p, Bonferroni p, BH q, survival flags)
        ├── paper1_extended_per_pah.csv         27 × (name, CID, log10PEF, 17 descriptor values)
        ├── paper1_extended_log.txt             ranked summary + K-alone / Hosoya-Z diagnostic
        ├── paper1_methyl_sensitivity.csv       5 predictors × (N=27 stats, N=24 stats, Δρ)
        ├── paper1_methyl_sensitivity_log.txt   methyl-excluded interpretation + §3.3 template
        ├── README.md                           English scope + reproduce block
        └── README_中文解读.md                  中文详解 + manuscript edit plan
```

---

## 5. Data provenance

`data/paper1_table1.csv` is the canonical v4 (v6.1 after 2026-04-14 CID cleaning) dataset. Every value has a source tier documented in `paper1_JCIM/provenance_v6.md`:

- **Tier 1** — Nisbet & LaGoy 1992 PEFs for 14 of 16 EPA Priority PAHs (acenaphthene and fluorene excluded: sp³ carbon in the ring system).
- **Tier 2** — Collins et al. 1998 for four dibenzopyrene isomers.
- **Tier 3** — OEHHA 2011 cancer potency database for methylated PAHs, converted to PEF via `PEF = CSF_compound / CSF_BaP` with `CSF_BaP = 12 (mg/kg/day)⁻¹` (OEHHA 1993 Proposition 65 default).
- **Tier 4** — 0.001 assigned to molecules with no published PEF and no evidence of carcinogenicity above background.

`K` and `|Aut|` values in `data/paper1_table1.csv` are independently verified by `code/compute_k_aut_v2.py` on the SMILES strings embedded in that script; all 27/27 match.

---

## 6. Verification checklist for reviewers

```bash
# 1. Fresh clone
git clone https://github.com/mahpmice/pah-kekule-descriptor.git
cd pah-kekule-descriptor

# 2. Isolated environment
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 3. Reproduce headline numbers
python3 code/paper1_audit.py | tee audit.log
grep "K/|Aut|  " audit.log              # expect ρ = +0.745, Bonferroni ✓
grep "K/|Aut|.*bay" audit.log           # expect ρ = +0.834
grep "LOO"  audit.log                   # expect LOO ρ = +0.708

# 4. Reproduce the N=17 extended audit
python3 code/paper1_extended_predictors.py
head -5 paper1_JCIM/extended_predictors_N17/paper1_extended_results.csv
```

Exact agreement with §3.1–§3.3 of the manuscript constitutes successful reproduction.

---

## 7. Archival

A tagged release at submission will be snapshotted to Zenodo for a citable DOI (`https://doi.org/10.5281/zenodo.XXXXXXX` — placeholder until the Zenodo hook is completed). The Data Availability statement in the manuscript points to both the GitHub URL and the Zenodo DOI once the snapshot exists.

---

*Every number in Paper 1 must trace back to this repository. If a scripted run disagrees with the manuscript, the run is right and the manuscript is wrong — please file an issue and cite the line of `paper1_JCIM/provenance_v6.md` that the failing number came from.*
