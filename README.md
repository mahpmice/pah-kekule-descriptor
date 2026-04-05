# pah-kekule-descriptor

**Graph Symmetry and Kekulé Structure Localization Predict CYP1-Mediated Bioactivation and Carcinogenic Potency of Polycyclic Aromatic Hydrocarbons**

Zhiwei Liu — Independent Researcher

---

## Overview

This repository contains all code and data for the paper cited above. The central claim is that the ratio *K*/|Aut(*G*)| — Kekulé structure count divided by graph automorphism group order — predicts PAH carcinogenic potency with Spearman ρ = +0.745 (*p* = 8.38 × 10⁻⁶) across 27 PAHs, while *K* alone does not survive Bonferroni correction.

Every number in the paper traces back to `code/statistical_audit.py`. Run that script to reproduce all reported statistics.

---

## Repository Structure

```
pah-kekule-descriptor/
├── README.md                         # this file
├── LICENSE                           # MIT
├── requirements.txt                  # Python dependencies
│
├── data/
│   ├── pah_dataset.csv               # 27-PAH dataset: SMILES, K, |Aut|, PEF, classification
│   └── pah_descriptors.csv           # descriptor table with provenance annotations
│
└── code/
    ├── compute_k_aut.py              # Kekulé count K and |Aut(G)| for all 26+1 PAHs
    │                                 #   — recursive perfect-matching enumeration for K
    │                                 #   — NetworkX VF2 algorithm for |Aut(G)|
    │                                 #   — validates against Cyvin & Gutman (1988) values
    │
    ├── compute_pauling_bond_orders.py # Exact Pauling bond orders for coronene and BaP
    │                                  #   — p(u,v) = K(G-{u,v}) / K(G)
    │                                  #   — verifies mechanistic claim in §3.5
    │
    ├── statistical_audit.py          # Single source of truth for all paper statistics
    │                                 #   — Spearman correlations, Bonferroni correction
    │                                 #   — Fisher exact tests, ROC AUC, LOO-CV
    │                                 #   — Bootstrap 95% CI for Spearman ρ
    │                                 #   — EPA-only validation (n=14)
    │
    └── generate_figures.py           # Reproduces Figures 1–4 and Figure S1
                                      #   — requires: matplotlib, rdkit, scipy, sklearn
```

---

## Reproducing the Results

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Verify K and |Aut(G)| values

```bash
python code/compute_k_aut.py
```

Expected output: all 27 PAHs with K and |Aut(G)|, verified against literature values. No mismatches should be reported.

### 3. Reproduce all paper statistics

```bash
python code/statistical_audit.py
```

This reproduces every number in the paper in order of appearance (§3.1 through §4.4), plus the Bonferroni correction summary table and bootstrap confidence intervals.

### 4. Compute Pauling bond orders (§3.5, §4.1)

```bash
python code/compute_pauling_bond_orders.py
```

Reproduces the coronene/BaP bond-order analysis cited in §3.5:
- Coronene: three equivalence classes at *p* = 0.30, 0.40, 0.70 (range 0.40)
- BaP: continuous distribution 0.11–0.89 (range 0.78); maximum at bay-region bond

### 5. Regenerate figures

```bash
python code/generate_figures.py
```

Outputs PNG and PDF versions of Figures 1–4 and Figure S1 in the current directory.

---

## Data

`data/pah_dataset.csv` — 27 PAHs with the following columns:

| Column | Description |
|--------|-------------|
| `name` | Compound name (IUPAC or common) |
| `pubchem_cid` | PubChem Compound ID (SMILES source) |
| `smiles` | Canonical SMILES from PubChem |
| `n_pi` | π-system atom count (exocyclic methyls excluded) |
| `K` | Kekulé structure count (perfect matchings) |
| `Aut` | Graph automorphism group order \|Aut(*G*)\| |
| `K_over_Aut` | Symmetry-corrected localization index |
| `PEF` | Potency Equivalency Factor relative to BaP (PEF = 1.0) |
| `log10_PEF` | log₁₀(PEF) — regression endpoint |
| `bay_fjord` | Bay or fjord region present (Y/N) |
| `methylated` | Contains exocyclic methyl group(s) (Y/N) |
| `PEF_source` | Citable source for PEF value |
| `classification_5_0` | Binary classification at threshold K/\|Aut\| ≥ 5.0 (TP/TN/FP/FN) |

**PEF source hierarchy:**
1. Nisbet & LaGoy (1992) — 14 EPA Priority PAHs
2. Collins et al. (1998) — four dibenzopyrene isomers
3. OEHHA (2011) cancer slope factors / BaP slope factor — methylated PAHs
4. Tier 4 assigned (0.001) — no published PEF, no evidence of carcinogenicity

**Note on benzo[c]phenanthrene (CID 9136):** No cancer PEF exists. The value 0.023 is a mutagenic equivalency factor (MEF) from Durant et al. (1996), used as a proxy and flagged explicitly in sensitivity analyses.

---

## Requirements

```
networkx>=3.0
scipy>=1.9
numpy>=1.23
scikit-learn>=1.1
matplotlib>=3.6
rdkit>=2022.09
```

Python ≥ 3.9 required.

---

## Key Results (from `statistical_audit.py`)

| Descriptor | Spearman ρ | *p* (raw) | Bonferroni-adjusted *p* | Survives α = 0.01 |
|------------|------------|-----------|------------------------|-------------------|
| *K* alone | +0.420 | 0.0291 | 0.146 | No |
| \|Aut(*G*)\| | −0.688 | 7.38 × 10⁻⁵ | 3.69 × 10⁻⁴ | Yes |
| *K*/\|Aut(*G*)\| | **+0.745** | **8.38 × 10⁻⁶** | **4.19 × 10⁻⁵** | **Yes** |
| Bay/fjord region | +0.802 | 4.83 × 10⁻⁷ | 2.42 × 10⁻⁶ | Yes |
| *K*/\|Aut\| × Bay | +0.834 | 6.68 × 10⁻⁸ | 3.34 × 10⁻⁷ | Yes |

Bootstrap 95% CI for *K*/|Aut(*G*)|: ρ = +0.745 (95% CI +0.546 to +0.855, *n* = 10 000 resamples).

Binary classification at threshold *K*/|Aut(*G*)| ≥ 5.0: AUC = 0.95, Fisher OR = 33.0 (*p* = 4.23 × 10⁻⁴).

---

## Citation

> Liu, Z. Graph Symmetry and Kekulé Structure Localization Predict CYP1-Mediated Bioactivation and Carcinogenic Potency of Polycyclic Aromatic Hydrocarbons. *J. Chem. Inf. Model.* **2026** (submitted).

---

## License

MIT License. See `LICENSE` for details.

---

## Correspondence

Zhiwei Liu — liuzhiwei [at] [institution]
