# pah-kekule-descriptor

Code and data for **Liu Z. (2026). *Graph-algebraic invariants predict PAH carcinogenic potency without machine learning.* Journal of Chemical Information and Modeling** (submitted).

**What it does.** On 27 polycyclic aromatic hydrocarbons (PAHs) with established cancer potency factors, the symmetry-discounted Kekulé ratio *K*/|Aut(*G*)| predicts log₁₀ PEF at Spearman ρ = +0.745 (Bonferroni-corrected, bootstrap 95% CI [+0.546, +0.855]). Adding the binary bay-region indicator to form *K*/|Aut(*G*)| × bay lifts Spearman ρ to +0.834 (Bonferroni p × 17 = 1.14 × 10⁻⁶). No training, no machine learning, no fitted parameters.

**Why it matters.** *K*/|Aut(*G*)| is computed in closed form from the π-electron graph of a PAH. Every step is auditable on paper; every number in the manuscript is reproduced by one of the five scripts in `code/`. This repository is the reference implementation cited by the Data Availability statement of the paper.

## Headline result

| Predictor                 | Spearman ρ | Bonferroni p (× 17) | BH q (× 17) |
|---------------------------|-----------:|--------------------:|------------:|
| *K*/|Aut(*G*)| × bay      | **+0.834** | 1.14 × 10⁻⁶         | 1.14 × 10⁻⁶ |
| bay-region indicator      | +0.802     | 1.68 × 10⁻⁵         | 8.4 × 10⁻⁶  |
| *K*/|Aut(*G*)|            | **+0.745** | 1.42 × 10⁻⁴         | 4.7 × 10⁻⁵  |
| MolecularWeight (comparator) | +0.690  | 1.5 × 10⁻³          | 4.3 × 10⁻⁴  |
| XLogP3 (comparator)       | +0.625     | 9.8 × 10⁻³          | 2.4 × 10⁻³  |
| Hosoya *Z*                | +0.490     | 0.235               | 0.063       |
| *K* alone                 | +0.420     | 0.495               | 0.115       |

The graph-algebraic composite wins by a margin (Δρ = +0.144 over the strongest physicochemical comparator, MW). *K* alone is **not** significant under Bonferroni × 17 — the symmetry quotient |Aut(*G*)| is doing the work. See `paper1_JCIM/extended_predictors_N17/README.md` and `paper1_extended_log.txt` for the full ranked table and methodological discussion.

## Reproduce in three commands

```bash
git clone https://github.com/mahpmice/pah-kekule-descriptor.git && cd pah-kekule-descriptor
python3 -m pip install -r requirements.txt
python3 code/paper1_audit.py
```

Full reproduction (audit + extended-predictor stress test + SI tables + figures + manuscript PDF) is documented in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

## What's in the box

- `code/` — eight canonical modules: `paper1_canonical_data.py` (27-PAH data SSOT), `compute_k_aut_v2.py` (K and |Aut|), `paper1_audit.py` (primary statistics), `paper1_extended_predictors.py` (N=17 Bonferroni×17), `paper1_methyl_sensitivity.py` (N=24 subset), `generate_paper1_si.py` (SI tables), `generate_paper1_figures.py` (figures), `build_paper1_pdf.py` (PDF build).
- `data/paper1_table1.csv` — 27 PAHs × (name, CID, SMILES, K, |Aut|, K/|Aut|, bay, log₁₀PEF). Every value has a source tier documented in `REPRODUCIBILITY.md` §5.
- `paper1_JCIM/manuscript_Liu2026_PAH_Kekule.pdf` — the submitted manuscript.
- `paper1_JCIM/extended_predictors_N17/` — post-hoc N=17 Bonferroni × 17 + BH q<0.05 stress test (17 predictors, including 10 graph-theoretic and 2 physicochemical comparators).

## Citation

> Liu Z. (2026). Graph-algebraic invariants predict PAH carcinogenic potency without machine learning. *Journal of Chemical Information and Modeling* (submitted). Code archive: https://github.com/mahpmice/pah-kekule-descriptor

## Patent notice

Methods described in this work are covered by **U.S. Provisional Patent Application No. 64/037,961** (filed 2026-04-13). Use of the code in this repository is governed by the LICENSE file at repository root; the license does not grant any patent rights beyond its own terms. For commercial licensing please contact the author.

## License

See [`LICENSE`](LICENSE). The repository is released under a permissive open-source license to support reviewer verification and scholarly re-use; patent licensing is handled separately.

## Contact

Zhiwei Liu (Independent Researcher) — mahpmiceliu@gmail.com — ORCID: 0009-0004-3926-0720
