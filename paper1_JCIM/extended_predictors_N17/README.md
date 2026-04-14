# Paper 1 (JCIM) — N=17 Extended-Predictor Audit

**Purpose.** Stress-test the five pre-registered graph-algebraic invariants of the manuscript (K, |Aut(G)|, K/|Aut(G)|, bay-region indicator, K/|Aut(G)|×bay) against twelve additional descriptors under Bonferroni (×17) and Benjamini–Hochberg (q<0.05) corrections on the 27-PAH dataset. Provides a post-hoc multiplicity audit across an expanded 17-descriptor panel, declared explicitly in the manuscript §2.4.

**Endpoint.** log₁₀ PEF (Collins 1998 + Durant 1996 + OEHHA, `data/paper1_table1.csv`).

## Reproduce

From the repository root:

```bash
python3 -m pip install -r requirements.txt
python3 code/paper1_extended_predictors.py
```

Outputs are written to this directory (`paper1_JCIM/extended_predictors_N17/`):

| file | description |
|------|-------------|
| `paper1_extended_per_pah.csv`  | 27 PAH × (name, CID, log₁₀PEF, 17 descriptor values) |
| `paper1_extended_results.csv`  | 17 predictors × (ρ, raw p, Bonferroni p, BH q, survival flags) |
| `paper1_extended_log.txt`      | Ranked summary with K-alone / Hosoya-Z diagnostic |

Network requirement: live HTTPS access to `https://pubchem.ncbi.nlm.nih.gov` for XLogP3 and MolecularWeight on each CID (Group B, 2 descriptors). The other 15 predictors need no network.

## Predictor set

Seventeen predictors total.

| # | Predictor | Group | Computed on |
|---|-----------|-------|-------------|
| 1 | K | pre-reg | π-graph |
| 2 | \|Aut\| | pre-reg | π-graph |
| 3 | K/\|Aut\| | pre-reg | π-graph |
| 4 | Bay (0/1) | pre-reg | molecule |
| 5 | K/\|Aut\| × Bay | pre-reg | π-graph × indicator |
| 6 | N_vertices | A: graph-theoretic | π-graph |
| 7 | N_edges | A | π-graph |
| 8 | N_rings | A | π-graph cyclomatic |
| 9 | Randić ¹χ | A | π-graph |
| 10 | Wiener W | A | π-graph distances |
| 11 | Hosoya Z | A | π-graph matchings |
| 12 | Balaban J | A | π-graph |
| 13 | Zagreb M1 | A | π-graph degrees |
| 14 | Zagreb M2 | A | π-graph edges |
| 15 | λ₁ | A | π-graph adjacency |
| 16 | XLogP3 | B: physicochem | whole molecule |
| 17 | MolecularWeight | B | whole molecule |

π-graph = ring-atom subgraph after methyl-carbon exclusion (see `code/compute_k_aut_v2.py::extract_pi_system`). SMILES are reused verbatim from `compute_k_aut_v2.py` (K/|Aut| audit source of truth).

## Headline result

K/|Aut(G)| × Bay remains rank-1 among all 17 predictors at Spearman ρ = +0.834 (Bonferroni p × 17 = 1.14 × 10⁻⁶, BH q = 1.14 × 10⁻⁶). K/|Aut(G)| alone is rank-3 at ρ = +0.745 (Bonferroni p × 17 = 1.4 × 10⁻⁴). K alone drops out under Bonferroni × 17 (p = 0.495) — consistent with and strengthening the §3.1 argument that K is necessary but not sufficient without the symmetry quotient. Hosoya Z (ρ = +0.490) and whole-molecule comparators MW (ρ = +0.690) and XLogP3 (ρ = +0.625) remain below K/|Aut(G)|, supporting the mechanistic-legibility positioning of §4.1 rather than a generic-descriptor scaling argument.

See `paper1_extended_log.txt` for the full ranked table and ranked diagnostic.

## Reference

This audit is logged in `paper1_JCIM/provenance_v6.md` as version v6.1 → v6.2 extension.
