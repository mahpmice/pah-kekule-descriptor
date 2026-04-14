#!/usr/bin/env python3
"""
Generate Supporting Information tables for Paper 1.
"""

from __future__ import annotations

import csv
import contextlib
import importlib.util
from math import log10, sqrt
from pathlib import Path
import io

import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parent


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


audit = load_module(ROOT / "paper1_audit.py", "paper1_audit_mod")
compute = load_module(ROOT / "compute_k_aut_v2.py", "paper1_compute_mod")


DATA = audit.DATA
names = [d[0] for d in DATA]
K_vals = np.array([d[1] for d in DATA], dtype=float)
Aut_vals = np.array([d[2] for d in DATA], dtype=float)
PEF_vals = np.array([d[3] for d in DATA], dtype=float)
bay = np.array([d[4] for d in DATA])
meth = np.array([d[5] for d in DATA])
K_Aut = K_vals / Aut_vals
logPEF = np.log10(PEF_vals)


CID_MAP = {
    "Benzene": 241,
    "Naphthalene": 931,
    "Acenaphthylene": 9161,
    "Fluoranthene": 9154,
    "Anthracene": 8418,
    "Phenanthrene": 995,
    "Pyrene": 31423,
    "Triphenylene": 9170,
    "Chrysene": 9171,
    "Benz[a]anthracene": 5954,
    "Benzo[c]phenanthrene": 9136,
    "Benzo[a]pyrene": 2336,
    "3-Methylcholanthrene": 1674,
    "5-Methylchrysene": 19427,
    "DMBA": 6001,
    "Dibenz[a,h]anthracene": 5889,
    "Dibenzo[a,l]pyrene": 9119,
    "Benzo[b]fluoranthene": 9153,
    "Benzo[k]fluoranthene": 9158,
    "Perylene": 9142,
    "Benzo[ghi]perylene": 9117,
    "Coronene": 9115,
    "Indeno[1,2,3-cd]pyrene": 9131,
    "Dibenzo[a,e]pyrene": 9126,
    "Dibenzo[a,h]pyrene": 9108,
    "Dibenzo[a,i]pyrene": 9106,
    "Benzo[e]pyrene": 9128,
}


PUBCHEM_CONFIRMED = {
    "5-Methylchrysene": {
        "cid": 19427,
        "title": "5-Methylchrysene",
        "smiles": "Cc1cc2ccccc2c2ccc3ccccc3c12",
        "note": "PubChem 2D record reproduces n_pi=18, K=8, |Aut|=2.",
    },
    "DMBA": {
        "cid": 6001,
        "title": "7,12-Dimethylbenz(a)anthracene",
        "smiles": "Cc1c2ccccc2c(C)c2c1ccc1ccccc12",
        "note": "PubChem 2D record reproduces n_pi=18, K=7, |Aut|=1.",
    },
    "Benzo[k]fluoranthene": {
        "cid": 9158,
        "title": "Benzo[k]fluoranthene",
        "smiles": "c1ccc2cc3c(cc2c1)-c1cccc2cccc-3c12",
        "note": "pynauty cross-check reproduces |Aut|=2.",
    },
    "Dibenzo[a,e]pyrene": {
        "cid": 9126,
        "title": "DIBENZO(a,e)PYRENE",
        "smiles": "c1ccc2c(c1)cc1c3ccccc3c3cccc4ccc2c1c43",
        "note": "PubChem correction: CID 9126, not 5765; reproduces K=17, |Aut|=1.",
    },
    "Dibenzo[e,l]pyrene": {
        "cid": 9122,
        "title": "Dibenzo(e,l)pyrene",
        "smiles": "c1ccc2c(c1)c1cccc3c4ccccc4c4cccc2c4c13",
        "note": "PubChem qualitative isomer outside main dataset; gives n_pi=24, K=20, |Aut|=4, K/|Aut|=5.0.",
    },
}


def write_csv(path: Path, header, rows):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def compute_loo():
    predicted = np.zeros(len(K_Aut))
    for i in range(len(K_Aut)):
        x_train = np.delete(K_Aut, i)
        y_train = np.delete(logPEF, i)
        slope = np.polyfit(x_train, y_train, 1)
        predicted[i] = np.polyval(slope, K_Aut[i])
    residuals = logPEF - predicted
    return predicted, residuals


def main():
    # Table S1
    s1_rows = []
    for name, meta in PUBCHEM_CONFIRMED.items():
        if name == "Dibenzo[e,l]pyrene":
            G = compute.extract_pi_system(compute.parse_smiles_to_graph(meta["smiles"]))
            K = compute.count_perfect_matchings(G)
            Aut = compute.count_automorphisms(G)
            n_pi = G.number_of_nodes()
            ratio = K / Aut
        else:
            idx = names.index(name)
            K = int(K_vals[idx])
            Aut = int(Aut_vals[idx])
            n_pi = int(
                compute.extract_pi_system(compute.parse_smiles_to_graph(meta["smiles"])).number_of_nodes()
            )
            ratio = K_Aut[idx]
        s1_rows.append(
            [name, meta["cid"], meta["title"], meta["smiles"], n_pi, K, Aut, f"{ratio:.2f}", meta["note"]]
        )
    write_csv(
        ROOT / "paper1_si_table_s1.csv",
        ["Name", "CID", "PubChem_title", "Canonical_SMILES", "n_pi", "K", "Aut", "K_over_Aut", "Verification_note"],
        s1_rows,
    )

    # Table S2
    s2_rows = []
    for i, d in enumerate(DATA, start=1):
        smiles = next(sm for nm, sm, *_ in compute.MOLECULES if nm == d[0])
        n_pi = compute.extract_pi_system(compute.parse_smiles_to_graph(smiles)).number_of_nodes()
        s2_rows.append(
            [
                i,
                d[0],
                CID_MAP[d[0]],
                n_pi,
                int(K_vals[i - 1]),
                int(Aut_vals[i - 1]),
                f"{K_Aut[i - 1]:.2f}",
                f"{PEF_vals[i - 1]:.3f}",
                f"{logPEF[i - 1]:+.2f}",
                "Y" if bay[i - 1] else "N",
                "Y" if meth[i - 1] else "N",
                d[6],
            ]
        )
    write_csv(
        ROOT / "paper1_si_table_s2.csv",
        ["#", "Name", "CID", "n_pi", "K", "Aut", "K_over_Aut", "PEF", "log10PEF", "Bay", "Methylated", "PEF_source"],
        s2_rows,
    )

    # Table S3
    floors = [0.01, 0.001, 0.0001]
    s3_rows = []
    for floor in floors:
        logp = np.log10(np.maximum(PEF_vals, floor))
        rho, p = stats.spearmanr(K_Aut, logp)
        s3_rows.append([f"Floor={floor}", f"{rho:+.3f}", f"{p:.2e}", "PEF floor sensitivity"])
    p_nl = np.array([audit.NL_VALUES.get(d[0], 0.001) for d in DATA], dtype=float)
    p_main = PEF_vals.copy()
    p_cons = PEF_vals.copy()
    p_cons[names.index("Benzo[c]phenanthrene")] = 0.001
    for label, arr in [
        ("Set A: N&L only", p_nl),
        ("Set B: Main hierarchy", p_main),
        ("Set C: BcP conservative", p_cons),
    ]:
        rho, p = stats.spearmanr(K_Aut, np.log10(arr))
        auc = roc_auc_score((arr >= 0.1).astype(int), K_Aut)
        s3_rows.append([label, f"{rho:+.3f}", f"{p:.2e}", f"AUC={auc:.2f}"])
    write_csv(ROOT / "paper1_si_table_s3.csv", ["Analysis", "rho", "p", "Notes"], s3_rows)

    # Table S4
    predicted, residuals = compute_loo()
    s4_rows = []
    for i, name in enumerate(names):
        true_class = "carcinogen" if PEF_vals[i] >= 0.1 else "non-carcinogen"
        pred_class = "carcinogen" if predicted[i] >= log10(0.1) else "non-carcinogen"
        s4_rows.append(
            [
                name,
                f"{K_Aut[i]:.2f}",
                f"{logPEF[i]:+.2f}",
                f"{predicted[i]:+.2f}",
                f"{residuals[i]:+.2f}",
                "Y" if meth[i] else "N",
                true_class,
                pred_class,
            ]
        )
    write_csv(
        ROOT / "paper1_si_table_s4.csv",
        ["Name", "K_over_Aut", "Observed_log10PEF", "LOO_predicted_log10PEF", "Residual", "Methylated", "True_class", "Predicted_class"],
        s4_rows,
    )

    # Table S5
    dbp = [
        ("Dibenzo[a,l]pyrene", 9119, None),
        ("Dibenzo[a,e]pyrene", 9126, PUBCHEM_CONFIRMED["Dibenzo[a,e]pyrene"]["smiles"]),
        ("Dibenzo[a,h]pyrene", 9108, None),
        ("Dibenzo[a,i]pyrene", 9106, None),
        ("Dibenzo[e,l]pyrene", 9122, PUBCHEM_CONFIRMED["Dibenzo[e,l]pyrene"]["smiles"]),
    ]
    s5_rows = []
    for name, cid, smiles in dbp:
        if name in names:
            i = names.index(name)
            pef = f"{PEF_vals[i]:.1f}"
            status = "In primary dataset"
            kval = int(K_vals[i])
            aval = int(Aut_vals[i])
            ratio = K_Aut[i]
        else:
            G = compute.extract_pi_system(compute.parse_smiles_to_graph(smiles))
            kval = compute.count_perfect_matchings(G)
            aval = compute.count_automorphisms(G)
            ratio = kval / aval
            pef = "not assigned"
            status = "Qualitative comparator only"
        s5_rows.append([name, cid, kval, aval, f"{ratio:.2f}", pef, status])
    write_csv(ROOT / "paper1_si_table_s5.csv", ["Isomer", "CID", "K", "Aut", "K_over_Aut", "PEF", "Status"], s5_rows)

    # Supporting Information markdown
    md = f"""# Supporting Information

## For

Graph Symmetry and Kekule Structure Localization Predict CYP1-Mediated Bioactivation and Carcinogenic Potency of Polycyclic Aromatic Hydrocarbons

Zhiwei Liu

### Contents

- Table S1. PubChem confirmations and graph-verification notes for corrected or high-risk structures.
- Table S2. Complete 26-PAH descriptor table used in the main analysis.
- Table S3. Sensitivity analyses for floor assignment and benzo[c]phenanthrene proxy assignment.
- Table S4. Leave-one-out predictions and residuals for every molecule.
- Table S5. Dibenzopyrene isomer comparison, including the qualitative comparator dibenzo[e,l]pyrene.
- Figure S1. ROC curve already exported as `figS1_ROC.pdf` and `figS1_ROC.png`.

### Key SI Notes

- 5-Methylchrysene is confirmed against PubChem CID 19427.
- DMBA is confirmed against PubChem CID 6001.
- Benzo[k]fluoranthene remains a real boundary false negative after `pynauty` confirmation of `|Aut| = 2`.
- Dibenzo[a,e]pyrene is corrected to PubChem CID 9126 with `K = 17` and `|Aut| = 1`.
- Dibenzo[e,l]pyrene maps to PubChem CID 9122 and gives `K = 20`, `|Aut| = 4`, `K/|Aut| = 5.0`.

### Exported Files

- [Table S1](paper1_si_table_s1.csv)
- [Table S2](paper1_si_table_s2.csv)
- [Table S3](paper1_si_table_s3.csv)
- [Table S4](paper1_si_table_s4.csv)
- [Table S5](paper1_si_table_s5.csv)
- [Figure S1 PDF](figS1_ROC.pdf)
- [Figure S1 PNG](figS1_ROC.png)

### LOO Summary

- LOO Spearman rho = +0.708 (95% bootstrap CI: +0.480 to +0.834)
- LOO RMSE = 1.25 log10 units
- RMSE unsubstituted (n=24) = 1.15
- RMSE methylated (n=3) = 1.88
- LOO classification accuracy = 19/27 (70.4%)

### Binary Summary

- Threshold `K/|Aut| >= 5.0`
- TP = 11, TN = 11, FP = 2, FN = 2
- AUC = 0.95 (full dataset); AUC = 0.95 (EPA Priority PAHs only, n=14)

### Bootstrap 95% CI

- Primary ρ = +0.745 (95% CI: +0.546 to +0.855, n=10000 resamples, seed=42)
- EPA-only (n=14): ρ = +0.759, p = 0.0016, AUC = 0.95

"""
    (ROOT / "paper1_supporting_information.md").write_text(md)


if __name__ == "__main__":
    main()
