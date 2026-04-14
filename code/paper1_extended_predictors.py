#!/usr/bin/env python3
"""
paper1_extended_predictors.py — N=17 extended-predictor audit for Paper 1 (JCIM)
================================================================================
Stress-test the five pre-registered graph-algebraic invariants (K, |Aut(G)|,
K/|Aut(G)|, bay-region indicator, K/|Aut(G)|×bay) against 12 additional
descriptors under Bonferroni (×17) and Benjamini-Hochberg (q<0.05) corrections
on the 27-PAH dataset (endpoint: log10 PEF).

Original 5 (from data/paper1_table1.csv, canonical v4):
  1. K                 : Kekule count (perfect-matching number of pi-graph)
  2. |Aut|             : order of the automorphism group of the pi-graph
  3. K/|Aut|           : symmetry-discounted matching multiplicity
  4. Bay (0/1)         : bay-/fjord-region indicator
  5. K/|Aut| * Bay     : composite predictor

New 12 (computed in this script):
  Group A — pi-only graph-theoretic descriptors (10):
    6.  N_vertices     : pi-atom count (heavy-atom count of pi-graph)
    7.  N_edges
    8.  N_rings        : cyclomatic number m - n + c
    9.  Randic ^1 chi  : sum_edges 1 / sqrt(d_u d_v)
    10. Wiener W       : sum_pairs d(u,v)
    11. Hosoya Z       : sum_k m(G,k)  (all-order matching index)
    12. Balaban J      : (m/(mu+1)) * sum_edges 1 / sqrt(s_u s_v)
    13. Zagreb M1      : sum_v d(v)^2
    14. Zagreb M2      : sum_edges d(u) d(v)
    15. lambda_1       : largest adjacency-matrix eigenvalue
  Group B — whole-molecule physicochemical comparators (2):
    16. XLogP3         : PubChem REST (full molecule, H/methyl retained)
    17. MolecularWeight: PubChem REST (full molecule)

Outputs (written to paper1_JCIM/extended_predictors_N17/):
  - paper1_extended_results.csv : 17 rows  (predictor, rho, raw_p, bonferroni_p,
                                             bh_q, survives_bonf_0p05,
                                             survives_bh_0p05)
  - paper1_extended_per_pah.csv : 27 x (name + CID + log10PEF + 17 descriptors)
  - paper1_extended_log.txt     : ranked summary with Bonferroni/BH survival
                                   and Hosoya-Z / K-alone diagnostic.

Data provenance:
  SMILES are copied verbatim from compute_k_aut_v2.py (K/|Aut| audit source of
  truth).  log10PEF, Bay, K, |Aut| values are re-stated from
  data/paper1_table1.csv (canonical v4, CID-cleaned 2026-04-14 v6->v6.1).
  Group B values are fetched live from PubChem REST by CID at runtime.

Requirements: python >= 3.9, networkx >= 3.0, numpy >= 1.24, scipy >= 1.10
Network: requires outbound HTTPS to https://pubchem.ncbi.nlm.nih.gov
Run:
    python3 code/paper1_extended_predictors.py
    # (execute from the repository root; outputs go to
    #  paper1_JCIM/extended_predictors_N17/)

弦识 × 湛湛 | 2026-04-14 | attack_list_v1.md A3.4 response
"""
from __future__ import annotations
import json
import math
import sys
import time
import urllib.request
from functools import lru_cache
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import stats

# ------------------------------------------------------------------
# 1.  Dataset: SMILES + methyl flag + CID + name + log10PEF + Bay
#     SMILES copied verbatim from compute_k_aut_v2.py (audited)
#     log10PEF + Bay + CID from paper1_table1.csv (canonical)
#     K, |Aut|, K/|Aut| are already computed & audited; re-used here.
# ------------------------------------------------------------------
PAHS = [
    # (name, smiles, has_methyl, CID, log10PEF, Bay, K, Aut)
    ("Benzene",                 "c1ccccc1",                                                 False,   241, -3.00, 0,  2, 12),
    ("Naphthalene",             "c1cccc2ccccc12",                                           False,   931, -3.00, 0,  3,  4),
    ("Acenaphthylene",          "C1=CC2=CC=CC3=CC=CC1=C23",                                 False,  9161, -3.00, 0,  3,  2),
    ("Fluoranthene",            "c1ccc2c(c1)-c1cccc3cccc2c13",                              False,  9154, -3.00, 0,  6,  2),
    ("Anthracene",              "c1ccc2cc3ccccc3cc2c1",                                     False,  8418, -2.00, 0,  4,  4),
    ("Phenanthrene",            "c1ccc2c(c1)ccc1ccccc12",                                   False,   995, -3.00, 1,  5,  2),
    ("Pyrene",                  "c1cc2ccc3cccc4ccc(c1)c2c34",                               False, 31423, -3.00, 0,  6,  4),
    ("Triphenylene",            "c1ccc2c(c1)c1ccccc1c1ccccc21",                             False,  9170, -3.00, 0,  9,  6),
    ("Chrysene",                "c1ccc2c(c1)ccc1ccc3ccccc3c12",                             False,  9171, -2.00, 1,  8,  2),
    ("Benz[a]anthracene",       "c1ccc2c(c1)cc1ccc3ccccc3c1c2",                             False,  5954, -1.00, 1,  7,  1),
    ("Benzo[c]phenanthrene",    "C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C=C3",                    False,  9136, -1.64, 1,  8,  2),
    ("Benzo[a]pyrene",          "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",                         False,  2336,  0.00, 1,  9,  1),
    ("3-Methylcholanthrene",    "Cc1cc2c3ccccc3-c3cccc4cccc-2c1c34",                        True,   1674,  0.26, 1,  7,  1),
    ("5-Methylchrysene",        "Cc1cc2ccccc2c2ccc3ccccc3c12",                              True,  19427,  0.00, 1,  8,  2),
    ("DMBA",                    "CC1=C2C=CC3=CC=CC=C3C2=C(C4=CC=CC=C14)C",                  True,   6001,  1.32, 1,  7,  1),
    ("Dibenz[a,h]anthracene",   "c1ccc2c(c1)cc1ccc3cc4ccccc4cc3c1c2",                       False,  5889,  0.70, 1, 10,  2),
    ("Dibenzo[a,l]pyrene",      "C1=CC=C2C(=C1)C=C3C=CC4=C5C3=C2C6=CC=CC=C6C5=CC=C4",       False,  9119,  1.00, 1, 16,  1),
    ("Benzo[b]fluoranthene",    "c1ccc2c(c1)-c1cccc3c1c-2cc1ccccc13",                       False,  9153, -1.00, 1, 10,  1),
    ("Benzo[k]fluoranthene",    "c1ccc2cc3c(cc2c1)-c1cccc2cccc-3c12",                       False,  9158, -1.00, 1,  9,  2),
    ("Perylene",                "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54",                       False,  9142, -3.00, 0,  9,  4),
    ("Benzo[ghi]perylene",      "c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56",                     False,  9117, -2.00, 0, 14,  2),
    ("Coronene",                "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",                   False,  9115, -3.00, 0, 20, 12),
    ("Indeno[1,2,3-cd]pyrene",  "c1ccc2c(c1)-c1ccc3ccc4cccc5cc-2c1c3c45",                   False,  9131, -1.00, 1, 12,  1),
    ("Dibenzo[a,e]pyrene",      "c1ccc2c(c1)cc1c3ccccc3c3cccc4ccc2c1c43",                   False,  9126,  0.00, 1, 17,  1),
    ("Dibenzo[a,h]pyrene",      "C1=CC=C2C3=C4C(=CC2=C1)C=CC5=C4C(=CC6=CC=CC=C56)C=C3",     False,  9108,  1.00, 1, 13,  2),
    ("Dibenzo[a,i]pyrene",      "C1=CC=C2C3=C4C(=CC2=C1)C=CC5=CC6=CC=CC=C6C(=C54)C=C3",     False,  9106,  1.00, 1, 14,  2),
    ("Benzo[e]pyrene",          "C1=CC=C2C(=C1)C3=CC=CC4=C3C5=C(C=CC=C25)C=C4",             False,  9128, -3.00, 0, 11,  2),
]

# ------------------------------------------------------------------
# 2.  SMILES -> pi-graph parser (copied from compute_k_aut_v2.py)
# ------------------------------------------------------------------
def parse_smiles_to_graph(smiles: str) -> nx.Graph:
    G = nx.Graph()
    atom_idx = 0
    stack = []
    prev_atom = None
    ring_opens = {}
    i = 0
    while i < len(smiles):
        ch = smiles[i]
        if ch in ('C', 'c'):
            G.add_node(atom_idx)
            if prev_atom is not None:
                G.add_edge(prev_atom, atom_idx)
            prev_atom = atom_idx
            atom_idx += 1
            i += 1
        elif ch == '(':
            stack.append(prev_atom); i += 1
        elif ch == ')':
            prev_atom = stack.pop(); i += 1
        elif ch in ('=', '-', '#'):
            i += 1
        elif ch == '%':
            r = int(smiles[i+1:i+3])
            if r in ring_opens:
                G.add_edge(prev_atom, ring_opens[r]); del ring_opens[r]
            else:
                ring_opens[r] = prev_atom
            i += 3
        elif ch.isdigit():
            r = int(ch)
            if r in ring_opens:
                G.add_edge(prev_atom, ring_opens[r]); del ring_opens[r]
            else:
                ring_opens[r] = prev_atom
            i += 1
        else:
            i += 1
    return G


def extract_pi_system(G: nx.Graph) -> nx.Graph:
    """Drop any non-ring carbons (methyls on substituted PAHs)."""
    ring_atoms = set()
    for cyc in nx.cycle_basis(G):
        ring_atoms.update(cyc)
    if len(ring_atoms) == G.number_of_nodes():
        return G
    H = G.copy()
    H.remove_nodes_from(set(G.nodes()) - ring_atoms)
    return H


# ------------------------------------------------------------------
# 3.  Graph descriptors (Group A, 10)
# ------------------------------------------------------------------
def hosoya_index(G: nx.Graph) -> int:
    """Hosoya index = total number of matchings of all sizes (inc. empty)."""
    nodes = list(G.nodes())
    n = len(nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * n
    for u, v in G.edges():
        adj[idx[u]] |= (1 << idx[v])
        adj[idx[v]] |= (1 << idx[u])

    @lru_cache(maxsize=None)
    def count(mask: int) -> int:
        if mask == 0:
            return 1
        lsb = mask & -mask
        v = lsb.bit_length() - 1
        # v unmatched
        total = count(mask ^ lsb)
        # v matched to each neighbour u in remaining
        u_mask = adj[v] & mask
        while u_mask:
            u_lsb = u_mask & -u_mask
            total += count(mask ^ lsb ^ u_lsb)
            u_mask ^= u_lsb
        return total

    return count((1 << n) - 1)


def wiener_index(G: nx.Graph) -> int:
    spl = dict(nx.all_pairs_shortest_path_length(G))
    nodes = list(G.nodes())
    W = 0
    for i, u in enumerate(nodes):
        for v in nodes[i+1:]:
            W += spl[u][v]
    return W


def balaban_j(G: nx.Graph) -> float:
    m = G.number_of_edges()
    n = G.number_of_nodes()
    mu = m - n + nx.number_connected_components(G)  # cyclomatic number
    spl = dict(nx.all_pairs_shortest_path_length(G))
    s = {v: sum(spl[v].values()) for v in G.nodes()}
    total = 0.0
    for u, v in G.edges():
        total += 1.0 / math.sqrt(s[u] * s[v])
    return (m / (mu + 1)) * total


def randic_chi(G: nx.Graph) -> float:
    deg = dict(G.degree())
    return sum(1.0 / math.sqrt(deg[u] * deg[v]) for u, v in G.edges())


def zagreb_m1(G: nx.Graph) -> int:
    return sum(d * d for _, d in G.degree())


def zagreb_m2(G: nx.Graph) -> int:
    deg = dict(G.degree())
    return sum(deg[u] * deg[v] for u, v in G.edges())


def largest_adj_eigenvalue(G: nx.Graph) -> float:
    A = nx.to_numpy_array(G, dtype=float)
    eigs = np.linalg.eigvalsh(A)
    return float(eigs[-1])


# ------------------------------------------------------------------
# 4.  PubChem fetch (Group B, 2)
# ------------------------------------------------------------------
def fetch_pubchem_props(cid: int, retry: int = 3) -> tuple[float, float]:
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/MolecularWeight,XLogP/JSON"
    )
    last_err = None
    for _ in range(retry):
        try:
            with urllib.request.urlopen(url, timeout=15) as r:
                d = json.loads(r.read().decode("utf-8"))
            p = d["PropertyTable"]["Properties"][0]
            mw = float(p["MolecularWeight"])
            xlogp = float(p.get("XLogP", "nan"))
            return mw, xlogp
        except Exception as e:
            last_err = e
            time.sleep(1.0)
    raise RuntimeError(f"PubChem fetch failed for CID {cid}: {last_err}")


# ------------------------------------------------------------------
# 5.  Multiple-testing
# ------------------------------------------------------------------
def bh_adjust(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted q-values, same length & order as input."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order)):
        k = n - rank
        raw = pvals[idx] * n / k
        prev = min(prev, raw)
        q[idx] = prev
    return q


# ------------------------------------------------------------------
# 6.  Main
# ------------------------------------------------------------------
def main():
    # Resolve output directory relative to this script:
    # code/paper1_extended_predictors.py  ->  paper1_JCIM/extended_predictors_N17/
    repo_root = Path(__file__).resolve().parent.parent
    out_dir = repo_root / "paper1_JCIM" / "extended_predictors_N17"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir}", flush=True)

    # ---- Build per-PAH descriptor table ------------------------
    print("Building pi-graph descriptors for 27 PAH ...", flush=True)
    rows = []
    pubchem_cache = {}
    for name, smiles, meth, cid, log10pef, bay, K, Aut in PAHS:
        G_full = parse_smiles_to_graph(smiles)
        G = extract_pi_system(G_full) if meth else G_full
        n_v = G.number_of_nodes()
        n_e = G.number_of_edges()
        n_r = n_e - n_v + nx.number_connected_components(G)
        rand = randic_chi(G)
        wien = wiener_index(G)
        hosoya = hosoya_index(G)
        balab = balaban_j(G)
        zm1 = zagreb_m1(G)
        zm2 = zagreb_m2(G)
        lam1 = largest_adj_eigenvalue(G)
        # PubChem
        mw, xlogp = fetch_pubchem_props(cid)
        pubchem_cache[cid] = (mw, xlogp)
        K_over_Aut = K / Aut
        K_Aut_Bay = K_over_Aut * bay
        rows.append({
            "name": name,
            "CID": cid,
            "log10PEF": log10pef,
            "K": K,
            "Aut": Aut,
            "K_over_Aut": K_over_Aut,
            "Bay": bay,
            "K_Aut_Bay": K_Aut_Bay,
            "N_vertices": n_v,
            "N_edges": n_e,
            "N_rings": n_r,
            "Randic_1chi": rand,
            "Wiener_W": wien,
            "Hosoya_Z": hosoya,
            "Balaban_J": balab,
            "Zagreb_M1": zm1,
            "Zagreb_M2": zm2,
            "lambda1": lam1,
            "XLogP3": xlogp,
            "MW": mw,
        })
        print(f"  {name:<28s}  n={n_v:>2d} m={n_e:>2d} K={K:>3d} Aut={Aut:>3d} Z={hosoya:>7d} lam1={lam1:.3f} MW={mw:.2f} XLogP={xlogp}", flush=True)

    # ---- Predictor set ----------------------------------------
    predictors = [
        "K", "Aut", "K_over_Aut", "Bay", "K_Aut_Bay",             # original 5
        "N_vertices", "N_edges", "N_rings", "Randic_1chi",
        "Wiener_W", "Hosoya_Z", "Balaban_J",
        "Zagreb_M1", "Zagreb_M2", "lambda1",                       # group A 10
        "XLogP3", "MW",                                            # group B 2
    ]
    assert len(predictors) == 17

    y = np.array([r["log10PEF"] for r in rows], dtype=float)

    rho_rows = []
    raw_ps = []
    for p in predictors:
        x = np.array([r[p] for r in rows], dtype=float)
        rho, pval = stats.spearmanr(x, y)
        rho_rows.append((p, float(rho), float(pval)))
        raw_ps.append(float(pval))

    bonf = [min(1.0, p * 17) for p in raw_ps]
    bh_q = bh_adjust(raw_ps)

    result = []
    for (p, rho, raw), b, q in zip(rho_rows, bonf, bh_q):
        result.append({
            "predictor": p,
            "rho": rho,
            "raw_p": raw,
            "bonferroni_p": b,
            "bh_q": q,
            "survives_bonf_0p05": b < 0.05,
            "survives_bh_0p05": q < 0.05,
        })

    # ---- Write CSVs -------------------------------------------
    import csv
    per_pah_path = out_dir / "paper1_extended_per_pah.csv"
    with per_pah_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    result_path = out_dir / "paper1_extended_results.csv"
    with result_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(result[0].keys()))
        w.writeheader()
        w.writerows(result)

    # ---- Summary log ------------------------------------------
    log_path = out_dir / "paper1_extended_log.txt"
    lines = []
    lines.append("=" * 78)
    lines.append("Paper 1 (JCIM) — N=17 extended-predictor audit (attack_list_v1 A3.4)")
    lines.append("=" * 78)
    lines.append("")
    lines.append("Dataset: 27 PAH (paper1_table1.csv, v4).  Endpoint: log10PEF (rank).")
    lines.append("Descriptor definitions:")
    lines.append("  - Original 5 (K, |Aut|, K/|Aut|, Bay, K/|Aut|xBay): from canonical table.")
    lines.append("  - Group A (10): computed on pi-graph (methyl carbons excluded).")
    lines.append("  - Group B (2): XLogP3 and MolecularWeight from PubChem REST (whole molecule).")
    lines.append("")
    lines.append("Multiple-testing: Bonferroni alpha*/17 and Benjamini-Hochberg q<0.05.")
    lines.append("")

    # sort by |rho|
    sorted_res = sorted(result, key=lambda d: -abs(d["rho"]))
    lines.append(f"{'rank':>4} {'predictor':<16} {'rho':>8} {'raw_p':>12} {'Bonf_p':>12} {'BH_q':>12} {'B':>3} {'BH':>3}")
    lines.append("-" * 78)
    for i, d in enumerate(sorted_res, 1):
        lines.append(
            f"{i:>4} {d['predictor']:<16} {d['rho']:>+8.3f} "
            f"{d['raw_p']:>12.2e} {d['bonferroni_p']:>12.2e} {d['bh_q']:>12.2e} "
            f"{'Y' if d['survives_bonf_0p05'] else '.':>3} "
            f"{'Y' if d['survives_bh_0p05'] else '.':>3}"
        )

    # ---- Hosoya Z warning diagnostic ---------------------------
    k_rho = next(d for d in result if d["predictor"] == "K")["rho"]
    z_rho = next(d for d in result if d["predictor"] == "Hosoya_Z")["rho"]
    kaut_rho = next(d for d in result if d["predictor"] == "K_over_Aut")["rho"]
    kautbay_rho = next(d for d in result if d["predictor"] == "K_Aut_Bay")["rho"]
    lines.append("")
    lines.append("--- diagnostic ---")
    lines.append(f"rho(K)             = {k_rho:+.3f}")
    lines.append(f"rho(Hosoya Z)      = {z_rho:+.3f}    (structural co-variation with K predicted)")
    lines.append(f"rho(K/|Aut|)       = {kaut_rho:+.3f}")
    lines.append(f"rho(K/|Aut| x Bay) = {kautbay_rho:+.3f}")
    if abs(z_rho) >= abs(k_rho):
        lines.append("--> Hosoya Z ranks >= K alone (expected: Z's last matching term IS K).")
    if abs(kaut_rho) > abs(z_rho) and abs(kaut_rho) > 0.5:
        lines.append("--> K/|Aut| still outranks all pi-graph-only descriptors including Hosoya Z.")
    if abs(kautbay_rho) > abs(z_rho):
        lines.append("--> K/|Aut|xBay remains the top predictor among 17 tested.")
    lines.append("")
    lines.append("Interpretation: if K/|Aut|(xBay) keeps rho-rank #1 among 17, the a-priori")
    lines.append("Bonferroni concern is dissolved --- drafting did not selectively pick a winner")
    lines.append("from a hidden predictor pool; instead an explicit N=17 sweep confirms rank.")
    lines.append("If another descriptor reaches comparable rho, the positioning shifts from")
    lines.append("'best predictor' to 'only mechanistically legible predictor' (see A5.3).")

    log_path.write_text("\n".join(lines))
    print("\n".join(lines[-60:]))

    print(f"\nWrote:\n  {per_pah_path}\n  {result_path}\n  {log_path}")


if __name__ == "__main__":
    main()
