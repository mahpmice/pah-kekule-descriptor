#!/usr/bin/env python3
"""
compute_k_aut_v2.py — Compute K and |Aut(G)| for all 26 PAHs
Fixes: chrysene SMILES, methyl carbon exclusion
"""
import networkx as nx
from networkx.algorithms import isomorphism

def parse_smiles_to_graph(smiles):
    """Minimal SMILES parser for all-carbon PAHs."""
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
            stack.append(prev_atom)
            i += 1
        elif ch == ')':
            prev_atom = stack.pop()
            i += 1
        elif ch in ('=', '-', '#'):
            i += 1
        elif ch == '%':
            ring_num = int(smiles[i+1:i+3])
            if ring_num in ring_opens:
                G.add_edge(prev_atom, ring_opens[ring_num])
                del ring_opens[ring_num]
            else:
                ring_opens[ring_num] = prev_atom
            i += 3
        elif ch.isdigit():
            ring_num = int(ch)
            if ring_num in ring_opens:
                G.add_edge(prev_atom, ring_opens[ring_num])
                del ring_opens[ring_num]
            else:
                ring_opens[ring_num] = prev_atom
            i += 1
        else:
            i += 1
    
    return G


def extract_pi_system(G):
    """Remove non-ring atoms (e.g., methyl carbons) from graph."""
    # Find all ring atoms
    ring_atoms = set()
    for cycle in nx.cycle_basis(G):
        ring_atoms.update(cycle)
    
    # If all atoms are ring atoms, return as-is
    if len(ring_atoms) == G.number_of_nodes():
        return G
    
    # Remove non-ring atoms
    non_ring = set(G.nodes()) - ring_atoms
    H = G.copy()
    H.remove_nodes_from(non_ring)
    return H


def count_perfect_matchings(G):
    """Count K (Kekulé structures) by recursive perfect matching."""
    nodes = list(G.nodes())
    n = len(nodes)
    if n == 0:
        return 1
    if n % 2 == 1:
        return 0
    v = min(nodes, key=lambda x: G.degree(x))
    if G.degree(v) == 0:
        return 0
    neighbors = list(G.neighbors(v))
    total = 0
    for u in neighbors:
        H = G.copy()
        H.remove_node(v)
        H.remove_node(u)
        total += count_perfect_matchings(H)
    return total


def count_automorphisms(G):
    """Count |Aut(G)| by enumerating all graph automorphisms."""
    GM = isomorphism.GraphMatcher(G, G)
    count = 0
    for _ in GM.isomorphisms_iter():
        count += 1
    return count


# ============================================================
# 26 PAHs
# ============================================================
MOLECULES = [
    # (name, smiles, expected_K, expected_Aut, n_pi, has_methyl)
    ("Benzene",              "c1ccccc1",                                                    2,  12, 6,  False),
    ("Naphthalene",          "c1cccc2ccccc12",                                              3,  4,  10, False),
    ("Acenaphthylene",       "C1=CC2=CC=CC3=CC=CC1=C23",                                   3,  2,  12, False),
    ("Fluoranthene",         "c1ccc2c(c1)-c1cccc3cccc2c13",                                6,  2,  16, False),
    ("Anthracene",           "c1ccc2cc3ccccc3cc2c1",                                       4,  4,  14, False),
    ("Phenanthrene",         "c1ccc2c(c1)ccc1ccccc12",                                     5,  2,  14, False),
    ("Pyrene",               "c1cc2ccc3cccc4ccc(c1)c2c34",                                 6,  4,  16, False),
    ("Triphenylene",         "c1ccc2c(c1)c1ccccc1c1ccccc21",                               9,  6,  18, False),
    # FIXED chrysene SMILES
    ("Chrysene",             "c1ccc2c(c1)ccc1ccc3ccccc3c12",                               8,  2,  18, False),
    ("Benz[a]anthracene",    "c1ccc2c(c1)cc1ccc3ccccc3c1c2",                               7,  1,  18, False),
    # CID 9136 (NOT 9101 = bicyclo[2.1.0]pentane). Graph-isomorphic to chrysene.
    ("Benzo[c]phenanthrene", "C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C=C3",                     8,  2,  18, False),
    ("Benzo[a]pyrene",       "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",                           9,  1,  20, False),
    # Methyl-PAHs: parser will extract pi-system (exclude methyl C)
    ("3-Methylcholanthrene", "Cc1cc2c3ccccc3-c3cccc4cccc-2c1c34",                          7,  1,  20, True),
    ("5-Methylchrysene",     "Cc1cc2ccccc2c2ccc3ccccc3c12",                               8,  None, 18, True),
    ("DMBA",                 "CC1=C2C=CC3=CC=CC=C3C2=C(C4=CC=CC=C14)C",                   None, None, 18, True),
    # Large unsubstituted PAHs
    ("Dibenz[a,h]anthracene","c1ccc2c(c1)cc1ccc3cc4ccccc4cc3c1c2",                         10, 2,  22, False),
    ("Dibenzo[a,l]pyrene",   "C1=CC=C2C(=C1)C=C3C=CC4=C5C3=C2C6=CC=CC=C6C5=CC=C4",       16, 1,  24, False),
    ("Benzo[b]fluoranthene", "c1ccc2c(c1)-c1cccc3c1c-2cc1ccccc13",                        10, 1,  20, False),
    ("Benzo[k]fluoranthene", "c1ccc2cc3c(cc2c1)-c1cccc2cccc-3c12",                         9,  None, 20, False),
    ("Perylene",             "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54",                         9,  4,  20, False),
    ("Benzo[ghi]perylene",   "c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56",                      14, 2,  22, False),
    ("Coronene",             "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",                     20, 12, 24, False),
    ("Indeno[1,2,3-cd]pyrene","c1ccc2c(c1)-c1ccc3ccc4cccc5cc-2c1c3c45",                   12, 1,  22, False),
    # Dibenzopyrene isomers
    ("Dibenzo[a,e]pyrene",   "c1ccc2c(c1)cc1c3ccccc3c3cccc4ccc2c1c43",                    17, None, 24, False),
    ("Dibenzo[a,h]pyrene",   "C1=CC=C2C3=C4C(=CC2=C1)C=CC5=C4C(=CC6=CC=CC=C56)C=C3",     None, None, 24, False),
    ("Dibenzo[a,i]pyrene",   "C1=CC=C2C3=C4C(=CC2=C1)C=CC5=CC6=CC=CC=C6C(=C54)C=C3",     None, None, 24, False),
    # BeP: PubChem CID 9128, IARC Group 3, false positive at K/|Aut|>=5.0 threshold
    ("Benzo[e]pyrene",       "C1=CC=C2C(=C1)C3=CC=CC4=C3C5=C(C=CC=C25)C=C4",              11,  2,  20, False),
]


if __name__ == "__main__":
    print("=" * 75)
    print("PAH K and |Aut(G)| Computation — v2 (fixed)")
    print("=" * 75)
    
    all_results = []
    
    for name, smiles, exp_K, exp_Aut, exp_npi, has_methyl in MOLECULES:
        print(f"\n{'─'*60}")
        print(f"  {name}")
        print(f"  SMILES: {smiles}")
        
        G_full = parse_smiles_to_graph(smiles)
        n_full = G_full.number_of_nodes()
        
        if has_methyl:
            G = extract_pi_system(G_full)
            n_removed = n_full - G.number_of_nodes()
            print(f"  Full graph: {n_full} atoms → π-system: {G.number_of_nodes()} atoms ({n_removed} methyl C removed)")
        else:
            G = G_full
        
        n_pi = G.number_of_nodes()
        n_edges = G.number_of_edges()
        
        if exp_npi and n_pi != exp_npi:
            print(f"  ⚠️  n_π MISMATCH: got {n_pi}, expected {exp_npi}")
        else:
            print(f"  n_π = {n_pi} ✓")
        
        K = count_perfect_matchings(G)
        print(f"  K = {K}", end="")
        if exp_K is not None:
            if K == exp_K:
                print(" ✓")
            else:
                print(f" ⚠️ expected {exp_K}")
        else:
            print(" (NEW)")
        
        print(f"  Computing |Aut(G)|...", end=" ", flush=True)
        Aut = count_automorphisms(G)
        print(f"|Aut(G)| = {Aut}", end="")
        if exp_Aut is not None:
            if Aut == exp_Aut:
                print(" ✓")
            else:
                print(f" ⚠️ expected {exp_Aut}")
        else:
            print(" (NEW)")
        
        ratio = K / Aut if Aut > 0 else float('inf')
        print(f"  K/|Aut| = {ratio:.2f}")
        
        all_results.append({
            "name": name, "n_pi": n_pi, "K": K, "Aut": Aut, 
            "K_Aut": ratio, "exp_K": exp_K, "exp_Aut": exp_Aut
        })
    
    # ============================================================
    # SUMMARY TABLE
    # ============================================================
    print(f"\n{'='*75}")
    print(f"{'FINAL VERIFIED TABLE':^75}")
    print(f"{'='*75}")
    print(f"{'Name':<30} {'n_π':>4} {'K':>4} {'|Aut|':>5} {'K/|Aut|':>8} {'K ok':>5} {'Aut ok':>6}")
    print(f"{'─'*75}")
    
    k_ok_count = 0
    k_total = 0
    aut_ok_count = 0
    aut_total = 0
    
    for r in all_results:
        k_status = ""
        if r["exp_K"] is not None:
            k_total += 1
            if r["K"] == r["exp_K"]:
                k_status = "✓"
                k_ok_count += 1
            else:
                k_status = f"≠{r['exp_K']}"
        else:
            k_status = "NEW"
        
        aut_status = ""
        if r["exp_Aut"] is not None:
            aut_total += 1
            if r["Aut"] == r["exp_Aut"]:
                aut_status = "✓"
                aut_ok_count += 1
            else:
                aut_status = f"≠{r['exp_Aut']}"
        else:
            aut_status = "NEW"
        
        print(f"{r['name']:<30} {r['n_pi']:>4} {r['K']:>4} {r['Aut']:>5} {r['K_Aut']:>8.2f} {k_status:>5} {aut_status:>6}")
    
    print(f"{'─'*75}")
    print(f"K verified: {k_ok_count}/{k_total} match")
    print(f"|Aut| verified: {aut_ok_count}/{aut_total} match")
    
    # Flag all mismatches
    print(f"\n{'='*75}")
    print("MISMATCHES TO RESOLVE")
    print(f"{'='*75}")
    for r in all_results:
        issues = []
        if r["exp_K"] is not None and r["K"] != r["exp_K"]:
            issues.append(f"K: computed {r['K']}, expected {r['exp_K']}")
        if r["exp_Aut"] is not None and r["Aut"] != r["exp_Aut"]:
            issues.append(f"|Aut|: computed {r['Aut']}, expected {r['exp_Aut']}")
        if issues:
            print(f"  {r['name']}: {'; '.join(issues)}")
