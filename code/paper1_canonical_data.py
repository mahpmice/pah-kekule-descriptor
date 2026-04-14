"""
paper1_canonical_data.py — Paper 1 (PAH Kekulé descriptor) 的 **唯一机器数据源**

规则（读这行之前不要动代码）：
-----------------------------------------------------------------------
1. 本模块的 `DATA` tuple 是 27 条 PAH 条目的权威机器可读形态：
   (name, K, Aut, PEF, bay, meth, PEF_source)
   — bay / meth 是 Python bool，不是字符串
   — K / Aut 是 int，K/|Aut| 在下游用精确除法重建（不读截断小数）
   — PEF 是 float（linear），log10PEF 在下游取

2. 下游 **所有** Paper 1 审校 / 敏感性 / 扩展脚本必须：
       from paper1_canonical_data import (
           DATA, names, K_vals, Aut_vals, PEF_vals,
           bay, meth, K_Aut, logPEF, n,
       )
   禁止再从 paper1_table1.csv 解析数据——CSV 是人读展示版本（含 "Y¹" 之类
   脚注上标），机器读取会在 bay 字段解析上分叉。已在 2026-04-14 被切过一次。

3. paper1_table1.csv 保留为 Supporting Information / 排版展示，不作为
   数据源。若 canonical DATA 与 CSV 数值分歧，canonical 胜出，CSV 需回修。

4. 新增 / 修改 PAH 条目只在本文件改。改动后下游重跑确认数字不漂移。

Why this file exists (2026-04-14 refactor): a previous methyl sensitivity script
从 CSV 解析 bay 字段，对 Benzo[c]phenanthrene 的 "Y¹"（上标脚注=Durant 1996
MEF proxy 标记）做了 `== "Y"` 比较，得到 bay=0，导致 K/|Aut|×Bay 相关系数
从 +0.834 掉到 +0.812。Debug 用了 8 个诊断脚本才定位到字符串解析。根本
症状是数据源分叉（audit.py 硬编码 vs CSV 解析），不是 "Y¹" parser 不够鲁棒。
唯一修复方式是单点数据源。
"""

from __future__ import annotations

import numpy as np

# ============================================================
# CANONICAL DATA — 27 PAHs
# ============================================================
# Column order: (name, K, Aut, PEF, bay, meth, PEF_source)
#   K, Aut: integers, verified by compute_k_aut_v2.py (2026-04-05)
#   PEF:    linear (not log); PEF=0.001 = Tier 4 assignment for
#           non-classified / IARC Group 3 compounds
#   bay:    True if the molecule contains a bay or fjord region
#   meth:   True if the molecule contains a methyl substituent (excluded
#           from K/|Aut| calculation — only the parent skeleton counts)
#   PEF_source: provenance string (audit trail, not used in analysis)
# ============================================================

DATA = [
    # name,                    K,  Aut,  PEF,      bay,   meth,  PEF_source
    ("Benzene",                2,  12,   0.001,    False, False, "Tier4: no carcin. evidence; assigned 0.001"),
    ("Naphthalene",            3,  4,    0.001,    False, False, "N&L 1992 Table 1"),
    ("Acenaphthylene",         3,  2,    0.001,    False, False, "N&L 1992 Table 1"),
    ("Fluoranthene",           6,  2,    0.001,    False, False, "N&L 1992 Table 1"),
    ("Anthracene",             4,  4,    0.01,     False, False, "N&L 1992 Table 1"),
    ("Phenanthrene",           5,  2,    0.001,    True,  False, "N&L 1992 Table 1"),
    ("Pyrene",                 6,  4,    0.001,    False, False, "N&L 1992 Table 1"),
    ("Triphenylene",           9,  6,    0.001,    False, False, "Tier4: IARC Group 3; assigned 0.001"),
    ("Chrysene",               8,  2,    0.01,     True,  False, "N&L 1992 Table 1"),
    ("Benz[a]anthracene",      7,  1,    0.1,      True,  False, "N&L 1992 Table 1"),
    # BcP: PubChem CID 9136 (NOT 9101 which is bicyclo[2.1.0]pentane).
    # K=8, |Aut|=2 verified by compute_k_aut_v2 (2026-04-05). Graph-isomorphic to chrysene.
    # K/|Aut|=4.00 < 5.0 threshold → correctly classified as TN (non-carcinogen).
    # No official cancer PEF. MEF=0.023 from Durant 1996 used as proxy.
    # CSV table1 marks this row's Bay cell as "Y¹" (¹ = PEF_source footnote),
    # which is a DISPLAY artifact — the bay-region geometry is True.
    ("Benzo[c]phenanthrene",   8,  2,    0.023,    True,  False, "Durant 1996 MEF (proxy; no IARC PEF)"),
    ("Benzo[a]pyrene",         9,  1,    1.0,      True,  False, "N&L 1992 reference compound (PEF≡1)"),
    # 3-MC: OEHHA CSF = 21.96 (mg/kg/day)⁻¹; PEF = 21.96/12 = 1.83
    ("3-Methylcholanthrene",   7,  1,    1.83,     True,  True,  "OEHHA 2011 CSF/BaP_CSF (21.96/12)"),
    # 5-MC: OEHHA CSF = 12 (mg/kg/day)⁻¹; PEF = 12/12 = 1.0
    ("5-Methylchrysene",       8,  2,    1.0,      True,  True,  "OEHHA 2011 CSF/BaP_CSF (12/12)"),
    # DMBA: OEHHA CSF = 249.96 (mg/kg/day)⁻¹; PEF = 249.96/12 = 20.83
    ("DMBA",                   7,  1,    20.83,    True,  True,  "OEHHA 2011 CSF/BaP_CSF (249.96/12)"),
    ("Dibenz[a,h]anthracene",  10, 2,    5.0,      True,  False, "N&L 1992 Table 1"),
    ("Dibenzo[a,l]pyrene",     16, 1,    10.0,     True,  False, "Collins et al. 1998 Table 2"),
    ("Benzo[b]fluoranthene",   10, 1,    0.1,      True,  False, "N&L 1992 Table 1"),
    ("Benzo[k]fluoranthene",   9,  2,    0.1,      True,  False, "N&L 1992 Table 1"),
    ("Perylene",               9,  4,    0.001,    False, False, "Tier4: IARC Group 3; assigned 0.001"),
    ("Benzo[ghi]perylene",     14, 2,    0.01,     False, False, "N&L 1992 Table 1"),
    ("Coronene",               20, 12,   0.001,    False, False, "Tier4: not classifiable; assigned 0.001"),
    ("Indeno[1,2,3-cd]pyrene", 12, 1,    0.1,      True,  False, "N&L 1992 Table 1"),
    ("Dibenzo[a,e]pyrene",     17, 1,    1.0,      True,  False, "Collins et al. 1998 Table 2"),
    ("Dibenzo[a,h]pyrene",     13, 2,    10.0,     True,  False, "Collins et al. 1998 Table 2"),
    ("Dibenzo[a,i]pyrene",     14, 2,    10.0,     True,  False, "Collins 1998"),
    # BeP: PubChem CID 9128; K=11, |Aut|=2 verified by compute_k_aut_v2.py (2026-04-05)
    # IARC Group 3; no bay region; false positive at K/|Aut|>=5.0 threshold
    ("Benzo[e]pyrene",         11, 2,    0.001,    False, False, "Tier4: IARC Group 3; assigned 0.001"),
]

# ============================================================
# EXTRACTED ARRAYS — downstream import these, not DATA tuples
# ============================================================
names    = [d[0] for d in DATA]
K_vals   = np.array([d[1] for d in DATA], dtype=float)
Aut_vals = np.array([d[2] for d in DATA], dtype=float)
PEF_vals = np.array([d[3] for d in DATA], dtype=float)
bay      = np.array([d[4] for d in DATA])            # dtype=bool
meth     = np.array([d[5] for d in DATA])            # dtype=bool

K_Aut    = K_vals / Aut_vals                         # precise float division
logPEF   = np.log10(PEF_vals)

n = len(DATA)

# ============================================================
# SELF-CHECK — run `python paper1_canonical_data.py` to verify integrity
# ============================================================
if __name__ == "__main__":
    assert n == 27, f"Expected 27 PAHs, got {n}"
    assert bay.dtype == bool, "bay must be a bool array"
    assert meth.dtype == bool, "meth must be a bool array"
    assert np.all(K_vals > 0), "All K values must be positive"
    assert np.all(Aut_vals > 0), "All |Aut| values must be positive"
    assert np.all(PEF_vals > 0), "All PEF values must be positive (log10 requires >0)"

    # Known canonical values — if these change, something broke.
    bcp_idx = names.index("Benzo[c]phenanthrene")
    assert bay[bcp_idx] == True, "BcP bay must be True (CSV 'Y¹' is display-only)"

    benzene_idx = names.index("Benzene")
    assert abs(K_Aut[benzene_idx] - 2/12) < 1e-10, "Benzene K/|Aut| precision check"

    n_bay = int(bay.sum())
    n_meth = int(meth.sum())
    print(f"paper1_canonical_data.py — integrity OK")
    print(f"  n = {n} PAHs")
    print(f"  bay-region PAHs: {n_bay} / {n}")
    print(f"  methylated PAHs: {n_meth} / {n}")
    print(f"  K/|Aut| range:   [{K_Aut.min():.3f}, {K_Aut.max():.3f}]")
    print(f"  log10(PEF) range: [{logPEF.min():+.2f}, {logPEF.max():+.2f}]")
