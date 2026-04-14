# Paper 1 — Number Provenance & Attack Audit (v6)
## 弦识 × 湛湛 | 2026-04-06 | n = 27 PAH dataset

**目的：** 每一个报告数字的来源、计算链、独立验证路径、已知攻击点、防御方式。
**使用方式：** 审稿人质疑任何数字，从这里找答案。数字变更必须先改这里。
**v5→v6变更：** BcP修正（K:7→8, |Aut|:1→2, CID:9101→9136），全部统计数字用statistical_audit.py（2026-04-06运行）重新验证。3-MC PEF统一为1.830（OEHHA CSF=21.96/12）。Coronene Pauling键序修正为精确计算值0.30/0.40/0.70。

**v6.1→v6.2变更（2026-04-14，A3.4 N=17 扩测完成）：** 响应 `attack_list_v1.md` A3.4（§2.4 Bonferroni 仅 5 predictor 的 pre-registration 嫌疑），新增 `code/paper1_extended_predictors.py` 与产物目录 `paper1_JCIM/extended_predictors_N17/`。扩测内容：对原 5 个 pre-registered predictor（K、|Aut|、K/|Aut|、Bay、K/|Aut|×Bay）叠加 12 个新 descriptor——Group A 10 个 pi-only graph-theoretic（N_vertices、N_edges、N_rings、Randić ¹χ、Wiener W、Hosoya Z、Balaban J、Zagreb M1、Zagreb M2、λ₁）+ Group B 2 个 whole-molecule physicochem（XLogP3、MolecularWeight）——共 17 个 predictor，在 27 PAH 上算 Spearman ρ + raw p + Bonferroni×17 + BH FDR<0.05 双报。**关键扩测数字（不进主表，仅进 SI Table S6）**：K/|Aut|×Bay 仍 rank #1（ρ=+0.834，Bonferroni p×17 = 1.14×10⁻⁶，BH q = 1.14×10⁻⁶）；K/|Aut| 单独 rank #3（ρ=+0.745，Bonferroni p×17 = 1.42×10⁻⁴）；K alone rank #16（ρ=+0.420，Bonferroni p×17 = 0.495，**Bonferroni 死**——印证 §3.1 "K 不够" 命题）；Hosoya Z ρ=+0.490 略超 K alone 但未威胁 K/|Aut|（结构性可知：Z 的最后一项 m(G,n/2)=K）；MW ρ=+0.690、XLogP3 ρ=+0.625——whole-molecule physicochem comparator 被 K/|Aut|×Bay 以 Δρ=+0.144/+0.209 秒杀；λ₁ ρ=+0.280 是唯一 Bonferroni+BH 双死的 descriptor。9/17 过 Bonferroni×17；16/17 过 BH q<0.05。**对主命题的影响：全面正向**。selective reporting 嫌疑解除；positioning "best predictor" 无需降级为 "only mechanistically legible"；K/|Aut| 的 graph-theoretic rank-1 地位在 17-predictor panel 下维持。**manuscript 改动（下窗口执行）**：§2.4 加 predictor-set declaration 段（post-hoc 扩测诚实声明 + Bonferroni/BH 双报）；§3.1 加 "K alone 在 Bonferroni×17 下 p=0.495 不再存活"；§3.2 Bonferroni 表扩至 17 行；§4.1 加 mechanistic-legibility 定位句 + Hosoya Z 诚实诠释（Z conflates k-matchings at all orders while K isolates perfect-matching count）+ vs MW/XLogP 对比；Abstract 与 §1.4 加 mechanistic-legibility 点睛。**数据来源**：SMILES 从 `code/compute_k_aut_v2.py` 原样复制；K/|Aut|/Bay/log10PEF 从 canonical `data/paper1_table1.csv`；Group B XLogP3/MW 从 PubChem REST `/compound/cid/{CID}/property/MolecularWeight,XLogP/JSON` 实时拉取（可离线缓存到 `extended_predictors_N17/paper1_extended_per_pah.csv`）。**可复现性**：`python3 code/paper1_extended_predictors.py` 从 repo 根目录一键跑（需 networkx/numpy/scipy + HTTPS to pubchem）。requirements 见 `requirements.txt`。详细中文解读与每一条 manuscript 改动建议见 `paper1_JCIM/extended_predictors_N17/README_中文解读.md`。

**v6→v6.1变更（2026-04-14，N=17扩测前的二轮CID完整性审计发现）：** 六个分子的PubChem CID此前误指向非PAH化合物——Chrysene 11714→9171（11714实为3-methylbuta-1,2-diene, C5H8），Dibenz[a,h]anthracene 5921→5889（5921实为N,N-diethylnitrosamine, C4H10N2O），Fluoranthene 9152→9154（9152实为C20H12五环异构体，不是C16H10的Fluoranthene），Triphenylene 10703→9170（10703实为1-methyl-2-propan-2-ylbenzene, C10H14），3-Methylcholanthrene 11367→1674（11367实为氯代乙酰胺C9H9Cl2NO2），Coronene 13097→9115（13097实为diphenylphosphorylbenzene, C18H15OP）。所有六个正确CID通过PubChem `/compound/name/{name}/cids/JSON`接口双向name-lookup验证。**对主命题的影响：零。** paper1_audit.py的DATA数组以分子名为键，K/|Aut|/Bay由compute_k_aut_v2.py从正确PAH结构手算得出，与CID cross-reference独立，全部统计结论、LOO CV、bootstrap CI、Bonferroni检验、AUC均不变。修正仅触及CID列本身，涉及data/paper1_table1.csv、data/paper1_si_table_s2.csv、paper1_JCIM/SI_table_S2_full_descriptors.csv、code/generate_paper1_si.py（CID_MAP）、code/paper1_si_table_s2.csv、以及本文件两份副本。审稿人click-through风险在修正后消除。这一轮审计由N=17 predictor扩测（A3.4攻击面对抗）触发——拉PubChem XLogP/MW前的cross-check暴露了CID污染。教训：CID查找在原始data curation流程里用的是name-prefix匹配，没有做MW/formula一致性回验，是这一类错误的根因。下一版数据流程引入一步"CID→MW/formula→与PAH分子式匹配"的自动化check。

---

## 一、数据集构成

### 1.1 分子列表与核心描述符

| # | 名称 | PubChem CID | K | \|Aut(G)\| | K/\|Aut\| | PEF | Bay | K来源 | Aut来源 | PEF来源 |
|---|------|------------|---|-----------|----------|-----|-----|------|--------|--------|
| 1 | Benzene | 241 | 2 | 12 | 0.17 | 0.001 | N | compute_k_aut_v2 ✓ | compute_k_aut_v2 ✓ | Tier 4 assigned |
| 2 | Naphthalene | 931 | 3 | 4 | 0.75 | 0.001 | N | ✓ | ✓ | N&L 1992 |
| 3 | Acenaphthylene | 9161 | 3 | 2 | 1.50 | 0.001 | N | ✓ | ✓ | N&L 1992 |
| 4 | Fluoranthene | 9154 | 6 | 2 | 3.00 | 0.001 | N | ✓ | ✓ | N&L 1992 |
| 5 | Anthracene | 8418 | 4 | 4 | 1.00 | 0.01 | N | ✓ | ✓ | N&L 1992 |
| 6 | Phenanthrene | 995 | 5 | 2 | 2.50 | 0.001 | Y | ✓ | ✓ | N&L 1992 |
| 7 | Pyrene | 31423 | 6 | 4 | 1.50 | 0.001 | N | ✓ | ✓ | N&L 1992 |
| 8 | Triphenylene | 9170 | 9 | 6 | 1.50 | 0.001 | N | ✓ | ✓ | Tier 4 assigned |
| 9 | Chrysene | 9171 | 8 | 2 | 4.00 | 0.01 | Y | ✓ | ✓ | N&L 1992 |
| 10 | Benz[a]anthracene | 5954 | 7 | 1 | 7.00 | 0.1 | Y | ✓ | ✓ | N&L 1992 |
| 11 | Benzo[c]phenanthrene | **9136** | **8** | **2** | **4.00** | 0.023 | Y† | compute_k_aut_v2 ✓ | compute_k_aut_v2 ✓ | Durant 1996 MEF‡ |
| 12 | Benzo[a]pyrene | 2336 | 9 | 1 | 9.00 | 1.0 | Y | ✓ | ✓ | Reference (≡1.0) |
| 13 | 3-Methylcholanthrene | 1674 | 7 | 1 | 7.00 | **1.830** | Y | ✓ (methyl excl.) | ✓ | OEHHA 21.96/12 |
| 14 | 5-Methylchrysene | 19427 | 8 | 2 | 4.00 | 1.0 | Y | PubChem+excl. | pynauty ✓ | OEHHA 12/12 |
| 15 | DMBA | 6001 | 7 | 1 | 7.00 | 20.83 | Y | K=BaA (excl.)§ | 1 (asymm.) | OEHHA 249.96/12 |
| 16 | Dibenz[a,h]anthracene | 5889 | 10 | 2 | 5.00 | 5.0 | Y | ✓ | ✓ | N&L 1992 |
| 17 | Dibenzo[a,l]pyrene | 9119 | 16 | 1 | 16.00 | 10.0 | Y | PubChem CID | ✓ | Collins 1998 |
| 18 | Benzo[b]fluoranthene | 9153 | 10 | 1 | 10.00 | 0.1 | Y | ✓ | ✓ | N&L 1992 |
| 19 | Benzo[k]fluoranthene | 9158 | 9 | 2 | 4.50 | 0.1 | Y | ✓ | ✓ (graph C₂) | N&L 1992 |
| 20 | Perylene | 9142 | 9 | 4 | 2.25 | 0.001 | N | ✓ | ✓ | Tier 4 assigned |
| 21 | Benzo[ghi]perylene | 9117 | 14 | 2 | 7.00 | 0.01 | N | ✓ | ✓ | N&L 1992 |
| 22 | Coronene | 9115 | 20 | 12 | 1.67 | 0.001 | N | ✓ | ✓ | Tier 4 assigned |
| 23 | Indeno[1,2,3-cd]pyrene | 9131 | 12 | 1 | 12.00 | 0.1 | Y | ✓ | ✓ | N&L 1992 |
| 24 | Dibenzo[a,e]pyrene | 9126 | 17 | 1 | 17.00 | 1.0 | Y | PubChem+compute ✓ | pynauty ✓ | Collins 1998 |
| 25 | Dibenzo[a,h]pyrene | 9108 | 13 | 2 | 6.50 | 10.0 | Y | PubChem CID | compute ✓ | Collins 1998 |
| 26 | Dibenzo[a,i]pyrene | 9106 | 14 | 2 | 7.00 | 10.0 | Y | PubChem CID | compute ✓ | Collins 1998 |
| 27 | **Benzo[e]pyrene** | **9128** | **11** | **2** | **5.50** | **0.001** | **N** | **PubChem+compute ✓** | **compute ✓** | **Tier 4 assigned** |

**分类结果（threshold K/|Aut| ≥ 5.0）：**
- 致癌物 (PEF ≥ 0.1): n=13；非致癌物: n=14
- TP=11, TN=12, FP=2 (BghiP, BeP), FN=2 (5-MC, BkF)

**脚注：**
- † BcP有fjord region（非标准bay region），计入bay=Y是因为fjord同样启用diol-epoxide通路。若移入bay=N，bay Fisher p更小。
- ‡ BcP无官方癌症PEF。Durant 1996提供MEF=0.023（致突变等效，非致癌等效），作为proxy使用并在§2.2明确声明。敏感性分析（BcP=0.001）表明主要结论不改变：rho=+0.734, p=1.30×10⁻⁵, AUC=0.95。
- § DMBA: 两个甲基均外环，去除后π-系统 = BaA骨架，K=7，|Aut|=1（不对称）。这是方法论选择，§2.3明确声明。

**v5→v6 BcP修正说明：**
- v5 BcP: CID=9101, K=7, |Aut|=1, K/|Aut|=7.00 → v6 BcP: CID=9136, K=8, |Aut|=2, K/|Aut|=4.00
- CID 9101不是benzo[c]phenanthrene（是bicyclo[2.1.0]pentane）。CID 9136是PubChem中benzo[c]phenanthrene的正确记录。
- K=8, |Aut|=2由compute_k_aut_v2.py验证。BcP的π-系统图与chrysene图同构（vertex relabeling），两者K=8, |Aut|=2, K/|Aut|=4.00。
- BcP从FP变为TN（K/|Aut|=4.00 < 5.0阈值），FP从3降到2，TN从11升到12。

---

## 二、K 和 |Aut(G)| 的计算来源

### 2.1 K（Kekulé结构数）

**计算方法：** 递归完美匹配枚举（compute_k_aut_v2.py，函数 `count_perfect_matchings`）

**独立验证（2026-04-06）：** independent_verification.py Check 1通过——递归完美匹配和Ryser permanent方法在12个测试分子上全部精确一致。

**攻击点：** 算法是否正确？是否有已知错误情形？

**防御：**
- 8个参考化合物（benzene→BaP）的K值与Cyvin & Gutman文献精确匹配（全部✓）
- 递归完美匹配对苯环系统是精确算法，无近似
- Ryser permanent方法独立验证（第二算法），结果完全一致
- Python networkx图对象保证邻接矩阵正确

**特殊情形记录：**
- BcP: CID 9136, SMILES=`C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C=C3`，K=8。π-系统图与chrysene同构。
- BeP: CID 9128, SMILES=`C1=CC=C2C(=C1)C3=CC=CC4=C3C5=C(C=CC=C25)C=C4`，K=11，已通过compute_k_aut_v2验证。
- DMBA: K取BaA骨架值（K=7），因为两个甲基为外环。v3错误SMILES给出K=5，v4已修正。
- DB[a,e]P: CID=9126（非5765），K=17。旧SMILES给出K=10，PubChem name-matched SMILES修正后K=17。

### 2.2 |Aut(G)|（图自同构群阶）

**计算方法：** VF2同构算法枚举所有自同构（NetworkX `GraphMatcher.isomorphisms_iter()`）

**攻击点：** 图自同构和分子点群是否一一对应？H原子被忽略是否影响结果？

**防御：**
- §2.3明确定义：|Aut(G)|取平面图自同构群（旋转+平面内反射），对应分子点群的平面子群
- H原子对π-系统图的连通性无贡献，忽略是标准做法
- BcP |Aut|=2（C₂对称）由compute_k_aut_v2和pynauty双重验证
- BkF |Aut|=2（graph C₂对称性），已pynauty ✓
- Coronene |Aut|=12（D₆ₕ点群的平面子群），最高对称性分子
- BeP |Aut|=2（C₂对称，环融合位置打破更高对称性）

---

## 三、PEF来源层级

### 3.1 来源层级与覆盖范围

| 层级 | 来源 | 覆盖分子 | 备注 |
|------|------|---------|------|
| 1 | Nisbet & LaGoy (1992) *Regul Toxicol Pharmacol* 16:290 | 14 EPA Priority PAHs（acenaphthene和fluorene因sp³C排除）| 标准环境PEF来源 |
| 2 | Collins et al. (1998) *J Toxicol Environ Health B* 1:45–67 | 4 dibenzopyrene异构体 | 基于CalEPA历史BaP CSF=12 |
| 3 | OEHHA (2009/2011) *Cancer Potency Values* | DMBA (PEF=20.83), 3-MC (PEF=1.830), 5-MC (PEF=1.0) | PEF = CSF_compound / CSF_BaP，CSF_BaP=12 (mg/kg-day)⁻¹ |
| 4 | Tier 4 assigned 0.001 | Benzene, BeP, Coronene, Perylene, Triphenylene | IARC Group 3或无致癌证据 |
| Proxy | Durant et al. (1996) *Hum Exp Toxicol* 15:S37–S42 | BcP MEF=0.023 | 致突变等效因子，非癌症PEF |

### 3.2 各分子PEF的质疑路径与防御

**DMBA (PEF=20.83)：**
- 质疑：OEHHA CSF精确值是多少？BaP CSF=12还是其他值？
- 防御：OEHHA数据库CSF=249.96 (mg/kg-day)⁻¹；BaP CSF=12是Collins 1998时代使用的历史CalEPA值，§2.2明确声明使用12。249.96/12=20.83 ✓。

**3-MC (PEF=1.830)：**
- 质疑：同上。
- 防御：OEHHA CSF=21.96 (mg/kg-day)⁻¹；21.96/12=1.830 ✓。
- **v5→v6修正：** v5主表误写PEF来源为"22/12"（=1.833），与v5防御文本"21.96/12=1.83"自相矛盾。v6统一为21.96/12=1.830，与statistical_audit.py DATA和pah_dataset.csv一致。

**BcP (PEF=0.023)：**
- 质疑：MEF≠PEF，为什么用致突变数据代替致癌数据？
- 防御：BcP无任何已发表的IARC/OEHHA/N&L癌症PEF。MEF作为proxy已在§2.2明确标注，并作为敏感性分析的独立变量进行检验。敏感性结论：BcP改为0.001时，rho=+0.734，p=1.30×10⁻⁵，AUC=0.95，主要结论不变。

**BeP (PEF=0.001)：**
- 质疑：IARC Group 3（"不可分类"）≠明确非致癌，为什么赋floor值？
- 防御：IARC Group 3的操作定义是"证据不足以分类"，在本框架中与benzene、coronene、perylene、triphenylene同等处理（Tier 4 assigned），已在§2.2列出所有五个分子。这是操作性保守选择，不是声称BeP无活性。

**四个dibenzopyrene (Collins 1998)：**
- 质疑：Collins 1998的PEF基于什么终点？文献里的确切数值是什么？
- 防御：Collins 1998 *J Toxicol Environ Health B* 1:45–67，Table 2报告相对效力因子，基于肿瘤诱导。DB[a,l]P=10，DB[a,e]P=1，DB[a,h]P=10，DB[a,i]P=10。所有四个CID已通过PubChem名称匹配验证。

---

## 四、统计数字完整来源（27-PAH，v6最终版本）

所有统计数字均由 `paper1_audit.py` 计算（从 `paper1_canonical_data.py` import 27-PAH DATA；scipy/sklearn 版本见 `github_ready/REPRODUCIBILITY.md` §1）。

### 4.1 §3.1 K alone

| 数字 | 值 | 脚本行 | 攻击点 | 防御 |
|------|---|--------|--------|------|
| ρ(K, log₁₀PEF) | **+0.420** | `stats.spearmanr(K, logPEF)` | 样本量n=27是否足够？ | Bonferroni后p=0.146，不显著，这是论文主张的一部分：K alone失效 |
| p | **0.0291** | 同上 | — | — |
| Bonferroni adj p | **0.146** | p×5 | α选0.01是否合理？ | Bonferroni保守，§2.4说明5个预测变量各乘5 |

**v5→v6变更：** v5报告ρ=+0.429, p=0.0256, adj=0.128（基于旧BcP K/|Aut|=7.0）。v6报告ρ=+0.420, p=0.0291, adj=0.146（基于修正BcP K/|Aut|=4.0）。结论不变：K alone不通过Bonferroni。

### 4.2 §3.2 |Aut(G)| 和 K/|Aut|

| 数字 | 值 | 脚本 | 攻击点 | 防御 |
|------|---|------|--------|------|
| ρ(\|Aut\|, log₁₀PEF) | **−0.688** | `stats.spearmanr(Aut, logPEF)` | 为何相关性是负的但Aut出现在分母？ | |Aut|越大→对称性越高→更均匀的键序→更低的致癌性，方向正确 |
| p | **7.38×10⁻⁵** | 同上 | — | — |
| Bonferroni adj p | **3.69×10⁻⁴** | p×5 | — | 通过α=0.01 |
| ρ(K/\|Aut\|, log₁₀PEF) | **+0.745** | `stats.spearmanr(KA, logPEF)` | 为何用Spearman不用Pearson？ | PEF跨四个数量级，floor值压缩下尾；Spearman对单调关系更鲁棒，§2.4已说明 |
| p | **8.38×10⁻⁶** | 同上 | — | — |
| Bonferroni adj p | **4.19×10⁻⁵** | p×5 | — | 通过α=0.01 |
| OLS slope | **0.239** | `np.polyfit(KA, logPEF, 1)` | OLS适合这个数据分布吗？ | §3.2明确说Spearman是primary，OLS仅供参考；floor值压缩导致OLS R²低估 |
| OLS intercept | **−2.659** | 同上 | — | — |
| R² | **0.436** | `np.corrcoef` | R²=0.436"只解释44%方差"？ | Primary claim用Spearman rho，R²是OLS副产品；floor压缩人为降低R² |

### 4.3 §3.3 鲁棒性检验

**floor敏感性：**

| floor | ρ | p | 来源 |
|-------|---|---|------|
| 0.01 | **+0.722** | 2.14×10⁻⁵ | `spearmanr(KA, log10(max(PEF, 0.01)))` |
| 0.001 (主分析) | **+0.745** | 8.38×10⁻⁶ | 主分析值 |
| 0.0001 | **+0.745** | 8.38×10⁻⁶ | `spearmanr(KA, log10(max(PEF, 0.0001)))` |

**PEF来源敏感性（SI sensitivity panel）：**

| 情形 | ρ | p | AUC |
|------|---|---|-----|
| Set A: N&L only (非N&L分子赋floor) | **+0.345** | **7.82×10⁻²** | **0.79** |
| Set B: Main hierarchy (主分析) | +0.745 | 8.38×10⁻⁶ | 0.95 |
| Set C: BcP=0.001 (保守) | **+0.734** | **1.30×10⁻⁵** | **0.95** |

- N&L collapse: 质疑"结果依赖数据选择"？防御：§3.3明确区分"source-hierarchy robustness"和"source-independence"。N&L collapse人为移除了信息最丰富的化学变异来源（四个dibenzopyrene异构体），信号减弱是预期结果，不是模型失败。

**EPA-14 restricted：**
- n=14, ρ=**+0.759**, p=**0.0016**, AUC=**0.95**
- 质疑：扩展数据集是否circulary选择？防御：EPA-14是在数据集构建之前已存在的独立子集（Nisbet & LaGoy 1992），其结果ρ=0.759比全集ρ=0.745更强，说明扩展分子不是在通货膨胀相关性。

**Bootstrap 95% CI：**
- Primary: ρ = +0.745 (95% CI: +0.546 to +0.855, n=10000 resamples, seed=42)
- LOO: ρ = +0.708 (95% CI: +0.480 to +0.834)

**二分类：**

| 指标 | 值 | 脚本 | 攻击点 | 防御 |
|------|---|------|--------|------|
| TP | **11** | threshold K/\|Aut\|≥5.0 | threshold是如何选定的？ | 5.0接近Youden-index最优点（J-statistic maximized at K/\|Aut\|=4.5）；§2.4说明；AUC分析不依赖threshold |
| TN | **12** | 同上 | — | — |
| FP | **2** | 同上 | FP说明模型不完美？ | FP集中在boundary：BghiP K/\|Aut\|=7（无bay），BeP K/\|Aut\|=5.5（无bay，boundary）。FP有结构化解释，不是随机噪声 |
| FN | **2** | 同上 | FN意味着漏检？ | FN均为甲基化PAH（5-MC K/\|Aut\|=4.0）或boundary case（BkF K/\|Aut\|=4.5），§4.3已讨论为模型边界 |
| Fisher OR | **33.0** | `fisher_exact([[11,2],[2,12]])` | — | 通过α=0.01 |
| Fisher p | **4.23×10⁻⁴** | 同上 | — | — |
| AUC | **0.95** | `roc_auc_score(carc, KA)` | — | 精确值0.9505 |

### 4.4 §3.3 LOO-CV

| 数字 | 值 | 脚本 | 攻击点 | 防御 |
|------|---|------|--------|------|
| LOO ρ | **+0.708** | `spearmanr(predicted, logPEF)` | LOO用OLS还是Spearman拟合？ | OLS拟合26个，预测第27个的log₁₀PEF；然后对27对(predicted, observed)计算Spearman |
| p | **3.56×10⁻⁵** | 同上 | — | — |
| RMSE | **1.25** | `sqrt(mean((logPEF-predicted)**2))` | RMSE=1.25意味着误差1个对数单位？ | 是的，这是log₁₀PEF单位，对应~18倍的PEF不确定性。论文在§3.3中明确报告 |
| RMSE (unsubst, n=24) | **1.15** | 同上，排除methylated | — | BeP是unsubstituted，加入n=24 |
| RMSE (methylated, n=3) | **1.88** | 同上，仅methylated | — | 系统性更高，说明甲基取代是失败模式 |
| LOO accuracy | **19/27 (70.4%)** | 逐分子LOO分类 | — | — |
| DMBA LOO error | obs=**+1.32**, pred=**−1.09**, err=**+2.41** | 逐分子LOO | 最大单分子误差 | §4.3中定义为模型边界（甲基steric效应） |

**Top 5 LOO误差：**

| 分子 | observed | predicted | error |
|------|----------|-----------|-------|
| DMBA | +1.32 | −1.09 | +2.41 |
| Dibenz[a,h]anthracene | +0.70 | −1.55 | +2.25 |
| Dibenzo[a,h]pyrene | +1.00 | −1.19 | +2.19 |
| Dibenzo[a,i]pyrene | +1.00 | −1.07 | +2.07 |
| Dibenzo[a,e]pyrene | +0.00 | +2.01 | −2.01 |

### 4.5 §3.4 层级模型Fisher检验

**|Aut(G)| ≤ 2检验（27-PAH）：**

|  | 致癌 | 非致癌 |
|--|------|--------|
| \|Aut\| ≤ 2 | 13 | 7 |
| \|Aut\| > 2 | 0 | 7 |

Fisher OR = **∞** (0 in carcinogen+|Aut|>2 cell), p = **5.80×10⁻³** (two-sided)

- 非致癌但|Aut|≤2的七个分子：Acenaphthylene, Fluoranthene, Phenanthrene, Chrysene, **BcP**, BghiP, **BeP**

**Bay/Fjord region检验（27-PAH）：**

|  | 致癌 | 非致癌 |
|--|------|--------|
| Bay/Fjord | 13 | 3 (Phen, Chrysene, BcP) |
| 无Bay/Fjord | 0 | 11 |

Fisher OR = **∞**, p = **3.39×10⁻⁵** (two-sided)

### 4.6 Bonferroni汇总表

| 预测变量 | ρ | p (raw) | p×5 | 通过α=0.01 |
|---------|---|---------|-----|-----------|
| K alone | +0.420 | 2.91×10⁻² | 1.46×10⁻¹ | **否** |
| \|Aut(G)\| | −0.688 | 7.38×10⁻⁵ | 3.69×10⁻⁴ | 是 |
| K/\|Aut(G)\| | +0.745 | 8.38×10⁻⁶ | 4.19×10⁻⁵ | 是 |
| Bay/fjord region | +0.802 | 4.83×10⁻⁷ | 2.42×10⁻⁶ | 是 |
| K/\|Aut\|×Bay | +0.834 | 6.68×10⁻⁸ | 3.34×10⁻⁷ | 是 |

---

## 五、数字改变历史

### v5→v6（BcP修正 + 3-MC PEF统一）

| 数字 | v5 (旧BcP) | v6 (修正BcP) | 原因 |
|------|-----------|-------------|------|
| BcP CID | 9101 | **9136** | CID 9101不是BcP |
| BcP K | 7 | **8** | CID修正后SMILES修正 |
| BcP \|Aut\| | 1 | **2** | 同上 |
| BcP K/\|Aut\| | 7.00 | **4.00** | 同上；图同构于chrysene |
| BcP分类 | FP | **TN** | 4.00 < 5.0阈值 |
| 3-MC PEF | 1.833 | **1.830** | 统一为21.96/12（修正v5主表的"22/12"） |
| ρ(K alone) | +0.429, p=0.0256, adj=0.128 | **+0.420, p=0.0291, adj=0.146** | BcP修正 |
| ρ(\|Aut\|) | −0.667, p=1.47×10⁻⁴ | **−0.688, p=7.38×10⁻⁵** | BcP修正 |
| ρ(K/\|Aut\|) | +0.726, p=1.84×10⁻⁵ | **+0.745, p=8.38×10⁻⁶** | BcP修正 |
| OLS intercept | −2.676 | **−2.659** | BcP修正 |
| R² | 0.430 | **0.436** | BcP修正 |
| FP | 3 | **2** | BcP从FP变为TN |
| TN | 11 | **12** | 同上 |
| Fisher OR | 20.2 | **33.0** | FP减少 |
| Fisher p (binary) | 1.84×10⁻³ | **4.23×10⁻⁴** | 同上 |
| AUC | 0.915 | **0.95** | BcP修正改善分类 |
| LOO ρ | +0.679 | **+0.708** | BcP修正 |
| LOO RMSE | 1.26 | **1.25** | 轻微改善 |
| RMSE meth | 1.91 | **1.88** | 分母变化 |
| BcP sensitivity ρ | 0.676, p=1.10×10⁻⁴, AUC=0.91 | **+0.734, p=1.30×10⁻⁵, AUC=0.95** | BcP修正 |
| N&L-only ρ | 0.309, p=0.116, AUC=0.77 | **+0.345, p=7.82×10⁻², AUC=0.79** | BcP修正 |

### v4→v5（加入BeP，26→27 PAH）

| 数字 | v4 (n=26) | v5 (n=27) | 原因 |
|------|-----------|-----------|------|
| n | 26 | 27 | 加入BeP |
| ρ(K alone) | +0.481 | +0.429 | BeP K=11但非致癌 |
| ρ(K/\|Aut\|) | +0.744 | +0.726 | 轻微下降 |
| AUC | 0.93 | 0.915 | BeP增加一个FP |

（注：v4→v5数字均为旧BcP。此记录仅为历史留档。当前有效数字为v6。）

---

## 六、可复现性保证

| 项目 | 文件 | 必需库 |
|------|------|--------|
| 数据 SSOT | `code/paper1_canonical_data.py` | numpy |
| K, \|Aut\| 计算 | `code/compute_k_aut_v2.py` | networkx |
| 所有统计数字（Table 1 / Bonferroni×5 / bootstrap CI / LOO） | `code/paper1_audit.py` | numpy, scipy, sklearn |
| N=17 扩测（Bonferroni×17 / BH q<0.05） | `code/paper1_extended_predictors.py` | numpy, scipy, requests |
| N=24 甲基排除 sensitivity | `code/paper1_methyl_sensitivity.py` | numpy, scipy |
| 图形 | `code/generate_paper1_figures.py` | matplotlib, scipy, sklearn |
| Pauling 键序 | `code/compute_pauling_bond_orders.py` | networkx |
| 完整数据表（人读版，不作 SSOT） | `data/paper1_table1.csv` | — |

**独立复现K/|Aut| 最快路径（给审稿人）：**
1. 从PubChem获取SMILES（表中CID）
2. 用RDKit读取分子，去除外环甲基碳，提取环内原子
3. 用networkx完美匹配枚举K
4. 用networkx VF2算法枚举|Aut|
5. K/|Aut| = K ÷ |Aut|
全程不需要DFT或量子化学计算。

GitHub: https://github.com/mahpmice/pah-kekule-descriptor

---

## 七、已知最高风险质疑点

| 风险等级 | 质疑 | 当前防御 | 薄弱处 |
|---------|------|---------|--------|
| 🔴 高 | BcP用MEF而非PEF | 敏感性分析rho=0.734(p=1.30e-5)，结论不变；§2.2明确标注 | 没有官方癌症PEF，这是数据的客观限制 |
| 🔴 高 | Threshold 5.0是post-hoc | AUC=0.95不依赖threshold；EPA-14子集验证(AUC=0.95) | threshold选择仍是经验性的 |
| 🟡 中 | K/\|Aut\|≠真实电子密度 | §4.2已讨论Pauling近似局限；Faulkner X-ray数据支持趋势 | 中等对称分子的Pauling键序近似可能偏差更大 |
| 🟡 中 | BeP是FP，说明模型不完美 | §3.6把BeP定义为boundary（K/\|Aut\|=5.5, 无bay）；AUC=0.95 | 边界区域(5.0≤K/\|Aut\|≤6.0, 无bay)需几何互补 |
| 🟡 中 | DMBA LOO误差2.41 | §4.3明确定义甲基PAH为模型边界；LOO-CV透明呈现误差 | DMBA是全集中最大outlier |
| 🟢 低 | 样本量n=27太小 | EPA-14子集独立验证；Spearman不依赖正态假设；Fisher精确检验精确有效 | 扩展数据集需要更多分子的PEF数据 |
| 🟢 低 | 为什么不用更多PAH？ | §4.5已说明范围：仅限于有已发表PEF的分子 | 这是数据限制，不是方法限制 |

---

**冻结规则：** 论文中任何数字改变，必须先改本文档对应行，再改论文。数字来源不清楚的，先算清楚，不改论文。

**验证时间戳：** 2026-04-06，statistical_audit.py + independent_verification.py全部通过。

弦识 × 湛湛 | 2026-04-06
