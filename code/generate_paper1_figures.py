#!/usr/bin/env python3
"""
Generate updated manuscript figures for Paper 1.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from rdkit import Chem
from rdkit.Chem import Draw
from scipy import stats
from sklearn.metrics import roc_curve, roc_auc_score


ROOT = Path(__file__).resolve().parent


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


audit = load_module(ROOT / "paper1_audit.py", "paper1_audit_mod")

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "bold",
        "figure.dpi": 180,
        "savefig.dpi": 300,
    }
)


names = [d[0] for d in audit.DATA]
K = np.array([d[1] for d in audit.DATA], dtype=float)
Aut = np.array([d[2] for d in audit.DATA], dtype=float)
PEF = np.array([d[3] for d in audit.DATA], dtype=float)
bay = np.array([d[4] for d in audit.DATA], dtype=bool)
meth = np.array([d[5] for d in audit.DATA], dtype=bool)
KA = K / Aut
logPEF = np.log10(PEF)


SMILES = {
    "Coronene": "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
    "Benzo[a]pyrene": "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",
    "Benzo[e]pyrene": "C1=CC=C2C(=C1)C3=CC=CC4=C3C5=C(C=CC=C25)C=C4",
}


PALETTE = {
    "base": "#1f3c88",
    "methyl": "#d1495b",
    "coronene": "#2a9d8f",
    "false": "#e76f51",
    "line": "#222222",
    "bg": "#faf8f3",
}


def point_style(name: str):
    if name == "Coronene":
        return dict(color=PALETTE["coronene"], marker="D", s=70, edgecolor="black", linewidth=0.6, zorder=5)
    if meth[names.index(name)]:
        return dict(color=PALETTE["methyl"], marker="^", s=65, edgecolor="black", linewidth=0.6, zorder=4)
    return dict(color=PALETTE["base"], marker="o", s=48, edgecolor="white", linewidth=0.5, zorder=3)


def annotate_selected(ax, xs, ys, labels_offsets):
    for label, (dx, dy) in labels_offsets.items():
        i = names.index(label)
        ax.annotate(
            label.replace("Benzo[", "B[").replace("Dibenzo[", "DB["),
            (xs[i], ys[i]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8,
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.85),
        )


def save(fig: plt.Figure, stem: str):
    fig.savefig(ROOT / f"{stem}.png", bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.pdf", bbox_inches="tight")


def fig1_main_results():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), facecolor=PALETTE["bg"])
    fig.patch.set_facecolor(PALETTE["bg"])

    # Panel A
    ax = axes[0]
    for i, name in enumerate(names):
        ax.scatter(K[i], logPEF[i], **point_style(name))
    rho_k, p_k = stats.spearmanr(K, logPEF)
    ax.set_title("A. Kekule Count Alone")
    ax.set_xlabel("K")
    ax.set_ylabel("log10(PEF)")
    ax.text(
        0.03,
        0.97,
        f"Spearman rho = {rho_k:+.3f}\np = {p_k:.2e}\nBonferroni: not significant",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc"),
    )
    annotate_selected(
        ax,
        K,
        logPEF,
        {
            "Coronene": (6, -6),
            "DMBA": (6, 4),
            "Dibenzo[a,e]pyrene": (6, -14),
            "Benzo[a]pyrene": (6, 4),
        },
    )

    # Panel B
    ax = axes[1]
    false_neg = {"5-Methylchrysene", "Benzo[k]fluoranthene"}
    false_pos = {"Benzo[ghi]perylene", "Benzo[e]pyrene"}
    for i, name in enumerate(names):
        style = point_style(name)
        if name in false_neg or name in false_pos:
            style["edgecolor"] = PALETTE["false"]
            style["linewidth"] = 1.6
            style["s"] = 76
        ax.scatter(KA[i], logPEF[i], **style)
    rho_ka, p_ka = stats.spearmanr(KA, logPEF)
    slope, intercept = np.polyfit(KA, logPEF, 1)
    xs = np.linspace(0, 18, 200)
    ax.plot(xs, slope * xs + intercept, color=PALETTE["line"], linewidth=1.4, alpha=0.8)
    ax.axvline(5.0, color="#777777", linestyle="--", linewidth=1.1)
    ax.axhline(-1.0, color="#999999", linestyle=":", linewidth=1.0)
    ax.set_title("B. Symmetry-Corrected Localization Index")
    ax.set_xlabel("K/|Aut|")
    ax.set_ylabel("log10(PEF)")
    ax.text(
        0.03,
        0.97,
        f"Spearman rho = {rho_ka:+.3f}\np = {p_ka:.2e}\nAUC = {roc_auc_score((PEF >= 0.1).astype(int), KA):.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc"),
    )
    annotate_selected(
        ax,
        KA,
        logPEF,
        {
            "5-Methylchrysene": (6, 4),
            "Benzo[k]fluoranthene": (6, 4),
            "Benzo[c]phenanthrene": (6, 4),
            "Benzo[ghi]perylene": (6, 4),
            "Coronene": (6, -12),
            "Dibenzo[a,e]pyrene": (6, 4),
        },
    )

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["base"], markeredgecolor="white", markersize=7, label="Unsubstituted PAH"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=PALETTE["methyl"], markeredgecolor="black", markersize=7, label="Methylated PAH"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=PALETTE["coronene"], markeredgecolor="black", markersize=7, label="Coronene"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=PALETTE["false"], markeredgewidth=1.6, markersize=8, label="Boundary misclassification"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, "fig1_main_results_v4")


def fig2_sensitivity_loo():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), facecolor=PALETTE["bg"])
    fig.patch.set_facecolor(PALETTE["bg"])

    # Panel A: sensitivity
    ax = axes[0]
    labels = ["Floor 0.01", "Floor 0.001", "Floor 0.0001", "BcP=0.001", "N&L only"]
    values = []
    for floor in [0.01, 0.001, 0.0001]:
        values.append(stats.spearmanr(KA, np.log10(np.maximum(PEF, floor)))[0])
    p_cons = PEF.copy()
    p_cons[names.index("Benzo[c]phenanthrene")] = 0.001
    values.append(stats.spearmanr(KA, np.log10(p_cons))[0])
    p_nl = np.array([audit.NL_VALUES.get(d[0], 0.001) for d in audit.DATA], dtype=float)
    values.append(stats.spearmanr(KA, np.log10(p_nl))[0])
    colors = [PALETTE["base"], PALETTE["base"], PALETTE["base"], "#457b9d", "#b56576"]
    ax.bar(range(len(labels)), values, color=colors, edgecolor="white")
    ax.set_ylim(0, 0.85)
    ax.set_ylabel("Spearman rho")
    ax.set_title("A. Sensitivity of the Main Correlation")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=28, ha="right")
    for i, v in enumerate(values):
        ax.text(i, v + 0.02, f"{v:.3f}", ha="center", va="bottom", fontsize=8)

    # Panel B: LOO
    ax = axes[1]
    predicted = np.zeros(len(KA))
    for i in range(len(KA)):
        x_train = np.delete(KA, i)
        y_train = np.delete(logPEF, i)
        coeff = np.polyfit(x_train, y_train, 1)
        predicted[i] = np.polyval(coeff, KA[i])
    for i, name in enumerate(names):
        ax.scatter(logPEF[i], predicted[i], **point_style(name))
    lims = [-3.2, 2.2]
    ax.plot(lims, lims, linestyle="--", color="#777777", linewidth=1.0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Observed log10(PEF)")
    ax.set_ylabel("LOO predicted log10(PEF)")
    ax.set_title("B. Leave-One-Out Predictions")
    rho_loo, p_loo = stats.spearmanr(predicted, logPEF)
    rmse = np.sqrt(np.mean((logPEF - predicted) ** 2))
    ax.text(
        0.03,
        0.97,
        f"rho = {rho_loo:+.3f}\np = {p_loo:.2e}\nRMSE = {rmse:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc"),
    )
    annotate_selected(
        ax,
        logPEF,
        predicted,
        {
            "DMBA": (6, 4),
            "Dibenzo[a,e]pyrene": (6, -14),
            "Dibenzo[a,h]pyrene": (6, 4),
        },
    )

    fig.tight_layout()
    save(fig, "fig2_sensitivity_loo_v4")


def fig4_coronene_bap():
    """
    Layout: 3-column
      col 0: Coronene structure + descriptor label block
      col 1: BaP structure + descriptor label block
      col 2: bar chart (K vs K/|Aut|) + interpretive caption
    Structure images are placed with tight spacing, descriptor text
    goes *below* each image to avoid title overlap.
    """
    fig = plt.figure(figsize=(13.0, 5.2), facecolor=PALETTE["bg"])
    fig.patch.set_facecolor(PALETTE["bg"])

    # 3 columns: mol images + text, mol images + text, bar chart
    gs = fig.add_gridspec(
        3, 3,
        width_ratios=[1.0, 1.0, 1.1],
        height_ratios=[0.12, 1.0, 0.55],
        hspace=0.08, wspace=0.10,
    )

    mol_info = [
        ("Coronene",        PALETTE["coronene"], 0),
        ("Benzo[a]pyrene",  PALETTE["base"],     1),
    ]

    for name, color, col in mol_info:
        i = names.index(name)

        # Title row
        ax_title = fig.add_subplot(gs[0, col])
        ax_title.axis("off")
        ax_title.text(
            0.5, 0.5, name,
            ha="center", va="center",
            fontsize=13, fontweight="bold", color=color,
            transform=ax_title.transAxes,
        )

        # Molecule image
        mol = Chem.MolFromSmiles(SMILES[name])
        img = Draw.MolToImage(mol, size=(520, 380))
        ax_img = fig.add_subplot(gs[1, col])
        ax_img.imshow(img)
        ax_img.axis("off")

        # Descriptor block below image
        ax_desc = fig.add_subplot(gs[2, col])
        ax_desc.axis("off")
        descriptor_lines = [
            f"K = {int(K[i])}    |Aut(G)| = {int(Aut[i])}",
            f"K/|Aut(G)| = {KA[i]:.2f}",
            f"PEF = {PEF[i]:.3g}   (IARC {'Group 1' if name == 'Benzo[a]pyrene' else 'Group 3'})",
        ]
        ax_desc.text(
            0.5, 0.80,
            descriptor_lines[0],
            ha="center", va="top",
            fontsize=11, color=PALETTE["line"],
            transform=ax_desc.transAxes,
        )
        ax_desc.text(
            0.5, 0.48,
            descriptor_lines[1],
            ha="center", va="top",
            fontsize=13, fontweight="bold", color=color,
            transform=ax_desc.transAxes,
        )
        ax_desc.text(
            0.5, 0.14,
            descriptor_lines[2],
            ha="center", va="top",
            fontsize=10, color="#555555",
            transform=ax_desc.transAxes,
        )

    # Right column: bar chart spanning all rows
    ax_bar = fig.add_subplot(gs[:, 2])
    labels  = ["Coronene\nK", "Coronene\nK/|Aut|", "BaP\nK", "BaP\nK/|Aut|"]
    values  = [20, 20 / 12, 9, 9.0]
    colors  = [PALETTE["coronene"], "#8ecae6", PALETTE["base"], "#90be6d"]
    bars = ax_bar.bar(range(4), values, color=colors, edgecolor="white", width=0.62)
    ax_bar.set_xticks(range(4))
    ax_bar.set_xticklabels(labels, fontsize=10)
    ax_bar.set_ylabel("Value", fontsize=11)
    ax_bar.set_title(
        "K counts possibilities;\nsymmetry decides which stay distinct",
        fontsize=11, fontweight="bold", pad=8,
    )
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    for idx, (bar, v) in enumerate(zip(bars, values)):
        label = f"{v:.2f}" if idx in {1} else f"{int(v):.0f}"
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2,
            v + 0.25,
            label,
            ha="center", va="bottom", fontsize=10, fontweight="bold",
        )
    # Threshold line
    ax_bar.axhline(5.0, color="#d62828", linestyle="--", linewidth=1.2, alpha=0.7)
    ax_bar.text(3.55, 5.15, "threshold 5.0", fontsize=8.5, color="#d62828", ha="right")
    ax_bar.set_ylim(0, 24)
    ax_bar.set_facecolor(PALETTE["bg"])

    save(fig, "fig4_coronene_bap_v4")


def fig3_hierarchy():
    """
    Clean left-to-right flowchart with 4 boxes in a single row.
    Box 4 (Bay/Fjord) sits below box 3, connected by a downward arrow.
    All text is manually wrapped to avoid matplotlib wrap=True artifacts.
    Font sizes: box title 12pt bold, box subtitle 10pt.
    """
    fig, ax = plt.subplots(figsize=(12.5, 5.2), facecolor=PALETTE["bg"])
    fig.patch.set_facecolor(PALETTE["bg"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Box definitions: (center_x, center_y, half_w, half_h, title, subtitle_lines, color)
    BOX_HW = 0.115   # half-width
    BOX_HH = 0.165   # half-height
    boxes = [
        (0.10, 0.68, BOX_HW, BOX_HH,
         "Even π-atom count",
         ["Structural precondition:", "Kekulé matching requires", "even-vertex graph"],
         "#2c6e49"),
        (0.36, 0.68, BOX_HW, BOX_HH,
         "|Aut(G)| ≤ 2",
         ["Symmetry gate:", "low symmetry permits", "bond inequivalence"],
         PALETTE["base"]),
        (0.62, 0.68, BOX_HW, BOX_HH,
         "K/|Aut(G)| ≥ 5.0",
         ["Localization index:", "matching multiplicity", "survives symmetry discount"],
         "#c77dff"),
        (0.62, 0.26, BOX_HW, BOX_HH,
         "Bay/Fjord region",
         ["Geometric gate:", "diol-epoxide pathway", "is structurally accessible"],
         "#e76f51"),
    ]

    for cx, cy, hw, hh, title, sub_lines, color in boxes:
        # Box background
        rect = FancyBboxPatch(
            (cx - hw, cy - hh), 2 * hw, 2 * hh,
            boxstyle="round,pad=0.015",
            facecolor="white",
            edgecolor=color,
            linewidth=2.0,
            transform=ax.transAxes,
            zorder=2,
        )
        ax.add_patch(rect)
        # Title
        ax.text(
            cx, cy + hh * 0.52,
            title,
            ha="center", va="center",
            fontsize=12, fontweight="bold", color=color,
            transform=ax.transAxes, zorder=3,
        )
        # Subtitle lines
        line_y = cy + hh * 0.05
        for line in sub_lines:
            ax.text(
                cx, line_y,
                line,
                ha="center", va="center",
                fontsize=10, color="#333333",
                transform=ax.transAxes, zorder=3,
            )
            line_y -= hh * 0.38

    # Arrows between boxes
    arrow_props = dict(arrowstyle="-|>", color="#555555", lw=1.8,
                       mutation_scale=16, connectionstyle="arc3,rad=0.0")
    def add_arrow(x0, y0, x1, y1):
        ax.annotate(
            "", xy=(x1, y1), xytext=(x0, y0),
            xycoords=ax.transAxes, textcoords=ax.transAxes,
            arrowprops=arrow_props, zorder=1,
        )

    # Box 1 → Box 2 (horizontal)
    add_arrow(0.10 + BOX_HW + 0.008, 0.68, 0.36 - BOX_HW - 0.008, 0.68)
    # Box 2 → Box 3 (horizontal)
    add_arrow(0.36 + BOX_HW + 0.008, 0.68, 0.62 - BOX_HW - 0.008, 0.68)
    # Box 3 → Box 4 (downward)
    add_arrow(0.62, 0.68 - BOX_HH - 0.008, 0.62, 0.26 + BOX_HH + 0.008)

    # Arrow labels
    ax.text(0.23, 0.715, "necessary", ha="center", va="bottom", fontsize=9,
            color="#555555", transform=ax.transAxes, style="italic")
    ax.text(0.49, 0.715, "necessary", ha="center", va="bottom", fontsize=9,
            color="#555555", transform=ax.transAxes, style="italic")
    ax.text(0.645, 0.47, "required", ha="left", va="center", fontsize=9,
            color="#555555", transform=ax.transAxes, style="italic")

    # Outcome labels
    outcome_props = dict(fontsize=11, fontweight="bold", transform=ax.transAxes, zorder=3)
    ax.text(0.88, 0.26, "→  CARCINOGENIC", color="#d62828", **outcome_props)
    ax.text(0.88, 0.68, "→  Not sufficient\n    alone", color="#777777",
            fontsize=10, transform=ax.transAxes, zorder=3, va="center")

    # Bottom caption
    caption = (
        "Low symmetry (|Aut| ≤ 2) is necessary but not sufficient.  "
        "K/|Aut(G)| ≥ 5.0 adds graded localization potential.  "
        "Bay/fjord geometry provides the mechanistic exit to the diol-epoxide pathway."
    )
    ax.text(
        0.02, 0.04, caption,
        ha="left", va="bottom", fontsize=9.5, color="#444444",
        transform=ax.transAxes, zorder=3,
    )

    save(fig, "fig3_hierarchy_v4")


def figS1_roc():
    y = (PEF >= 0.1).astype(int)
    fpr, tpr, _ = roc_curve(y, KA)
    auc = roc_auc_score(y, KA)
    fig, ax = plt.subplots(figsize=(4.6, 4.4), facecolor=PALETTE["bg"])
    fig.patch.set_facecolor(PALETTE["bg"])
    ax.plot(fpr, tpr, color=PALETTE["base"], linewidth=2.0, label=f"AUC = {auc:.2f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#999999")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Figure S1. ROC for K/|Aut|")
    ax.legend(frameon=False, loc="lower right")
    save(fig, "figS1_ROC_v4")


def main():
    fig1_main_results()
    fig2_sensitivity_loo()
    fig3_hierarchy()
    fig4_coronene_bap()
    figS1_roc()


if __name__ == "__main__":
    main()
