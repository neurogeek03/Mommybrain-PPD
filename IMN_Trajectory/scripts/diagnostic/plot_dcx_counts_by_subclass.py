"""
Stacked bar plots showing Dcx+ / Dcx- cell counts per sample for three subclasses,
with samples coloured by treatment (cort vs oil).

Data are hard-coded from docs/dcx_cells.md.

Usage:
    conda run -n anndata_env python scripts/diagnostic/plot_dcx_counts_by_subclass.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── data ───────────────────────────────────────────────────────────────────────
CORT_SAMPLES = {"BC28", "BC3", "BC9"}
OIL_SAMPLES  = {"BC15", "BC14", "BC13"}

# ordered: oil first, then cort (consistent across all panels)
SAMPLE_ORDER = ["BC13", "BC14", "BC15", "BC3", "BC9", "BC28"]

DATA = {
    "038 DG-PIR\nEx IMN": {
        "BC3":  {"Dcx+": 0,  "Dcx-": 2},
        "BC9":  {"Dcx+": 0,  "Dcx-": 3},
        "BC13": {"Dcx+": 1,  "Dcx-": 3},
        "BC14": {"Dcx+": 1,  "Dcx-": 6},
        "BC28": {"Dcx+": 0,  "Dcx-": 2},
        # BC15 absent from this subclass
    },
    "037 DG\nGlut": {
        "BC3":  {"Dcx+": 15, "Dcx-": 967},
        "BC9":  {"Dcx+": 20, "Dcx-": 1133},
        "BC13": {"Dcx+": 21, "Dcx-": 555},
        "BC14": {"Dcx+": 8,  "Dcx-": 1031},
        "BC15": {"Dcx+": 99, "Dcx-": 547},
        "BC28": {"Dcx+": 30, "Dcx-": 1645},
    },
    "045 OB-STR-CTX\nInh IMN": {
        "BC3":  {"Dcx+": 0,  "Dcx-": 57},
        "BC9":  {"Dcx+": 1,  "Dcx-": 18},
        "BC13": {"Dcx+": 2,  "Dcx-": 13},
        "BC14": {"Dcx+": 0,  "Dcx-": 30},
        "BC15": {"Dcx+": 1,  "Dcx-": 33},
        "BC28": {"Dcx+": 0,  "Dcx-": 16},
    },
}

# ── style ──────────────────────────────────────────────────────────────────────
COLOR_DCX_POS = "#C0392B"   # red
COLOR_DCX_NEG = "#AED6F1"   # light blue

COLOR_OIL   = "#E67E22"     # orange — used for tick labels
COLOR_CORT  = "#27AE60"     # green  — used for tick labels

OUT_DIR = "out"
os.makedirs(OUT_DIR, exist_ok=True)

# ── figure ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    2, 3,
    figsize=(13, 8),
    gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
)

subclass_keys = list(DATA.keys())

for col, subclass in enumerate(subclass_keys):
    ax_bar = axes[0, col]
    ax_pct = axes[1, col]

    samples = [s for s in SAMPLE_ORDER if s in DATA[subclass]]
    dcx_pos = [DATA[subclass][s]["Dcx+"] for s in samples]
    dcx_neg = [DATA[subclass][s]["Dcx-"] for s in samples]
    totals  = [p + n for p, n in zip(dcx_pos, dcx_neg)]
    pct_pos = [100 * p / t if t > 0 else 0 for p, t in zip(dcx_pos, totals)]

    x = np.arange(len(samples))
    bar_w = 0.55

    # stacked bars
    ax_bar.bar(x, dcx_neg, width=bar_w, color=COLOR_DCX_NEG, label="Dcx\u2212")
    ax_bar.bar(x, dcx_pos, width=bar_w, bottom=dcx_neg, color=COLOR_DCX_POS, label="Dcx+")

    # annotate Dcx+ count on top of each bar
    for xi, (pos, neg) in enumerate(zip(dcx_pos, dcx_neg)):
        if pos > 0:
            ax_bar.text(
                xi, neg + pos + 0.5, str(pos),
                ha="center", va="bottom", fontsize=7.5, color=COLOR_DCX_POS, fontweight="bold"
            )

    ax_bar.set_title(subclass, fontsize=10, fontweight="bold", pad=6)
    ax_bar.set_ylabel("Cell count", fontsize=9)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([])   # labels shown in bottom panel
    ax_bar.tick_params(axis="x", bottom=False)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    # % Dcx+ panel
    bar_colors = [COLOR_DCX_POS if p > 0 else "#EAECEE" for p in pct_pos]
    ax_pct.bar(x, pct_pos, width=bar_w, color=bar_colors, edgecolor="white")
    ax_pct.set_ylabel("%Dcx+", fontsize=9)
    ax_pct.set_xticks(x)
    ax_pct.set_xticklabels(samples, fontsize=9)
    ax_pct.spines["top"].set_visible(False)
    ax_pct.spines["right"].set_visible(False)

    # colour x-tick labels by treatment
    for tick_label, sample in zip(ax_pct.get_xticklabels(), samples):
        if sample in CORT_SAMPLES:
            tick_label.set_color(COLOR_CORT)
            tick_label.set_fontweight("bold")
        else:
            tick_label.set_color(COLOR_OIL)
            tick_label.set_fontweight("bold")

    # add a thin vertical separator between oil and cort groups
    # oil: first positions, cort: remaining
    n_oil = sum(1 for s in samples if s in OIL_SAMPLES)
    if 0 < n_oil < len(samples):
        for ax in [ax_bar, ax_pct]:
            ax.axvline(n_oil - 0.5, color="#CCCCCC", linewidth=1, linestyle="--", zorder=0)

    # treatment bracket labels at top of bar panel
    n_cort = len(samples) - n_oil
    if n_oil > 0:
        ax_bar.annotate(
            "Oil", xy=(n_oil / 2 - 0.5, 1.0), xycoords=("data", "axes fraction"),
            xytext=(0, 6), textcoords="offset points",
            ha="center", va="bottom", fontsize=8.5, color=COLOR_OIL, fontweight="bold"
        )
    if n_cort > 0:
        ax_bar.annotate(
            "Cort", xy=(n_oil + n_cort / 2 - 0.5, 1.0), xycoords=("data", "axes fraction"),
            xytext=(0, 6), textcoords="offset points",
            ha="center", va="bottom", fontsize=8.5, color=COLOR_CORT, fontweight="bold"
        )

# ── shared legend ──────────────────────────────────────────────────────────────
dcx_neg_patch = mpatches.Patch(color=COLOR_DCX_NEG, label="Dcx\u2212")
dcx_pos_patch = mpatches.Patch(color=COLOR_DCX_POS, label="Dcx+")
oil_patch  = mpatches.Patch(color=COLOR_OIL,  label="Oil (BC13, BC14, BC15)")
cort_patch = mpatches.Patch(color=COLOR_CORT, label="Cort (BC3, BC9, BC28)")

fig.legend(
    handles=[dcx_neg_patch, dcx_pos_patch, oil_patch, cort_patch],
    loc="lower center",
    ncol=4,
    fontsize=9,
    frameon=False,
    bbox_to_anchor=(0.5, -0.01),
)

fig.suptitle("Dcx+ vs Dcx\u2212 cells by subclass and treatment", fontsize=12, fontweight="bold", y=1.02)

out_path = os.path.join(OUT_DIR, "dcx_counts_by_subclass.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")
plt.close()
