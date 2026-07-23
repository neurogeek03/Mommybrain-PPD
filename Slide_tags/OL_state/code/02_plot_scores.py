"""
Combined grid figure: all UCell signatures × cell types in one PNG.

Rows = signatures (grouped by section A–E from docs/gene_signatures.md).
Columns = cell types.
Each subplot is a CORT-vs-OIL KDE on a shared x-axis per row.

Run with:
    uv run python code/02_plot_scores.py
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

CSV = "data/ucell_scores_v4.csv"
OUT_PATH = "data/ucell_combined_grid.png"

CELL_TYPE_COL = "subclass_name"
TREATMENT_COL = "treatment"

CELL_TYPE_ORDER = [
    "327 Oligo NN",
    "326 OPC NN",
    "318 Astro-NT NN",
    "319 Astro-TE NN",
    "333 Endo NN",
    "037 DG Glut",
    "052 Pvalb Gaba",
]

CELL_TYPE_LABELS = {
    "327 Oligo NN":     "Oligodendrocyte",
    "326 OPC NN":       "OPC",
    "318 Astro-NT NN":  "Astro-NT",
    "319 Astro-TE NN":  "Astro-TE",
    "333 Endo NN":      "Endothelial",
    "037 DG Glut":      "DG granule",
    "052 Pvalb Gaba":   "Pvalb int.",
}

COLORS = {"CORT": "#d62728", "OIL": "#1f77b4"}

# Signatures grouped by section. Each tuple: (score_col, short_label, expected_direction)
SIGNATURE_ROWS = [
    # A. Stress / GR
    ("GR_activation_UCell",              "A1 GR activation",        "UP"),
    ("GR_brake_UCell",                   "A2 GR anti-inflam brake", "UP"),
    ("ISR_hypoxia_UCell",                "A3 ISR / hypoxia",        "UP"),
    # B. Metabolic
    ("Metabolic_reprogramming_UCell",    "B1 Metabolic reprog",     "UP"),
    ("Lipid_anabolism_loss_UCell",       "B1 Lipid anabolism",      "DOWN"),
    ("OXPHOS_UCell",                     "B2 OXPHOS (Cx I–V)",      "DOWN"),
    # C. Identity
    ("Myelin_identity_UCell",            "C1 Mature myelin",        "DOWN"),
    ("OPC_proliferation_UCell",          "C2 OPC proliferation",    "DOWN"),
    ("Astro_pan_reactive_UCell",         "C3 Astro pan-reactive",   "UP"),
    ("Astro_A1_neurotoxic_UCell",        "C3 Astro A1 neurotoxic",  "UP"),
    ("Astro_A2_protective_UCell",        "C3 Astro A2 protective",  "UP"),
    ("BBB_transport_UCell",              "C4 BBB transport",        "UP"),
    ("BBB_vascular_identity_UCell",      "C4 BBB vascular id",      "DOWN"),
    # D. Synapse
    ("Astro_glutamate_support_UCell",    "D1 Astro glutamate",      "DOWN"),
    ("Astro_GABA_support_UCell",         "D2 Astro GABA (PPD)",     "DOWN"),
    # E. Cross-cutting
    ("Glial_cyto_remodel_UCell",         "E1 Cyto remodeling",      "UP"),
    ("Clock_genes_UCell",                "E2 Circadian clock",      "DOWN"),
    ("Adenosine_activity_sensing_UCell", "E3 Adenosine sensing",    "UP"),
]


def main():
    df = pd.read_csv(CSV, index_col="barcode")

    n_rows = len(SIGNATURE_ROWS)
    n_cols = len(CELL_TYPE_ORDER)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(2.6 * n_cols, 1.5 * n_rows),
        sharex="row",
    )

    for i, (score_col, label, direction) in enumerate(SIGNATURE_ROWS):
        # Per-row x-limit from observed range (drop top 0.5% for legibility)
        row_vals = df[score_col].dropna()
        xmax = row_vals.quantile(0.995) if len(row_vals) else 1.0

        for j, ct in enumerate(CELL_TYPE_ORDER):
            ax = axes[i, j]
            sub = df[df[CELL_TYPE_COL] == ct]
            for treatment in ("OIL", "CORT"):
                vals = sub.loc[sub[TREATMENT_COL] == treatment, score_col].dropna()
                if len(vals) == 0:
                    continue
                sns.kdeplot(
                    vals,
                    ax=ax,
                    color=COLORS[treatment],
                    fill=True,
                    alpha=0.3,
                    linewidth=1.2,
                    clip=(0, 1),
                )

            # Column titles only on top row
            if i == 0:
                ax.set_title(CELL_TYPE_LABELS[ct], fontsize=11, fontweight="bold")

            # Row labels only on left column
            if j == 0:
                color = "#d62728" if direction == "UP" else "#1f77b4"
                ax.set_ylabel(
                    f"{label}\n({direction})",
                    fontsize=9, fontweight="bold", color=color,
                    rotation=0, ha="right", va="center", labelpad=8,
                )
            else:
                ax.set_ylabel("")

            ax.set_xlim(0, xmax)
            ax.set_xlabel("")
            ax.set_yticks([])
            ax.tick_params(axis="x", labelsize=7)
            sns.despine(ax=ax, left=True)

    # Single shared legend in upper right
    handles = [
        plt.Line2D([0], [0], color=COLORS["OIL"], lw=4, alpha=0.5, label="OIL"),
        plt.Line2D([0], [0], color=COLORS["CORT"], lw=4, alpha=0.5, label="CORT"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
               fontsize=11, frameon=True)

    fig.suptitle(
        f"UCell signature scores — CORT vs OIL across {n_cols} cell types\n"
        f"(rows: {n_rows} literature-curated signatures; expected direction shown in red=UP / blue=DOWN)",
        fontsize=12, y=1.002,
    )
    fig.tight_layout(rect=[0.05, 0, 0.98, 0.99])
    fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT_PATH}")


if __name__ == "__main__":
    main()
