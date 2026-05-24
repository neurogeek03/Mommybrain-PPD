#!/usr/bin/env python3
"""
Generate per-cell-type DEG summary reports for the top 10 cell types (lib_size_50k).
Output: out/runs/lib_size_50k/docs/top10_celltype_DEG_reports.md
"""

import pandas as pd
from pathlib import Path

RESULTS_DIR = Path("/scratch/mfafouti/Mommybrain-PPD/Slide_tags/EdgeR/out/runs/lib_size_50k/edger_lrt/OIL_vs_CORT")
OUT_FILE = Path("/scratch/mfafouti/Mommybrain-PPD/Slide_tags/EdgeR/out/runs/lib_size_50k/docs/top10_celltype_DEG_reports.md")
FDR_THRESH = 0.05
N_GENES = 15  # genes to show per table

# Top 10 cell types by DEG count (FDR < 0.05), ranked
TOP10 = [
    ("327_Oligo_NN",               "Oligodendrocytes",                608),
    ("333_Endo_NN",                "Endothelium",                     376),
    ("318_Astro-NT_NN",            "Astrocytes (NT subtype)",         313),
    ("319_Astro-TE_NN",            "Astrocytes (TE subtype)",         266),
    ("030_L6_CT_CTX_Glut",         "L6 Corticothalamic Neurons",      170),
    ("326_OPC_NN",                 "Oligodendrocyte Precursors (OPC)", 168),
    ("004_L6_IT_CTX_Glut",         "L6 Intratelencephalic Neurons",   168),
    ("151_TH_Prkcd_Grin2c_Glut",   "Thalamus (Prkcd/Grin2c+)",       127),
    ("007_L2_3_IT_CTX_Glut",       "L2/3 Intratelencephalic Neurons", 125),
    ("330_VLMC_NN",                "Vascular Leptomeningeal Cells",    111),
]


def fmt_fdr(val):
    if val < 1e-10:
        return f"{val:.2e}"
    elif val < 0.001:
        return f"{val:.3e}"
    else:
        return f"{val:.4f}"


def fmt_fc(val):
    return f"{val:+.2f}"


def make_section(cell_id, cell_label, n_degs_expected):
    tsv = RESULTS_DIR / f"{cell_id}_edgeR_results.tsv"
    df = pd.read_csv(tsv, sep="\t", index_col=0)
    df.index.name = "gene"
    df = df.reset_index()

    total_tested = len(df)
    sig = df[df["FDR"] < FDR_THRESH].copy()
    n_sig = len(sig)
    n_up = (sig["logFC"] > 0).sum()
    n_down = (sig["logFC"] < 0).sum()

    # Sort subsets
    by_sig = sig.sort_values("FDR").head(N_GENES)
    by_up = sig[sig["logFC"] > 0].sort_values("logFC", ascending=False).head(N_GENES)
    by_down = sig[sig["logFC"] < 0].sort_values("logFC").head(N_GENES)

    lines = []
    cell_id_fmt = cell_id.replace("_", "\\_")
    lines.append(f"## {cell_id_fmt} — {cell_label}")
    lines.append("")
    lines.append(f"**Genes tested:** {total_tested:,} | "
                 f"**Significant (FDR < 0.05):** {n_sig} | "
                 f"**Up in CORT:** {n_up} | **Down in CORT:** {n_down}")
    lines.append("")

    def gene_table(subset, sort_col, title):
        lines.append(f"### {title}")
        lines.append("")
        lines.append("| Gene | logFC | logCPM | FDR |")
        lines.append("|---|---|---|---|")
        for _, row in subset.iterrows():
            lines.append(
                f"| {row['gene']} | {fmt_fc(row['logFC'])} | {row['logCPM']:.2f} | {fmt_fdr(row['FDR'])} |"
            )
        lines.append("")

    gene_table(by_sig,  "FDR",   f"Top {N_GENES} Most Significant (by FDR)")
    gene_table(by_up,   "logFC", f"Top {N_GENES} Most Up-regulated in CORT (FDR < 0.05)")
    gene_table(by_down, "logFC", f"Top {N_GENES} Most Down-regulated in CORT (FDR < 0.05)")

    lines.append("---")
    lines.append("")
    return "\n".join(lines)


header = """# Top 10 Cell Types — DEG Reports (lib_size_50k)

**Comparison:** OIL (control) vs CORT (corticosterone-treated)
**Model:** EdgeR LRT | Library size filter: 50k
**Significance threshold:** FDR < 0.05
**Date:** 2026-04-28

Positive logFC = higher in CORT. Negative logFC = higher in OIL/control.
Cell types ranked by total number of significant DEGs (descending).

---

"""

sections = [header]
for cell_id, label, n in TOP10:
    print(f"Processing {cell_id}...")
    sections.append(make_section(cell_id, label, n))

OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
OUT_FILE.write_text("\n".join(sections))
print(f"\nDone. Written to {OUT_FILE}")
