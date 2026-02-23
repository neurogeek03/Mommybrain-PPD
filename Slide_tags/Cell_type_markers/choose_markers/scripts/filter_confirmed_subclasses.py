"""
filter_confirmed_subclasses.py

Reads final_markers.csv and keeps only subclasses where at least one
marker gene was confirmed in both platforms (gene_1 is not empty).

Output: confirmed_marker_subclasses.csv — used as input for dotplots.
"""

import pandas as pd

IN_CSV  = "final_markers.csv"
OUT_CSV = "confirmed_marker_subclasses.csv"

df = pd.read_csv(IN_CSV)

# Keep rows where at least gene_1 was assigned
confirmed = df[df["gene_1"].notna() & (df["gene_1"] != "")].copy()

# Collect all non-empty genes into a single list column for easy downstream use
def collect_genes(row):
    return [g for g in [row["gene_1"], row["gene_2"], row["gene_3"]]
            if pd.notna(g) and g != ""]

confirmed["markers"] = confirmed.apply(collect_genes, axis=1)
confirmed["n_markers"] = confirmed["markers"].apply(len)

# Drop subclasses with fewer than 2 markers
enough = confirmed[confirmed["n_markers"] >= 2].copy()

enough.to_csv(OUT_CSV, index=False)

print(f"Total subclasses in final_markers.csv:    {len(df)}")
print(f"With at least 1 gene confirmed:           {len(confirmed)}")
print(f"With at least 2 genes confirmed (saved):  {len(enough)}\n")

print(enough[["subclass_name", "broad_class", "gene_1", "gene_2", "gene_3"]].to_string(index=False))

# Also report excluded subclasses
excluded = df[df["gene_1"].isna() | (df["gene_1"] == "")]
print(f"\nExcluded ({len(excluded)}) — no shared markers found:")
for _, r in excluded.iterrows():
    print(f"  {r['subclass_name']}: {r['notes']}")
