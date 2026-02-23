"""
choose_markers.py

Loads both AnnData objects, runs sc.tl.rank_genes_groups on each (restricted
to shared subclasses), then selects 2-3 final marker genes per subclass using:

  Priority 1 — marker_1 from common_markers.csv in top-N of BOTH platforms
  Priority 2 — marker_2, then marker_3, same check
  Priority 3 — fallback: top genes ranked in BOTH platforms (by mean rank)

Outputs final_markers.csv.

NOTE: adata.X must be log-normalised before running this script.
"""

import re
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

sc.settings.verbosity = 1

# ============ PARAMETERS ============
project_path = Path(__file__).resolve().parents[1]
SLIDE_TAGS_H5AD  = "/scratch/mfafouti/Mommybrain/Slide_tags/Cell_type_markers/subclasses/data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
SLIDE_SEQ_H5AD   = "/scratch/mfafouti/Mommybrain/Slide_tags/Cell_type_markers/subclasses/data/filtered_combined_20260222_173834.h5ad"
SLIDE_TAGS_COL   = "subclass_name"
SLIDE_SEQ_COL    = "RCTD_first_type_rat"

COMMON_CSV       = project_path / "common_markers.csv"
OUT_CSV          = project_path / "updated_final_markers.csv"

TOP_N  = 50   # genes pulled from rank_genes_groups per subclass
N_GENES = 3   # target markers per subclass

# ============ HELPERS ============
def normalize(name):
    """Spaces, hyphens, slashes → underscores for name matching."""
    return re.sub(r"[ /-]", "_", str(name))

def extract_ranked_genes(adata, top_n):
    """
    Returns dict: normalized_subclass_key → list of top-N gene names (ranked).
    Assumes rank_genes_groups has already been run.
    """
    result = {}
    rgg = adata.uns["rank_genes_groups"]
    groups = rgg["names"].dtype.names
    for grp in groups:
        genes = list(rgg["names"][grp][:top_n])
        result[normalize(grp)] = genes
    return result

# ============ LOAD & SUBSET TO SHARED SUBCLASSES ============
common = pd.read_csv(COMMON_CSV)
shared_names_st = set(common["subclass_name"])                       # spaces/hyphens
shared_names_ss = {normalize(n) for n in common["subclass_name"]}   # underscores

print("Loading slide-tags...")
st = sc.read_h5ad(SLIDE_TAGS_H5AD)
st = st[st.obs[SLIDE_TAGS_COL].isin(shared_names_st)].copy()
st.var_names = st.var["gene_symbol"].astype(str).values
st.var_names_make_unique()
print(f"  {st.n_obs:,} cells across {st.obs[SLIDE_TAGS_COL].nunique()} subclasses")

print("Loading slide-seq...")
ss = sc.read_h5ad(SLIDE_SEQ_H5AD)
# slide-seq uses underscores — match against normalised shared names
ss = ss[ss.obs[SLIDE_SEQ_COL].apply(normalize).isin(shared_names_ss)].copy()
ss.var_names = ss.var["name"].astype(str).values
ss.var_names_make_unique()
print(f"  {ss.n_obs:,} cells across {ss.obs[SLIDE_SEQ_COL].nunique()} subclasses")

# ============ NORMALIZE & LOG-TRANSFORM ============
print("\nNormalizing and log-transforming slide-tags...")
sc.pp.normalize_total(st)
sc.pp.log1p(st)
print("Preprocessing complete.")

print("\nNormalizing and log-transforming slide-seq...")
sc.pp.normalize_total(ss)
sc.pp.log1p(ss)
print("Preprocessing complete.")

# ============ RANK GENES GROUPS ============
print("\nRunning rank_genes_groups on slide-tags...")
sc.tl.rank_genes_groups(
    st,
    groupby=SLIDE_TAGS_COL,
    method="wilcoxon",
    n_genes=TOP_N,
    key_added="rank_genes_groups",
)

print("Running rank_genes_groups on slide-seq...")
sc.tl.rank_genes_groups(
    ss,
    groupby=SLIDE_SEQ_COL,
    method="wilcoxon",
    n_genes=TOP_N,
    key_added="rank_genes_groups",
)

# ============ EXTRACT RANKED GENE LISTS ============
# dict: normalised_subclass_key → [gene1, gene2, ..., geneN]  (rank order)
st_ranked = extract_ranked_genes(st, TOP_N)
ss_ranked = extract_ranked_genes(ss, TOP_N)

def rank_in(gene, ranked_list):
    """0-based rank of gene in list, or None if absent."""
    try:
        return ranked_list.index(gene)
    except ValueError:
        return None

def in_both_top_n(gene, key):
    """True if gene is in top-N for this subclass in both platforms."""
    st_r = rank_in(gene, st_ranked.get(key, []))
    ss_r = rank_in(gene, ss_ranked.get(key, []))
    return st_r is not None and ss_r is not None

def mean_rank(gene, key):
    """Mean rank across both platforms (lower = better). None if absent in either."""
    st_r = rank_in(gene, st_ranked.get(key, []))
    ss_r = rank_in(gene, ss_ranked.get(key, []))
    if st_r is None or ss_r is None:
        return None
    return (st_r + ss_r) / 2

# ============ SELECT MARKERS ============
rows = []

for _, row in common.iterrows():
    subclass = row["subclass_name"]
    key      = normalize(subclass)
    m1, m2, m3 = row["marker_1"], row["marker_2"], row["marker_3"]

    selected   = []   # (gene, reason)
    skip_notes = []

    for priority, gene in enumerate([m1, m2, m3], start=1):
        label = f"marker_{priority}"
        if pd.isna(gene) or gene == "":
            continue
        if in_both_top_n(gene, key):
            selected.append((gene, f"{label}: in top-{TOP_N} of both platforms"))
        else:
            st_r = rank_in(gene, st_ranked.get(key, []))
            ss_r = rank_in(gene, ss_ranked.get(key, []))
            where = []
            if st_r is not None: where.append("slide-tags")
            if ss_r is not None: where.append("slide-seq")
            loc = " & ".join(where) if where else "neither"
            skip_notes.append(f"{label} ({gene}): only in {loc} — skipped")

    # Fallback: genes in top-N of BOTH, ranked by mean rank
    if len(selected) < N_GENES:
        already   = {g for g, _ in selected}
        common_top = set(st_ranked.get(key, [])) & set(ss_ranked.get(key, []))
        candidates = common_top - already - {m1, m2, m3}
        ranked_candidates = sorted(
            candidates,
            key=lambda g: mean_rank(g, key)
        )
        for g in ranked_candidates:
            if len(selected) >= N_GENES:
                break
            selected.append((g, f"fallback: rank {mean_rank(g, key):.1f} mean across both"))

    # Build row
    genes   = [g for g, _ in selected[:N_GENES]]
    reasons = [n for _, n in selected[:N_GENES]] + skip_notes

    rows.append({
        "subclass_name": subclass,
        "broad_class":   row["broad_class"],
        "gene_1": genes[0] if len(genes) > 0 else "",
        "gene_2": genes[1] if len(genes) > 1 else "",
        "gene_3": genes[2] if len(genes) > 2 else "",
        "notes":  " | ".join(reasons),
    })

# ============ SAVE ============
out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)
print(f"\nSaved {len(out)} subclasses to: {OUT_CSV}")

# Summary
confirmed = out[out["gene_1"] != ""]
print(f"\nSubclasses with ≥1 marker confirmed in both: {len(confirmed)}")
print(f"Subclasses with no shared markers:           {len(out) - len(confirmed)}\n")

no_m1 = out[~out["notes"].str.contains(f"marker_1: in top-{TOP_N}")]
print(f"Subclasses where marker_1 was NOT in top-{TOP_N} of both ({len(no_m1)}):")
for _, r in no_m1.iterrows():
    print(f"  {r['subclass_name']}: {r['notes']}")
