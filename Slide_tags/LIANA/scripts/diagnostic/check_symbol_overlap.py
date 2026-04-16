"""
Diagnostic: check how well rat gene symbols from EdgeR DEA
match the mouseconsensus LR resource used in df_to_lr.

Key distinction: mouseconsensus uses underscore-joined names for receptor
complexes (e.g. Acvr1_Acvr2a). These are handled internally by LIANA by
splitting on '_' — they will never appear as single DEA row names. This
script separates simple genes from complexes before computing overlap.

Run from scripts/ with:
    conda run -n liana python check_symbol_overlap.py
"""

from pathlib import Path
import pandas as pd
import liana as li

project_path = Path.cwd().parents[0]

# ── 1. Load rat DEA gene symbols ────────────────────────────────────────────
dea_path = project_path / 'input' / 'edgeR_dge_input_liana.csv'
dea_df = pd.read_csv(dea_path, index_col=0)

rat_genes = set(dea_df.index.unique())
print(f"Unique rat gene symbols in EdgeR DEA: {len(rat_genes)}")

# ── 2. Load mouseconsensus and split simple genes from complexes ─────────────
resource = li.rs.select_resource('mouseconsensus')
mouse_lr_genes = set(resource['ligand'].unique()) | set(resource['receptor'].unique())

simple_lr  = {g for g in mouse_lr_genes if '_' not in g}
complex_lr = {g for g in mouse_lr_genes if '_' in g}

# Decompose complexes into their individual subunit symbols
complex_subunits = {sub for g in complex_lr for sub in g.split('_')}

print(f"\nmouseconsensus resource:")
print(f"  Total unique LR entries:    {len(mouse_lr_genes)}")
print(f"  Simple genes (no '_'):      {len(simple_lr)}")
print(f"  Complex entries (has '_'):  {len(complex_lr)}")
print(f"  Unique subunits in complexes: {len(complex_subunits)}")

# ── 3. Overlap for simple genes ─────────────────────────────────────────────
overlap_simple  = rat_genes & simple_lr
missing_simple  = simple_lr - rat_genes

print(f"\nSimple gene overlap (rat DEA ∩ mouseconsensus simple LR genes):")
print(f"  {len(overlap_simple)} / {len(simple_lr)} found  "
      f"({100 * len(overlap_simple) / len(simple_lr):.1f}%)")
print(f"  {len(missing_simple)} simple LR genes NOT in rat DEA")

# ── 4. Overlap for complex subunits ─────────────────────────────────────────
overlap_subunits = rat_genes & complex_subunits
missing_subunits = complex_subunits - rat_genes

print(f"\nComplex subunit overlap (rat DEA ∩ individual subunit symbols):")
print(f"  {len(overlap_subunits)} / {len(complex_subunits)} subunits found  "
      f"({100 * len(overlap_subunits) / len(complex_subunits):.1f}%)")
print(f"  {len(missing_subunits)} subunit symbols NOT in rat DEA")

# ── 5. Diagnose missing simple genes: naming issue or absent expression? ─────
# Use the RAW h5ad (all genes before HVG filtering), not the preprocessed one.
# The preprocessed h5ad only has ~2000 HVGs and gives false "absent" results.
import glob
import scanpy as sc

raw_h5ad_files = glob.glob(str(project_path / 'data' / '*.h5ad'))
h5ad_files     = glob.glob(str(project_path / 'out' / 'runs' / '*' / 'preprocessing' / '*.h5ad'))

if raw_h5ad_files:
    print(f"\nUsing raw h5ad for gene universe: {raw_h5ad_files[0]}")
    adata = sc.read_h5ad(raw_h5ad_files[0])
    adata_genes = set(adata.var_names)
    missing_but_in_adata = missing_simple & adata_genes
    missing_entirely     = missing_simple - adata_genes
    print(f"\nOf {len(missing_simple)} missing simple LR genes:")
    print(f"  {len(missing_but_in_adata)} ARE in adata.var_names "
          f"(filtered out of DEA due to low expression)")
    print(f"  {len(missing_entirely)} are NOT in adata at all "
          f"(genuine absence or rat/mouse name divergence)")
    print(f"\n  Absent from adata entirely (first 30 — potential naming divergence):")
    print(sorted(missing_entirely)[:30])
else:
    print(f"\nNo raw h5ad found in data/ — skipping adata cross-check.")
    print(f"Missing simple LR genes (first 30):")
    print(sorted(missing_simple)[:30])

# ── 6. Case-insensitive sanity check ────────────────────────────────────────
case_gains = len({g.upper() for g in rat_genes} & {g.upper() for g in simple_lr}) - len(overlap_simple)
print(f"\nCase-insensitive check on simple genes: +{case_gains} additional matches "
      f"(0 = no case mismatch between rat and mouse conventions)")

# ── 7. CCL/CCR chemokine axis check ─────────────────────────────────────────
# These receptors are relevant for CORT-driven neuroinflammation / microglial CCC.
# They have identical names in mouse and rat, so absence from adata suggests
# expression filtering (HVG) rather than a naming issue.
print("\n" + "="*60)
print("CCL/CCR chemokine axis (CORT neuroinflammation relevance)")
print("="*60)

chemokine_genes = [
    'Cx3cr1', 'Cx3cl1',   # fractalkine axis — neuron → microglia
    'Ccr2',   'Ccl2',     # monocyte/microglia recruitment
    'Ccr5',   'Ccl5',     # neuroinflammation
    'Ccr1',   'Ccl4',     # broad neuroinflammation
    'Ccl21',  'Ccl27',    # rat equivalents of mouse Ccl21b / Ccl27b
]

# Load HVG-filtered adata for the third column (if available)
hvg_genes = set()
if h5ad_files:
    adata_hvg = sc.read_h5ad(h5ad_files[0])
    hvg_genes = set(adata_hvg.var_names)

if raw_h5ad_files:
    print(f"{'Gene':<12} {'in raw':>8} {'in HVG':>8} {'in DEA':>8}  note")
    print("-" * 65)
    for g in chemokine_genes:
        in_raw = g in adata_genes
        in_hvg = g in hvg_genes
        in_dea = g in rat_genes
        if not in_raw:
            note = "genuinely absent from dataset — not a naming issue if in DEA"
        elif not in_dea:
            note = "expressed but below DEA threshold (won't appear in lr_res)"
        elif not in_hvg:
            note = "in DEA but cut by HVG selection — expr_prop still computed from raw"
        else:
            note = "fully covered"
        print(f"{g:<12} {str(in_raw):>8} {str(in_hvg):>8} {str(in_dea):>8}  {note}")
else:
    print("No raw h5ad found — checking DEA only.")
    print(f"{'Gene':<12} {'in DEA':>10}  note")
    print("-" * 40)
    for g in chemokine_genes:
        in_dea = g in rat_genes
        print(f"{g:<12} {str(in_dea):>10}")
