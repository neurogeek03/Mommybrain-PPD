import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import seaborn as sns

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
data_dir = project_path / 'data' / 'objects'
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# =================== PARAMS ===================
files = sorted(data_dir.glob("*.h5ad"))
umi_data = []
labels = []

# =================== INPUT ===================
for file in files:
    sample_name = file.name.split("_")[0]

    adata = sc.read_h5ad(file)

    # Total UMIs per cell
    umi_per_cell = np.array(adata.X.sum(axis=1)).ravel()

    umi_data.append(umi_per_cell)
    labels.append(sample_name)

# =================== OUTPUT ===================
plt.figure(figsize=(9, 6))

# plt.violinplot(
#     umi_data,
#     showmeans=True,
#     showmedians=True,
#     showextrema=True
# )

plt.boxplot(
    umi_data,
    showfliers=False,
    medianprops=dict(color="black", linewidth=2)
)

plt.xticks(
    range(1, len(labels) + 1),
    labels,
    rotation=45,
    ha="right"
)

plt.xticks(range(1, len(labels) + 1), labels, rotation=45, ha="right")
plt.yscale("log")
plt.ylabel("log UMIs per cell")
plt.title("UMI per-cell distribution across samples")

plt.tight_layout()
fig_path = output_base / f'violin_{stamp}_umi_dist_allsamples.png'
plt.savefig(fig_path, dpi=300)

# COLLECT DATA
# -----------------------------
rows = []

for file in files:
    sample_name = file.name.split("_")[0]

    adata = sc.read_h5ad(file)

    umi_per_cell = np.array(adata.X.sum(axis=1)).ravel()

    rows.append(
        pd.DataFrame({
            "sample": sample_name,
            "umis_per_cell": umi_per_cell
        })
    )

# -----------------------------
# WRITE CSV
# -----------------------------
df = pd.concat(rows, ignore_index=True)
output_csv = output_base / "umi_per_cell_by_sample.csv"
df.to_csv(output_csv, index=False)

