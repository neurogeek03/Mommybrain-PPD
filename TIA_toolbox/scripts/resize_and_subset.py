import argparse
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc

from skimage.transform import resize

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# =========== Arg parser ===========
parser = argparse.ArgumentParser(
        description="Resize mask array and subset adata object"
    )

parser.add_argument("-s", "--sample", required=True, help="Sample name (e.g. B01)")
parser.add_argument("-i", "--input", required=True, help="Path to adata object")
parser.add_argument("-m", "--masks", required=True, help="Path to mask folder")
parser.add_argument("-o", "--output", required=True, help="Path to mask folder")

args = parser.parse_args()
sample = args.sample
adata_path = Path(args.input)
mask_dir = Path(args.masks)
output_base = Path(args.output)

# =========== Paths ===========
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
# output_base = project_path / 'out'
# output_base.mkdir(exist_ok=True, parents=True)

# =========== adata ===========
# obj_path = input_dir / f'{sample}_anndata.h5ad'
adata = sc.read_h5ad(adata_path)

# mask
mask_array_path = mask_dir / f'{sample}_mask.npy'
mask_array = np.load(mask_array_path)
print(mask_array)
print(mask_array.shape)

print('resizing array...')
arr_big = resize(mask_array,(7557, 7602),order=0).astype(mask_array.dtype)
print(arr_big.size)

print(arr_big)
print(arr_big.shape)

# remove padding from mask, in the same way as we added it to the input image
unpadded = arr_big[500:-500, 500:-500]
print(unpadded)
print(unpadded.shape)
print(unpadded.size)

# =================== Build barcode matrix ===================
# coords & convert to integers 
coords = pd.DataFrame(adata.obsm["X_spatial"], columns=["x", "y"])
x = coords["x"].astype(int).to_numpy()
y = coords["y"].astype(int).to_numpy()
barcodes = adata.obs_names.to_numpy()

print(f'the shape of the coords df is {coords.shape}')

# Convert to arrays
x = np.asarray(x)
y = np.asarray(y)
b = np.asarray(barcodes)

# Get sorted unique coordinates
xs = np.sort(np.unique(x))
ys = np.sort(np.unique(y))

x_to_col = {v: i for i, v in enumerate(xs)}
y_to_row = {v: i for i, v in enumerate(ys)}

M = np.empty((len(ys), len(xs)), dtype=object)

for xc, yc, bc in zip(x, y, b):
    M[y_to_row[yc], x_to_col[xc]] = bc

rows, cols = np.where(M != None)
plt.scatter(xs[cols], ys[rows])

# building a color dict to sense check array
umis = adata.obs["log10_numReads"].to_numpy() 
umi_map = dict(zip(barcodes, umis))
bead_barcodes = M[rows, cols]            # extract barcodes in matrix order
bead_colors = np.array([umi_map[bc] for bc in bead_barcodes])
# color_grid = np.empty_like(M, dtype=float)

# for r in range(M.shape[0]):
#     for c in range(M.shape[1]):
#         bc = M[r, c]
#         color_grid[r, c] = umi_map.get(bc, np.nan)   # missing spots become nan

# =================== PLOT barcodes array ===================
# plt.scatter(xs[cols], ys[rows], c=bead_colors, s=3, cmap="viridis")
# plt.colorbar()
# plt.xlim(xs.min(), xs.max())
# plt.ylim(ys.min(), ys.max())
# plt.gca().invert_yaxis()
# path_array_plot = output_base / 'test_array_plot.png' #works!!
# plt.savefig(path_array_plot)

# =================== Resize mask to fit array ===================
# Resize with nearest-neighbor to preserve mask
mask_resized = resize(
    unpadded,
    (len(ys), len(xs)),
    order=0,
    preserve_range=True
).astype(bool)
print('shapes of 1)mask 2)our barcode array')
print(mask_resized.shape, M.shape)

mask_resized = np.flipud(mask_resized)

keep = mask_resized[rows, cols] # boolean array

filtered_barcodes = bead_barcodes[keep]
filtered_rows = rows[keep]
filtered_cols = cols[keep]

filtered_x = xs[filtered_cols]
filtered_y = ys[filtered_rows]
filtered_barcodes = bead_barcodes[keep]
print(f"Kept {len(filtered_barcodes)} / {adata.n_obs} beads inside tissue.")

# # =================== Plot to check ===================
# plt.figure(figsize=(9, 9))

# # All beads (background)
# plt.scatter(xs[cols], ys[rows], s=1, c="lightgray", label="All beads")

# # Mask-kept beads
# plt.scatter(xs[cols[keep]], ys[rows[keep]], 
#             s=1, c="red", label="Inside mask")

# plt.gca().invert_yaxis()    # spatial conventions
# plt.legend()
# plt.title("Mask sanity check on barcode grid")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.tight_layout()
# filtered_puck_path = output_base / 'filtered_puck.png'
# plt.savefig(filtered_puck_path)

# =================== Subset adata ===================
adata_filtered = adata[filtered_barcodes].copy()
adata_filtered.obsm["X_spatial"] = adata_filtered.obsm["X_spatial"].copy()

# extract plot data 
x = adata_filtered.obsm["X_spatial"][:, 0]
y = adata_filtered.obsm["X_spatial"][:, 1]
colors = adata_filtered.obs["log10_numReads"]

# 3. Plot using matplotlib
plt.figure(figsize=(8,6))
plt.scatter(x, y, c=colors, cmap="viridis", s=1)
plt.gca().invert_yaxis()  # optional, usually needed for spatial plots
plt.colorbar(label="log10_numReads")
plt.title("Filtered beads")

# 4. Save figure
fig_path = output_base / f"{sample}_{adata_filtered.n_obs}_masked_spatial.png"
plt.savefig(fig_path, dpi=300)

# 5. Save adata filtered
adata_fil_path = output_base / f'{sample}_{adata_filtered.n_obs}_masked.h5ad'
adata_filtered.write_h5ad(adata_fil_path)
