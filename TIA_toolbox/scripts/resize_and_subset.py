import argparse
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc

from skimage.transform import resize
import matplotlib.colors as mcolors
from PIL import Image

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
# create new metadata column w/ mask
# print("adata.n_obs:", adata.n_obs)
# sanity check
print("len(bead_barcodes):", len(bead_barcodes))
print("len(rows):", len(rows))
print("len(cols):", len(cols))
print("len(keep):", len(keep))
print('---------------------------------------')
missing = adata.obs_names.difference(bead_barcodes)
print('Below is the number of beads that do not have spatial info:')
print(len(missing))

# initialize everyone as false as some bead barcodes do not have spatial coords 
print('===================================')
# adata.obs["in_tissue_mask"] = False
# mask_bool = np.isin(adata.obs_names, filtered_barcodes)
# adata.obs['in_tissue_mask'] = mask_bool
# adata.obs.loc[bead_barcodes, "in_tissue_mask"] = keep
# sanity check: these should match
assert set(filtered_barcodes).issubset(set(adata.obs_names))

# initialize: everyone is outside tissue
adata.obs["in_tissue"] = False

# mark kept barcodes as inside tissue
adata.obs.loc[filtered_barcodes, "in_tissue"] = True
print(adata.obs["in_tissue"].value_counts())


print('Preview of adata file to saved:')
print(adata)
print('===================================')
# adata_filtered = adata[adata.obs["in_tissue_mask"]].copy()

# # Side-by side plot of unmasked vs masked 
# adatas = {
#     "All beads": adata,
#     "Filtered beads": adata_filtered,
# }

# fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

# for ax, (title, ad) in zip(axes, adatas.items()):
#     x = ad.obsm["X_spatial"][:, 0]
#     y = ad.obsm["X_spatial"][:, 1]
#     colors = ad.obs["log10_numReads"]

#     sc = ax.scatter(
#         x, y,
#         c=colors,
#         cmap="viridis",
#         s=1
#     )
#     fig.colorbar(sc, ax=ax, label="log10_numReads")

#     ax.invert_yaxis()
#     ax.set_title(f"{title}\n(n = {ad.n_obs})")
#     ax.set_aspect("equal")

# fig.suptitle(sample)
# # fig.tight_layout()

# fig_path = output_base / f"{sample}_spatial_before_after_mask.png"
# plt.savefig(fig_path, dpi=300)


# # ======================


# extract plot data 
# x = adata_filtered.obsm["X_spatial"][:, 0]
# y = adata_filtered.obsm["X_spatial"][:, 1]
# colors = adata_filtered.obs["log10_numReads"]

# # 3. Plot using matplotlib
# plt.figure(figsize=(8,6))
# plt.scatter(x, y, c=colors, cmap="viridis", s=1)
# plt.gca().invert_yaxis()  # optional, usually needed for spatial plots
# plt.colorbar(label="log10_numReads")
# plt.title("Filtered beads")

# Combined plot 
adata_filtered = adata[filtered_barcodes].copy()
adata_filtered.obsm["X_spatial"] = adata_filtered.obsm["X_spatial"].copy()
# ----------------- first plot -----------------
title2 = 'Masked sample'

adata_filtered.obsm["X_spatial"] = adata_filtered.obsm["X_spatial"].copy()
x = adata_filtered.obsm["X_spatial"][:, 0]
y = adata_filtered.obsm["X_spatial"][:, 1]
colors = adata_filtered.obs["log10_numReads"]
min_val = colors.min()
max_val = colors.max()
print(f"Masked sample color range: {min_val:.2f} – {max_val:.2f}")

fig2, ax2 = plt.subplots(figsize=(8,6))
sc2 = ax2.scatter(x, y, c=colors, cmap="viridis", s=1)
ax2.invert_yaxis()
fig2.colorbar(sc2, ax=ax2, label="log10_numReads")
ax2.set_title(f"{title2}\n(n = {adata_filtered.n_obs})")

png2_path = output_base / f"{sample}_masked_sample.png"
fig2.savefig(png2_path, dpi=300)
plt.close(fig2)

# ----------------- second plot -----------------
title1 = 'Full sample'
adata = adata[bead_barcodes].copy()
adata.obsm["X_spatial"] = adata.obsm["X_spatial"].copy()

x = adata.obsm["X_spatial"][:, 0]
y = adata.obsm["X_spatial"][:, 1]
# clip full sample colors to masked range
full_colors = np.clip(adata.obs["log10_numReads"].values, min_val, max_val)
# normalize to same range
norm = mcolors.Normalize(vmin=min_val, vmax=max_val)

fig1, ax1 = plt.subplots(figsize=(8,6))
sc1 = ax1.scatter(x, y, c=full_colors, cmap="viridis", s=1)
ax1.invert_yaxis()
fig1.colorbar(sc1, ax=ax1, label="log10_numReads")
ax1.set_title(f"{title1}\n(n = {adata.n_obs})")

png1_path = output_base / f"{sample}_full_sample.png"
fig1.savefig(png1_path, dpi=300)
plt.close(fig1)  # close figure to free memory

# ----------------- stitch images together -----------------
img1 = Image.open(png1_path)
img2 = Image.open(png2_path)

# create new image side by side
total_width = img1.width + img2.width
max_height = max(img1.height, img2.height)

new_img = Image.new('RGB', (total_width, max_height), (255, 255, 255))
new_img.paste(img1, (0, 0))
new_img.paste(img2, (img1.width, 0))

new_img_path = output_base / f"{sample}_combined_samples.png"
new_img.save(new_img_path)

# # ----------------- first plot -----------------
# title1 = 'Full sample'
# adata.obsm["X_spatial"] = adata.obsm["X_spatial"].copy()

# # extract plot data 
# x = adata.obsm["X_spatial"][:, 0]
# y = adata.obsm["X_spatial"][:, 1]
# colors = adata.obs["log10_numReads"]

# fig1, ax1 = plt.subplots(figsize=(8,6))
# sc1 = ax1.scatter(
#     x, y, c=colors, cmap="viridis", s=1
# )
# ax1.invert_yaxis()
# fig1.colorbar(sc1, ax=ax1, label="log10_numReads")
# ax1.set_title(f"{title1}\n(n = {adata.n_obs})")

# # ----------------- second plot -----------------
# title2 = 'Masked sample'
# adata_filtered = adata[filtered_barcodes].copy()
# adata_filtered.obsm["X_spatial"] = adata_filtered.obsm["X_spatial"].copy()

# # extract plot data 
# x = adata_filtered.obsm["X_spatial"][:, 0]
# y = adata_filtered.obsm["X_spatial"][:, 1]
# colors = adata_filtered.obs["log10_numReads"]

# fig1, ax1 = plt.subplots(figsize=(8,6))
# sc1 = ax1.scatter(
#     x, y, c=colors, cmap="viridis", s=1
# )
# ax1.invert_yaxis()
# fig1.colorbar(sc1, ax=ax1, label="log10_numReads")
# ax1.set_title(f"{title2}\n(n = {adata_filtered.n_obs})")


# # 4. Save figure
# fig_path = output_base / f"{sample}_{adata_filtered.n_obs}_masked_spatial.png"
# plt.savefig(fig_path, dpi=300)

# 5. Save adata filtered
# adata_fil_path = output_base / f'{sample}_{adata_filtered.n_obs}_masked.h5ad'
# adata_filtered.write_h5ad(adata_fil_path)

# 5. Save original adata with a new metadata column 
adata_meta_path = output_base / f'{sample}_full_in_mask_{adata_filtered.n_obs}.h5ad'
adata.write_h5ad(adata_meta_path)