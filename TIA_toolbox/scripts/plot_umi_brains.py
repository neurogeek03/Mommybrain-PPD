#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def main():
    parser = argparse.ArgumentParser(
        description="Make UMI-colored spatial plots + padded frames"
    )

    parser.add_argument("-s", "--sample", required=True,
                        help="Sample name (e.g. B01)")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to input .h5ad file")
    parser.add_argument("-o", "--output", required=True,
                        help="Directory for downstream pipeline outputs")
    parser.add_argument("-p", "--plots", required=True,
                        help="Directory where plots will be saved")
    parser.add_argument("-sc", "--scalar", required=True,
                        help="scalar for umi plot")

    args = parser.parse_args()

    sample = args.sample
    adata_path = Path(args.input)
    output_dir = Path(args.output)
    plots_dir = Path(args.plots)
    scalar = args.scalar

    timestamp = datetime.now().strftime("%M%S%H_%Y%m%d")
    print(f"------ Plotting sample {sample} at {timestamp} ------")
    print(f"Reading: {adata_path}")

    # =================== PARAMS ===================
    dpi = 600
    pad = 500
    pad_color = 255

    # =================== Read adata ===================
    adata = sc.read_h5ad(adata_path)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    coords = pd.DataFrame(adata.obsm["X_spatial"], columns=["x", "y"])
    color_variable = "nCount_RNA"

    # Calculate scale clipping maximum 
    umi_per_cell = np.array(adata.X.sum(axis=1)).ravel()
    median_umi = np.median(umi_per_cell)
    print("Median UMI per cell:", median_umi)

    max_umi = median_umi * scalar # for some samples I used 1.8 or 2.8, based on trial & error (manually inspect masks)
    # =================== UMI PLOT ===================
    coords[color_variable] = adata.obs[color_variable].values
    clipped_colors = np.clip(coords[color_variable], 0, max_umi)
    norm = mcolors.Normalize(vmin=0, vmax=max_umi)

    fig = plt.figure(figsize=(12, 12))
    ax = plt.gca()

    ax.scatter(
        coords["x"], coords["y"],
        c=clipped_colors,
        s=0.1,
        cmap="Reds",
        norm=norm
    )

    ax.set_title("")
    ax.set_aspect("equal")
    ax.axis("off")
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    spatial_png = plots_dir / f"{sample}_spatial_{timestamp}.png"
    plt.savefig(spatial_png, dpi=dpi, bbox_inches="tight", pad_inches=0.0)
    plt.close()

    print(f"Saved spatial plot at: {spatial_png}")

    # =================== Process w numpy ===================
    img = np.array(Image.open(spatial_png))
    print(f"Original image shape: {img.shape}")

    gray = img if img.ndim == 2 else img.mean(axis=2)
    threshold = 250

    rows = np.where((gray < threshold).any(axis=1))[0]
    cols = np.where((gray < threshold).any(axis=0))[0]

    cropped = img[rows.min():rows.max()+1, cols.min():cols.max()+1]
    print(f"Cropped shape: {cropped.shape}")

    cropped_path = plots_dir / f"{sample}_cropped_{timestamp}.png"
    Image.fromarray(cropped).save(cropped_path)

    # =================== Add padding ===================
    padded = np.pad(
        cropped,
        ((pad, pad), (pad, pad), (0, 0)),
        mode="constant",
        constant_values=pad_color
    )

    padded_path = output_dir / f"{sample}_padded.png"
    Image.fromarray(padded).save(padded_path)

    print(f"Padded image saved at: {padded_path}")
    print("Done.")

if __name__ == "__main__":
    main()
