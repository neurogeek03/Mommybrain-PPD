from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import requests

from tiatoolbox.tools.tissuemask import MorphologicalMasker
from tiatoolbox.utils import imwrite
from tiatoolbox.wsicore.wsireader import WSIReader

stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
print(f'------ Script started at {stamp} ------')

project_path = Path.cwd().parents[0]
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True,parents=True)
masks = output_base / "masks"
masks.mkdir(exist_ok=True,parents=True)
print(f'Current working dir  {project_path}')

# PARAMS
res = 0.7
img_path = project_path / 'data' / 'log10_reads.png'

# importing own data in PNG format
# see allowed formats at: https://github.com/TissueImageAnalytics/tiatoolbox/blob/master/tiatoolbox/wsicore/wsireader.py
wsi_object = WSIReader.open(input_img=img_path)
# access img metadata
image_array = wsi_object.img
print(f"Image shape: {image_array.shape}") # shown in (height, weight, RGB)

wsi_info = wsi_object.info.as_dict()
# Print one item per line
print(*list(wsi_info.items()), sep="\n")  # noqa: T201

wsi_thumb = wsi_object.slide_thumbnail(resolution=0.95, units="baseline")
plt.imshow(wsi_thumb)
plt.axis("off")
thumbnail_path = output_base / "wsi_thumbnail_plot.png"
plt.savefig(thumbnail_path, bbox_inches="tight", pad_inches=0)

# create tissue mask 
mask = wsi_object.tissue_mask(resolution=res, units="baseline")
mask_thumb = mask.slide_thumbnail(resolution=res, units="baseline")
# save side-by-side
def save_side_by_side(image_1: np.ndarray, image_2: np.ndarray, filepath: str) -> None:
    """Combines two images side-by-side and saves the resulting figure."""
    
    # 1. Create a new figure and a 1x2 grid of subplots
    plt.figure(figsize=(10, 5)) # Set a reasonable figure size for two images
    
    # 2. Plot Image 1 (Left side)
    plt.subplot(1, 2, 1)
    plt.imshow(image_1)
    plt.axis("off")
    
    # 3. Plot Image 2 (Right side)
    plt.subplot(1, 2, 2)
    plt.imshow(image_2)
    plt.axis("off")
    
    # 4. Save the combined figure to the specified filepath
    # 'bbox_inches="tight"' crops the image to the content, removing excess white space.
    plt.savefig(filepath, bbox_inches="tight", pad_inches=0)
    
    # 5. Close the figure to free up memory (important when processing many images)
    plt.close()
    
    print(f"Image successfully saved to: {filepath}")

# --- Example Usage ---

# Assuming wsi_thumb and mask_thumb are your NumPy arrays
mask_root = f"{res}_{stamp}"
output_path = output_base / f"{mask_root}_side_by_side_comparison.png"
save_side_by_side(wsi_thumb, mask_thumb, filepath=output_path)

# saving array 
mask_array = mask.img
print(mask_array)
print(mask_array.shape)
mask_path = masks / mask_root
np.save(mask_path, mask_array)

print(f'mask array saved in {mask_path}')