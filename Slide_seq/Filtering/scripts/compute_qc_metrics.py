import os 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

# Set paths
singlet_score_cutoff = 330

input_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/OUT/FINAL_RCTD_newgenelist/merged_metadata_csvs"
metadata_experimental = "/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR/slide_seq_metadata.csv"
output_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Filtering/OUT/FINAL_rctd_run"
os.makedirs(output_dir, exist_ok=True)

metadata_experimental = pd.read_csv(metadata_experimental)

early_samples = [
    "B08",
    "B19",
    "B33",
    "B21",
    "B01",
    "B18",
    "B37",
    "B42"
]

# Make a list of all files 
dfs = [] 
for file in os.listdir(input_dir):
    df = pd.read_csv(os.path.join(input_dir, file))
    sample = file.split("_")[0]
    df["sample"] = sample
    df["singlet"] = np.where(df["spot_class"] == "singlet", "singlet", "other")
    dfs.append(df)

metadata = pd.concat(dfs)
metadata = metadata.merge(metadata_experimental, on="sample", how="left")
metadata["group"] = metadata["pregnancy"] + "_" + metadata["day"] + "_" + metadata["treatment"]
print(metadata.head())
print(metadata["group"].unique())

sample_order = []

# Loop through your desired hue order and collect samples for each group
sample_order = []
for g in ['Nulliparous_PD8_NONE','Preg_GD20_NONE','Post_PD8_NONE', 'Post_PD8_OIL','Post_PD8_CORT', 'Post_PD23_OIL','Post_PD23_CORT']:
    samples_in_group = metadata.loc[metadata['group'] == g, 'sample'].unique()
    sample_order.extend(samples_in_group)

plt.figure(figsize=(14,8))

# Initial figure
sns.violinplot(
    data=metadata, 
    x="sample",
    y="singlet_score",
    hue="group",
    hue_order=['Nulliparous_PD8_NONE','Preg_GD20_NONE','Post_PD8_NONE', 'Post_PD8_OIL','Post_PD8_CORT', 'Post_PD23_OIL','Post_PD23_CORT'],
    palette="Set2",
    order=sample_order,
    # split = True,
    # inner="quart",
    cut=0
)
count = metadata.loc[(metadata['spot_class'] == 'singlet') & (metadata['singlet_score'] > 330)].shape[0]

# Draw a horizontal line at y=330
plt.axhline(y=330, color='red', linestyle='--', linewidth=1.5)

plt.legend(title="Pregnancy & Treatment groups", bbox_to_anchor=(1.01, 1), loc='upper left')
plt.xticks(rotation=45)

ax = plt.gca()
red_samples = early_samples  # samples you want in red
for tick in ax.get_xticklabels():
    if tick.get_text() in red_samples:
        tick.set_color("red")

plt.tight_layout()
plt.title(f"Singlet score per sample: cutoff = {singlet_score_cutoff}. Total singlets = {count}")
plt.savefig(os.path.join(output_dir, "groups_samples_violin.png") ,dpi=300, bbox_inches='tight')

