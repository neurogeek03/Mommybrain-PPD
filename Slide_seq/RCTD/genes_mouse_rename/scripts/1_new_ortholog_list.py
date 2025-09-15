import os
import pandas as pd
import numpy as np

rat_to_mouse = pd.read_csv('/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/rat_to_mouse.csv')
output_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/out"

pd.set_option("display.max_columns", None)
print(rat_to_mouse.head())

rat_to_mouse = rat_to_mouse.dropna(
    subset=rat_to_mouse.columns.difference(["Gene stable ID", "Gene name"]),
    how="all"
)
print(rat_to_mouse.head())
# print(len(rat_to_mouse))

# print(len(rat_to_mouse["Gene stable ID"].unique()))
# print(len(rat_to_mouse["Mouse gene stable ID"].unique()))

df_filtered = (
    rat_to_mouse.sort_values(
        by="%id. target Mouse gene identical to query gene", ascending=False
    )
    .drop_duplicates(subset="Gene stable ID", keep="first")
)


# print(len(df_filtered["Gene stable ID"].unique()))
# print(len(df_filtered["Mouse gene stable ID"].unique()))
# print(len(df_filtered))
# print(df_filtered.head(30))

df_filtered_backwards = (
    df_filtered.sort_values(
        by="%id. query gene identical to target Mouse gene", ascending=False
    )
    .drop_duplicates(subset="Mouse gene stable ID", keep="first")
)

print(len(df_filtered_backwards["Gene stable ID"].unique()))
print(len(df_filtered_backwards["Mouse gene stable ID"].unique()))
print(len(df_filtered_backwards))

# # Count unique mouse genes per rat gene
# counts = rat_to_mouse.groupby("Gene stable ID")["Mouse gene stable ID"].nunique()

# # Group by the number of unique mouse genes and count how many rat genes fall in each group
# summary = counts.value_counts().sort_index()

# print(summary)

df_filtered_backwards.to_csv(os.path.join(output_dir,"1-1_rat_to_mouse_filtered.csv"), index=False)