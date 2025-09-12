import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the CSV
df = pd.read_csv("/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/post_RCTD_evaluation/all_samples_coronal_run/d5_umi100_RCTD_first_type_mouse_sample_counts.csv")

# Sum counts across all samples for each RCTD type
total_counts = df.groupby("RCTD_first_type_mouse")["count"].sum().reset_index()

# Select the top 35 most prevalent types
top35_types = total_counts.nlargest(35, "count")["RCTD_first_type_mouse"]

# Filter original df to keep only top 35
df_top35 = df[df["RCTD_first_type_mouse"].isin(top35_types)]

# Order RCTD types by total counts
order = total_counts[total_counts["RCTD_first_type_mouse"].isin(top35_types)] \
        .sort_values("count", ascending=False)["RCTD_first_type_mouse"]

df_top35["RCTD_first_type_mouse"] = pd.Categorical(df_top35["RCTD_first_type_mouse"], categories=order, ordered=True)

# Plot using seaborn
plt.figure(figsize=(12, 8))
sns.barplot(
    data=df_top35,
    x="RCTD_first_type_mouse",
    y="count",
    hue="sample"
)
plt.xticks(rotation=90)
plt.xlabel("RCTD First Type (Mouse)")
plt.ylabel("Cell Counts")
plt.title("Top 35 Most Prevalent RCTD First Types Across Samples")
plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/post_RCTD_evaluation/all_samples_coronal_run/ordered_counts_coronal.png", dpi=300, bbox_inches='tight')
