import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import plotly.express as px

project_path = Path.cwd()
out_dir = project_path / 'out'
figs = out_dir / 'figures'
figs.mkdir(exist_ok=True, parents=True)

bar_path = figs / 'bar_sccoda_proportions.png'
pie_path = figs / 'pie_sccoda_proportions.png'
html_path = figs / 'cell_type_composition.html'
credible_thresh = 0.95
threshold = 0.007 # threshold to show labels on pie chart

# ============ GET COLORS ============
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])
# color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)
# Create mapping from label number (prefix before _) to hex color
def get_num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except:
        return -1
color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)
color_df = color_df.sort_values("num_prefix")
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))
label_to_hex['Other'] = '#CCCCCC' 


# Example data: replace with your full intercepts and effects CSVs
intercepts_path = out_dir / "sccoda_intercepts.csv"
intercepts = pd.read_csv(intercepts_path)

effects_path = out_dir / "sccoda_effects.csv"
effects = pd.read_csv(effects_path)

# Convert log-abundances to proportions (softmax)
def log_to_prop(log_values):
    exp_vals = np.exp(log_values)
    return exp_vals / exp_vals.sum()

# Baseline proportions
baseline_props = log_to_prop(intercepts["Final Parameter"].values)
# After treatment
after_treatment_props = log_to_prop(intercepts["Final Parameter"].values + effects["Final Parameter"].values)

df = pd.DataFrame({
    "Cell Type": intercepts["Cell Type"],
    "Baseline": baseline_props,
    "Treatment": after_treatment_props,
    "Inclusion": effects["Inclusion probability"]
})

# Sort by baseline proportion (high to low)
df = df.sort_values("Baseline", ascending=False).reset_index(drop=True)

# Plot
x = np.arange(len(intercepts))
width = 0.35

fig, ax = plt.subplots(figsize=(25,10))
ax.bar(x - width/2, df["Baseline"], width, label='Baseline', color='lightblue')
ax.bar(x + width/2, df["Treatment"], width, label='Treatment', color='gray')

# Add asterisks for credible changes
for i, p in enumerate(df["Inclusion"]):
    if p >= credible_thresh:
        ax.text(x[i] + width/2, df["Treatment"][i] + 0.01, "*", ha='center', va='bottom', color='red', fontsize=14)

ax.set_xticks(x)  # x = positions of the tick centers
ax.set_xticklabels(df["Cell Type"], rotation=90, ha="center")  # center labels
ax.set_ylabel("Predicted Proportion")
ax.set_title("Cell Type Proportions Before and After Treatment")
ax.margins(x=0)
ax.legend()
plt.tight_layout()
plt.savefig(bar_path)

# =================== PIE CHART ===================
# grouping into broader categories 
# df_grouped = df.copy()
# df_grouped["Group"] = df_grouped["Cell Type"].apply(lambda x: x.split()[-1])
# baseline_group = df_grouped.groupby("Group")["Baseline"].sum()
# treatment_group = df_grouped.groupby("Group")["Treatment"].sum()
# baseline_group = baseline_group / baseline_group.sum()
# treatment_group = treatment_group / treatment_group.sum()
# baseline_labels = [f"{grp}\n{p*100:.1f}%" for grp, p in zip(baseline_group.index, baseline_group)]
# treatment_labels = [f"{grp}\n{p*100:.1f}%" for grp, p in zip(treatment_group.index, treatment_group)]

# Normalize to proportions (in case they don’t already sum to 1 due to model averaging)
baseline_props = df["Baseline"] / df["Baseline"].sum()
treatment_props = df["Treatment"] / df["Treatment"].sum()
# Convert to percent labels 
baseline_labels = [f"{ct}\n{p*100:.1f}%" if p > threshold else ""
                   for ct, p in zip(df["Cell Type"], baseline_props)]

treatment_labels = [f"{ct}\n{p*100:.1f}%" if p > threshold else ""
                    for ct, p in zip(df["Cell Type"], treatment_props)]

fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Baseline pie
axes[0].pie(baseline_props, labels=baseline_labels, 
            textprops={'fontsize': 12}, 
            wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
axes[0].set_title("OIL - grouped cell type proportions")

# Treatment pie
axes[1].pie(treatment_props, labels=treatment_labels, 
            textprops={'fontsize': 12}, 
            wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
axes[1].set_title("CORT - grouped cell type proportions")

plt.tight_layout()
plt.show()
plt.savefig(pie_path, dpi=300)

# =================== PLOTLY ===================
# Build long-format DataFrame for Plotly
baseline_df = pd.DataFrame({
    'Cell Type': df['Cell Type'],
    'Proportion': baseline_props,
    'Condition': 'Baseline'
})

treatment_df = pd.DataFrame({
    'Cell Type': df['Cell Type'],
    'Proportion': treatment_props,
    'Condition': 'Treatment'
})

plot_df = pd.concat([baseline_df, treatment_df], ignore_index=True)
# Optional: group very small slices into "Other"
def group_small(df, threshold=0.01):
    mask = df['Proportion'] < threshold
    if mask.any():
        small_sum = df.loc[mask, 'Proportion'].sum()
        df = df.loc[~mask].copy()
        df = pd.concat([df, pd.DataFrame([{
            'Cell Type':'Other', 
            'Proportion': small_sum,
            'Condition': df['Condition'].iloc[0]
        }])], ignore_index=True)
    return df

baseline_grouped = group_small(baseline_df, threshold)
treatment_grouped = group_small(treatment_df, threshold)
plot_df_grouped = pd.concat([baseline_grouped, treatment_grouped], ignore_index=True)
# Plot with Plotly Express
fig = px.pie(
    plot_df_grouped,
    values='Proportion',
    names='Cell Type',
    facet_col='Condition',    # Baseline vs Treatment side by side
    color='Cell Type',         # consistent coloring
    color_discrete_map=label_to_hex,
    title='Cell Type Composition: Baseline vs Treatment'
)

# Show percentages on the wedges
fig.update_traces(textinfo='percent')
# Legend always visible
fig.write_html(html_path)
