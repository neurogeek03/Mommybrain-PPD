import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import plotly.express as px
from adjustText import adjust_text

# ============ params ============
comparison_folder = 'rostral_cort_vs_non_cort'
comparison = 'condition[T.NONE]'

credible_thresh = 0.95
threshold = 0.007 # threshold to show labels on pie chart

# ============ paths ============
project_path = Path.cwd()
out_dir = project_path / 'out'/'slide_seq' /comparison_folder
figs = out_dir / 'figures' 
figs.mkdir(exist_ok=True, parents=True)

effects_path = out_dir / "sccoda_effects.csv"
intercepts_path = out_dir / "sccoda_intercepts.csv"

bar_path = figs / 'bar_sccoda_proportions.png'
pie_path = figs / 'pie_sccoda_proportions.png'
html_path = figs / 'cell_type_composition.html'
volcano_path = figs / f'error_bar_volcano_sccoda_{comparison}.png'

# ============ get ABC atlas colors (WMB) ============
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


# read in data 
intercepts = pd.read_csv(intercepts_path)
effects = pd.read_csv(effects_path)

# # subset effects result
# effects = effects[effects["Covariate"] == "coronal_section[T.rostral]"].copy()
# effects = effects.reset_index(drop=True)

print(effects.head())
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

# # =================== PIE CHART ===================
# # grouping into broader categories 
# # df_grouped = df.copy()
# # df_grouped["Group"] = df_grouped["Cell Type"].apply(lambda x: x.split()[-1])
# # baseline_group = df_grouped.groupby("Group")["Baseline"].sum()
# # treatment_group = df_grouped.groupby("Group")["Treatment"].sum()
# # baseline_group = baseline_group / baseline_group.sum()
# # treatment_group = treatment_group / treatment_group.sum()
# # baseline_labels = [f"{grp}\n{p*100:.1f}%" for grp, p in zip(baseline_group.index, baseline_group)]
# # treatment_labels = [f"{grp}\n{p*100:.1f}%" for grp, p in zip(treatment_group.index, treatment_group)]

# # Normalize to proportions (in case they don’t already sum to 1 due to model averaging)
# baseline_props = df["Baseline"] / df["Baseline"].sum()
# treatment_props = df["Treatment"] / df["Treatment"].sum()
# # Convert to percent labels 
# baseline_labels = [f"{ct}\n{p*100:.1f}%" if p > threshold else ""
#                    for ct, p in zip(df["Cell Type"], baseline_props)]

# treatment_labels = [f"{ct}\n{p*100:.1f}%" if p > threshold else ""
#                     for ct, p in zip(df["Cell Type"], treatment_props)]

# fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# # Baseline pie
# axes[0].pie(baseline_props, labels=baseline_labels, 
#             textprops={'fontsize': 12}, 
#             wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
# axes[0].set_title("OIL - grouped cell type proportions")

# # Treatment pie
# axes[1].pie(treatment_props, labels=treatment_labels, 
#             textprops={'fontsize': 12}, 
#             wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
# axes[1].set_title("CORT - grouped cell type proportions")

# plt.tight_layout()
# plt.show()
# plt.savefig(pie_path, dpi=300)

# # =================== PLOTLY ===================
# # Build long-format DataFrame for Plotly
# baseline_df = pd.DataFrame({
#     'Cell Type': df['Cell Type'],
#     'Proportion': baseline_props,
#     'Condition': 'Baseline'
# })

# treatment_df = pd.DataFrame({
#     'Cell Type': df['Cell Type'],
#     'Proportion': treatment_props,
#     'Condition': 'Treatment'
# })

# plot_df = pd.concat([baseline_df, treatment_df], ignore_index=True)
# # Optional: group very small slices into "Other"
# def group_small(df, threshold=0.01):
#     mask = df['Proportion'] < threshold
#     if mask.any():
#         small_sum = df.loc[mask, 'Proportion'].sum()
#         df = df.loc[~mask].copy()
#         df = pd.concat([df, pd.DataFrame([{
#             'Cell Type':'Other', 
#             'Proportion': small_sum,
#             'Condition': df['Condition'].iloc[0]
#         }])], ignore_index=True)
#     return df

# baseline_grouped = group_small(baseline_df, threshold)
# treatment_grouped = group_small(treatment_df, threshold)
# plot_df_grouped = pd.concat([baseline_grouped, treatment_grouped], ignore_index=True)
# # Plot with Plotly Express
# fig = px.pie(
#     plot_df_grouped,
#     values='Proportion',
#     names='Cell Type',
#     facet_col='Condition',    # Baseline vs Treatment side by side
#     color='Cell Type',         # consistent coloring
#     color_discrete_map=label_to_hex,
#     title='Cell Type Composition: Baseline vs Treatment'
# )

# # Show percentages on the wedges
# fig.update_traces(textinfo='percent')
# # Legend always visible
# fig.write_html(html_path)


# =================== VOLCANO ===================
df = effects

ip_threshold = 0.8
effect_threshold = 0  # zero: no change direction

# Highlight significant effects
df['significant'] = (df['Inclusion probability'] > ip_threshold)

xerr = [
    df['Final Parameter'] - df['HDI 3%'],
    df['HDI 97%'] - df['Final Parameter']
]

# Colors for points
colors = df['significant'].map({True: 'red', False: 'gray'}).tolist()

# plt.figure(figsize=(10,12))

# # Scatter for points with per-point colors
# plt.scatter(
#     df['Final Parameter'],
#     df['Inclusion probability'],
#     s=120,
#     alpha=0.8,
#     c=colors
# )

# # Plot horizontal error bars for HDI
# for xi, yi, xe in zip(df['Final Parameter'], df['Inclusion probability'], zip(*xerr)):
#     plt.errorbar(xi, yi, xerr=[[xe[0]], [xe[1]]], fmt='none', ecolor='gray', elinewidth=1, capsize=3)

# # Threshold lines
# plt.axhline(0.8, linestyle='--', color='black', alpha=0.6)
# plt.axvline(0, linestyle='--', color='black', alpha=0.6)

fig = plt.figure(figsize=(10,6))
gs = GridSpec(2, 1, height_ratios=[3, 1])  # top:bottom ratio

ax1 = fig.add_subplot(gs[0])  # top axis
ax2 = fig.add_subplot(gs[1], sharex=ax1)  # bottom axis

# Plot same points in both axes
ax1.scatter(df['Final Parameter'], df['Inclusion probability'], c=colors)
ax2.scatter(df['Final Parameter'], df['Inclusion probability'], c=colors)

# Set y-limits
ax1.set_ylim(0.8, 1.05)   # zoomed-in top
ax2.set_ylim(0, 0.8)      # zoomed-out bottom

# Hide spines between the axes for “break”
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(bottom=False)
ax2.xaxis.tick_bottom()

# Add diagonal break marks if you want
d = .010
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d,+d),(-d,+d), **kwargs)
ax1.plot((1-d,1+d),(-d,+d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d,+d),(1-d,1+d), **kwargs)
ax2.plot((1-d,1+d),(1-d,1+d), **kwargs)

# Plot horizontal error bars for HDI on ax1
for xi, yi, xe in zip(df['Final Parameter'], df['Inclusion probability'], zip(*xerr)):
    ax1.errorbar(xi, yi, xerr=[[xe[0]], [xe[1]]], fmt='none', ecolor='gray', elinewidth=1, capsize=3)

# Threshold lines on ax1
ax1.axhline(0.8, linestyle='--', color='black', alpha=0.6)
ax1.axvline(0, linestyle='--', color='black', alpha=0.6)

# for _, row in df[df['significant']].iterrows():
#     plt.text(row['Final Parameter'], row['Inclusion probability'] + 0.02, row['Cell Type'],
#              ha='center', fontsize=5)

texts = []
for _, row in df[df['significant']].iterrows():
    texts.append(
        ax1.text(
            row['Final Parameter'], 
            row['Inclusion probability'] + 0.01, 
            row['Cell Type'],
            ha='center', 
            fontsize=6
        )
    )

adjust_text(texts, ax=ax1, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax2.set_xlabel("Effect size (posterior mean, log-ratio)")
ax1.set_ylabel("Inclusion Probability")
ax1.set_title("Effect size vs. Inclusion Probability, error bars: Highest Density Interval (HDI)")
# plt.ylim(0, 1.05)
plt.tight_layout()
plt.savefig(volcano_path)


# plt.figure(figsize=(10, 12))
# scatter = plt.scatter(
#     df['Final Parameter'],
#     df['Inclusion probability'],
#     s=120,
#     alpha=0.8,
#     c=df['significant'].map({True: 'red', False: 'gray'})
# )

# # Styling
# plt.axhline(ip_threshold, linestyle='--', color='black', alpha=0.6)
# plt.axvline(effect_threshold, linestyle='--', color='black', alpha=0.6)

# plt.xlabel("Effect size (posterior mean, log-ratio)")
# plt.ylabel("Inclusion Probability")
# plt.title("Volcano-style plot using Inclusion Probability")

# # Optional: label significant points
# for _, row in df[df['significant']].iterrows():
#     plt.text(row['Final Parameter'], row['Inclusion probability'] + 0.02, row['Cell Type'],
#              ha='center', fontsize=8)

# plt.ylim(0, 1.05)
# plt.tight_layout()
# plt.savefig(volcano_path)