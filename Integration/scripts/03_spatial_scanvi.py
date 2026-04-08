"""
03_spatial_scanvi.py
--------------------
Interactive spatial viewer for scANVI-predicted cell types on Slide-seq sections.

Loads integrated.h5ad (output of 01_scanvi.py), filters to Slide-seq cells,
and generates a two-panel interactive HTML with:
  - Side-by-side spatial sections (selectable via dropdown)
  - Cells colored by scANVI predicted label
  - Hover-to-highlight by cell type across both panels

Optionally filters low-confidence predictions (scanvi_max_prob < CONF_THRESHOLD).

Usage:
  python 03_spatial_scanvi.py [config_path]
"""

import os
import re
import json
import argparse
import pandas as pd
import scanpy as sc
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import yaml
from pathlib import Path

# =============================================================================
# CONFIG
# =============================================================================

parser = argparse.ArgumentParser()
parser.add_argument(
    "config",
    nargs="?",
    default=Path(__file__).parent.parent / "config" / "01_scanvi.conf",
)
args = parser.parse_args()

with open(args.config) as f:
    cfg = yaml.safe_load(f)

OUT_DIR     = Path(cfg["out_dir"])
CONFIG_NAME = cfg["config_name"]
REF_METHOD  = cfg.get("reference_method", "slide_tags")
COLOR_CSV   = "/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv"
META_CSV    = Path(__file__).parent.parent / "meta" / "slideseq_metadata.csv"

CONFIG_DIR  = OUT_DIR / CONFIG_NAME
INTEGRATED_PATH = CONFIG_DIR / "integrated.h5ad"
VIZ_DIR     = CONFIG_DIR / "viz_scanvi"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_HTML = VIZ_DIR / "08_spatial_scanvi_prediction.html"

# Confidence filter — set to None to disable
CONF_THRESHOLD = 0.8

# =============================================================================
# COLORS
# =============================================================================

color_df = pd.read_csv(COLOR_CSV, usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)


def _num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except Exception:
        return -1


color_df["_num"] = color_df["name"].apply(_num_prefix)
color_df = color_df.sort_values("_num")
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))


def get_color(ct):
    hex_val = label_to_hex.get(ct, "888888")
    return f"#{hex_val}" if not hex_val.startswith("#") else hex_val


# =============================================================================
# LOAD + FILTER
# =============================================================================

print(f"[{CONFIG_NAME}] Loading integrated object...")
adata = sc.read_h5ad(INTEGRATED_PATH)
print(f"  Shape: {adata.shape}")
print(f"  obsm keys: {list(adata.obsm.keys())}")

# Filter to Slide-seq cells only
adata = adata[adata.obs["method"] != REF_METHOD].copy()
print(f"  Slide-seq cells: {adata.n_obs}")

# Optional confidence filter
if CONF_THRESHOLD is not None:
    n_before = adata.n_obs
    adata = adata[adata.obs["scanvi_max_prob"] >= CONF_THRESHOLD].copy()
    print(f"  After conf≥{CONF_THRESHOLD} filter: {adata.n_obs} "
          f"({n_before - adata.n_obs} removed, "
          f"{100*(n_before - adata.n_obs)/n_before:.1f}%)")

assert "X_spatial" in adata.obsm, "Missing obsm['X_spatial'] — check merged_raw.h5ad"
assert "scanvi_prediction" in adata.obs.columns, "Missing obs['scanvi_prediction']"
assert "sample_id" in adata.obs.columns, "Missing obs['sample_id']"

# =============================================================================
# BUILD MAIN DATAFRAME
# =============================================================================

coords = adata.obsm["X_spatial"]
df = pd.DataFrame(coords, columns=["x", "y"], index=adata.obs_names)
df["prediction"] = adata.obs["scanvi_prediction"].astype(str).values
df["sample_id"]  = adata.obs["sample_id"].values
df["max_prob"]   = adata.obs["scanvi_max_prob"].values

# Join metadata
meta_df = pd.read_csv(META_CSV)
meta_df.columns = meta_df.columns.str.strip()
meta_df["Sample"] = meta_df["Sample"].str.strip()
meta_indexed = meta_df.set_index("Sample")
for col in meta_indexed.columns:
    df[col] = df["sample_id"].map(meta_indexed[col])

# =============================================================================
# SAMPLE LABELS (for dropdown)
# =============================================================================

available_samples = sorted(df["sample_id"].dropna().unique().tolist())
sample_label = {}
for sid in available_samples:
    row = meta_df[meta_df["Sample"] == sid]
    if row.empty:
        sample_label[sid] = sid
    else:
        r = row.iloc[0]
        sample_label[sid] = (
            f"{sid} — {r.get('Day', '')} | "
            f"{r.get('Treatment', '')} | "
            f"{r.get('Position', '')}"
        )
print(f"  Samples: {available_samples}")

# =============================================================================
# SERIALIZE SPATIAL DATA AS JSON (loaded by JS on dropdown change)
# =============================================================================

all_predictions = sorted(df["prediction"].dropna().unique().tolist())

sample_data = {}
for sid in available_samples:
    sdf = df[df["sample_id"] == sid]
    if sdf.empty:
        continue
    sample_data[sid] = {
        "x":  [round(v, 1) for v in sdf["x"].tolist()],
        "y":  [round(v, 1) for v in sdf["y"].tolist()],
        "ct": sdf["prediction"].tolist(),
    }

color_map = {ct: get_color(ct) for ct in all_predictions}

# =============================================================================
# BUILD FIGURE (empty spatial panels; JS populates on load)
# =============================================================================

fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=["Left Panel", "Right Panel"],
    column_widths=[0.5, 0.5],
    horizontal_spacing=0.03,
)

# Invisible dummy traces to force legend entries (one per predicted type)
for ct in all_predictions:
    fig.add_trace(
        go.Scattergl(
            x=[None], y=[None],
            mode="markers",
            marker=dict(size=6, color=get_color(ct)),
            name=ct,
            legendgroup=ct,
            showlegend=True,
            meta=ct,
            hoverinfo="skip",
        ),
        row=1, col=1,
    )

fig.update_layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    font=dict(color="#222222", size=11),
    showlegend=True,
    legend=dict(
        x=1.01, y=1,
        bgcolor="rgba(255,255,255,0)",
        bordercolor="#dddddd",
        borderwidth=1,
        font=dict(size=9, color="#333"),
        tracegroupgap=2,
        itemsizing="constant",
    ),
    margin=dict(l=10, r=180, t=55, b=10),
    height=650,
)

for col in [1, 2]:
    fig.update_xaxes(visible=False, showgrid=False, zeroline=False, row=1, col=col)
    fig.update_yaxes(
        visible=False, showgrid=False, zeroline=False,
        autorange="reversed",
        scaleanchor=("x" if col == 1 else "x2"),
        scaleratio=1,
        row=1, col=col,
    )

for ann in fig.layout.annotations:
    ann.font.color = "#333333"
    ann.font.size  = 11

# =============================================================================
# DROPDOWN HTML
# =============================================================================

default_left  = available_samples[0] if available_samples else ""
default_right = available_samples[1] if len(available_samples) > 1 else default_left

left_options = "\n".join(
    f'<option value="{s}"{" selected" if s == default_left else ""}>'
    f'{sample_label.get(s, s)}</option>'
    for s in available_samples
)
right_options = "\n".join(
    f'<option value="{s}"{" selected" if s == default_right else ""}>'
    f'{sample_label.get(s, s)}</option>'
    for s in available_samples
)

conf_note = f"Showing cells with prediction confidence ≥ {CONF_THRESHOLD}" if CONF_THRESHOLD else "All Slide-seq cells shown"

# =============================================================================
# JAVASCRIPT
# =============================================================================

js = f"""
const sampleData = {json.dumps(sample_data)};
const colorMap   = {json.dumps(color_map)};
const sampleMeta = {json.dumps(sample_label)};

let currentLeft  = '{default_left}';
let currentRight = '{default_right}';

const gd = document.getElementById('main_fig');

// After initial render, load default samples
gd.on('plotly_afterplot', function initView() {{
    gd.removeListener('plotly_afterplot', initView);
    updateView();
    attachHoverEvents();
}});

function buildSpatialTraces(sampleId, col) {{
    if (!sampleId || !sampleData[sampleId]) return [];
    const d = sampleData[sampleId];
    const n = d.x.length;
    const xax = col === 1 ? 'x'  : 'x2';
    const yax = col === 1 ? 'y'  : 'y2';
    const groups = {{}};
    for (let i = 0; i < n; i++) {{
        const ct = d.ct[i];
        if (!groups[ct]) groups[ct] = {{x: [], y: []}};
        groups[ct].x.push(d.x[i]);
        groups[ct].y.push(d.y[i]);
    }}
    const traces = [];
    for (const [ct, pts] of Object.entries(groups)) {{
        traces.push({{
            type: 'scattergl',
            x: pts.x, y: pts.y,
            mode: 'markers',
            marker: {{size: 3, color: colorMap[ct] || '#888', opacity: 0.85}},
            name: ct,
            legendgroup: ct,
            showlegend: false,
            meta: ct,
            xaxis: xax, yaxis: yax,
            hovertemplate: sampleId + ' | ' + ct + '<extra></extra>',
        }});
    }}
    return traces;
}}

function updateView() {{
    // Legend dummy traces (indices 0..N-1) + spatial traces
    const legendTraces = gd.data.slice(0, {len(all_predictions)});
    const leftTraces   = buildSpatialTraces(currentLeft,  1);
    const rightTraces  = buildSpatialTraces(currentRight, 2);
    const allTraces    = [...legendTraces, ...leftTraces, ...rightTraces];

    Plotly.react(gd, allTraces, gd.layout);

    const leftLabel  = sampleMeta[currentLeft]  || currentLeft;
    const rightLabel = sampleMeta[currentRight] || currentRight;
    Plotly.relayout(gd, {{
        'annotations[0].text': leftLabel,
        'annotations[1].text': rightLabel,
    }});
}}

document.getElementById('selectLeft').addEventListener('change', function() {{
    currentLeft = this.value;
    updateView();
}});
document.getElementById('selectRight').addEventListener('change', function() {{
    currentRight = this.value;
    updateView();
}});

function attachHoverEvents() {{
    gd.on('plotly_hover', function(eventData) {{
        const hovered = eventData.points[0].data.meta;
        if (!hovered) return;
        const match = [], other = [];
        gd.data.forEach(function(trace, i) {{
            if (trace.meta === hovered) match.push(i);
            else other.push(i);
        }});
        if (match.length) Plotly.restyle(gd, {{'marker.opacity': 1.0}},  match);
        if (other.length) Plotly.restyle(gd, {{'marker.opacity': 0.04}}, other);
    }});

    gd.on('plotly_unhover', function() {{
        Plotly.restyle(gd, {{'marker.opacity': 0.85}}, gd.data.map((_, i) => i));
    }});
}}
"""

# =============================================================================
# ASSEMBLE HTML
# =============================================================================

plotly_div_html = fig.to_html(
    full_html=False,
    include_plotlyjs=True,
    div_id="main_fig",
    config={"responsive": True, "scrollZoom": True},
)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>scANVI Spatial Viewer — {CONFIG_NAME}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
      margin: 0; padding: 14px 16px;
      background: #ffffff; color: #222222;
    }}
    h2 {{
      margin: 0 0 10px 0; font-size: 14px; font-weight: 600;
      color: #333333; letter-spacing: 0.3px;
    }}
    .controls {{
      display: flex; flex-wrap: wrap; gap: 20px; align-items: flex-end;
      padding: 10px 14px; background: #f7f7f7;
      border: 1px solid #e0e0e0; border-radius: 6px; margin-bottom: 8px;
    }}
    .ctrl-group {{ display: flex; flex-direction: column; gap: 4px; }}
    .ctrl-label {{
      font-size: 10px; font-weight: 600; color: #666;
      text-transform: uppercase; letter-spacing: 0.6px;
    }}
    select {{
      background: #fff; color: #222; border: 1px solid #ccc;
      border-radius: 4px; padding: 5px 8px; font-size: 12px;
      cursor: pointer; min-width: 240px;
    }}
    select:hover, select:focus {{ border-color: #888; outline: none; }}
    .hint {{ font-size: 10px; color: #999; margin-top: 4px; }}
    #main_fig {{ width: 100%; }}
  </style>
</head>
<body>
  <h2>scANVI Predicted Cell Types — Slide-seq Spatial ({CONFIG_NAME})</h2>
  <div class="controls">
    <div class="ctrl-group">
      <span class="ctrl-label">Left Panel</span>
      <select id="selectLeft">{left_options}</select>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">Right Panel</span>
      <select id="selectRight">{right_options}</select>
    </div>
    <div class="ctrl-group">
      <span class="hint">{conf_note}</span>
      <span class="hint">Hover a dot to highlight that cell type across both panels</span>
    </div>
  </div>

  {plotly_div_html}

  <script>
    {js}
  </script>
</body>
</html>
"""

with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
    f.write(html)

print(f"\nSaved: {OUTPUT_HTML}")
print(f"  File size: {os.path.getsize(OUTPUT_HTML) / 1e6:.1f} MB")
