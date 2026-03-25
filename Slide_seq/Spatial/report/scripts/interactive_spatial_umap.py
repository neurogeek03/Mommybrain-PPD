import os
import re
import json
import argparse
import pandas as pd
import scanpy as sc
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ========== ARGS ==========
parser = argparse.ArgumentParser(
    description="Interactive Plotly HTML: neurons-only spatial sections + UMAP."
)
parser.add_argument(
    "-i", "--input_adata",
    default="/scratch/mfafouti/Mommybrain/Slide_seq/Spatial/report/data/objects/merged_neurons_non_neurons_slide_seq_15.h5ad",
)
parser.add_argument(
    "-m", "--metadata",
    default="/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR/slide_seq_metadata.csv",
)
parser.add_argument(
    "-c", "--color_map",
    default="/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv",
)
parser.add_argument(
    "-o", "--output",
    default="/scratch/mfafouti/Mommybrain/Slide_seq/Spatial/report/out_NN_seq/interactive_spatial_umap.html",
)
parser.add_argument("--umap_subsample", type=int, default=50000)
args = parser.parse_args()

os.makedirs(os.path.dirname(args.output), exist_ok=True)

# ========== COLORS ==========
color_df = pd.read_csv(args.color_map, usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)


def _num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except Exception:
        return -1


color_df["_num"] = color_df["name"].apply(_num_prefix)
color_df = color_df.sort_values("_num")
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

SPOT_CLASS_COLORS = {
    "singlet":  "#2166ac",
    "doublet":  "#d6604d",
    "rejected": "#aaaaaa",
}


def normalize(name):
    return re.sub(r"[ /-]", "_", str(name))


def infer_broad_class(name):
    n = normalize(name)
    if n.endswith("_NN"):
        return "Non-neuronal"
    if "_Glut" in n or "_IMN" in n:
        return "Glutamatergic"
    if "_Gaba" in n:
        return "GABAergic"
    return "Unknown"


def get_celltype_color(ct):
    return label_to_hex.get(ct, "#888888")


# ========== LOAD ADATA ==========
print(f"Loading: {args.input_adata}")
adata = sc.read_h5ad(args.input_adata)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  obs columns: {list(adata.obs.columns)}")
print(f"  obsm keys:   {list(adata.obsm.keys())}")


def find_col(obs, candidates):
    for c in candidates:
        if c in obs.columns:
            return c
    return None


TYPE_COL       = "RCTD_first_type_rat"
SPOT_CLASS_COL = find_col(adata.obs, ["RCTD_spot_class_rat", "RCTD_spot_class", "spot_class"])

print(f"  type_col:        {TYPE_COL}")
print(f"  spot_class_col:  {SPOT_CLASS_COL}")
print(f"  X_umap present:  {'X_umap' in adata.obsm}")

has_umap       = "X_umap" in adata.obsm
has_spot_class = SPOT_CLASS_COL is not None

assert TYPE_COL in adata.obs.columns, f"Missing required column: {TYPE_COL}"
assert "X_spatial" in adata.obsm,     "Missing obsm['X_spatial']"
assert "sample" in adata.obs.columns, "Missing obs['sample']"

# ========== BUILD MAIN DATAFRAME ==========
# Preserve original obs index throughout
obs_cols = [TYPE_COL, "sample"]
if SPOT_CLASS_COL:
    obs_cols.append(SPOT_CLASS_COL)

coords = adata.obsm["X_spatial"]
df = pd.DataFrame(coords, columns=["x", "y"], index=adata.obs_names)
df = df.join(adata.obs[obs_cols].copy(), how="left")
df["celltype_str"] = df[TYPE_COL].astype(str)
df["broad_class"]  = df["celltype_str"].apply(infer_broad_class)

# ---- Keep neurons only ----
df = df[df["broad_class"].isin({"Glutamatergic", "GABAergic"})].copy()
df["celltype_plot"] = df["celltype_str"].copy()
print(f"  Neurons retained: {len(df)}")

# ---- Join metadata (preserve index via map) ----
meta_df = pd.read_csv(args.metadata)
meta_df.columns = meta_df.columns.str.strip()
meta_indexed = meta_df.set_index("sample")
for col in meta_indexed.columns:
    df[col] = df["sample"].map(meta_indexed[col])

# ========== UMAP DATAFRAME ==========
umap_df = None
if has_umap:
    umap_coords = adata.obsm["X_umap"]
    umap_full = pd.DataFrame(umap_coords, columns=["umap_x", "umap_y"], index=adata.obs_names)
    # Filter to neurons only using the preserved index
    umap_df = umap_full.loc[umap_full.index.isin(df.index)].copy()
    umap_df["celltype_plot"] = df.loc[umap_df.index, "celltype_plot"]

    if len(umap_df) > args.umap_subsample:
        frac = args.umap_subsample / len(umap_df)
        umap_df = (
            umap_df.groupby("celltype_plot", group_keys=False)
            .apply(lambda x: x.sample(frac=frac, random_state=42))
        ).head(args.umap_subsample)
    print(f"  UMAP neurons after subsample: {len(umap_df)}")

# ========== SAMPLE METADATA ==========
available_samples = sorted(df["sample"].dropna().unique().tolist())
sample_label = {}
for sid in available_samples:
    row = meta_df[meta_df["sample"] == sid]
    if row.empty:
        sample_label[sid] = sid
    else:
        r = row.iloc[0]
        sample_label[sid] = (
            f"{sid} — {r['pregnancy']} | {r['day']} | "
            f"{r['treatment']} | {r['coronal_section']}"
        )
print(f"  Samples: {available_samples}")

# ========== SERIALIZE SPATIAL DATA AS RAW JSON ==========
# Data stored once per sample; JS builds Plotly traces on demand.
all_celltypes = sorted(df["celltype_plot"].dropna().unique().tolist())

sample_data = {}
for sid in available_samples:
    sdf = df[df["sample"] == sid]
    if sdf.empty:
        continue
    entry = {
        "x":       [round(v, 1) for v in sdf["x"].tolist()],
        "y":       [round(v, 1) for v in sdf["y"].tolist()],
        "ct":      sdf["celltype_plot"].tolist(),
    }
    if has_spot_class:
        entry["sc"] = sdf[SPOT_CLASS_COL].astype(str).replace("nan", "unknown").tolist()
    sample_data[sid] = entry

color_map = {ct: get_celltype_color(ct) for ct in all_celltypes}

# ========== BUILD FIGURE (UMAP TRACES ONLY) ==========
n_cols = 3 if has_umap else 2
col_widths = [0.36, 0.36, 0.28] if has_umap else [0.5, 0.5]
subplot_titles = ["Left Panel", "Right Panel"]
if has_umap:
    subplot_titles.append("UMAP (neurons)")

fig = make_subplots(
    rows=1, cols=n_cols,
    subplot_titles=subplot_titles,
    column_widths=col_widths,
    horizontal_spacing=0.03,
)

umap_trace_count = 0
if has_umap:
    for ct in all_celltypes:
        ct_umap = umap_df[umap_df["celltype_plot"] == ct]
        if ct_umap.empty:
            continue
        fig.add_trace(
            go.Scattergl(
                x=ct_umap["umap_x"].round(3).tolist(),
                y=ct_umap["umap_y"].round(3).tolist(),
                mode="markers",
                marker=dict(size=3, color=get_celltype_color(ct), opacity=0.7),
                name=ct,
                legendgroup=ct,
                showlegend=True,
                meta=ct,
                hovertemplate="%{meta}<extra></extra>",
            ),
            row=1, col=n_cols,
        )
        umap_trace_count += 1

print(f"  UMAP traces: {umap_trace_count}")

# ========== LAYOUT ==========
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
    margin=dict(l=10, r=170, t=55, b=10),
    height=620,
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
if has_umap:
    fig.update_xaxes(visible=False, showgrid=False, zeroline=False, row=1, col=n_cols)
    fig.update_yaxes(visible=False, showgrid=False, zeroline=False, row=1, col=n_cols)

for ann in fig.layout.annotations:
    ann.font.color = "#333333"
    ann.font.size  = 11

# ========== DROPDOWN HTML ==========
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

spotclass_btn = (
    '<button class="mode-btn" id="btn-spotclass" '
    'onclick="setMode(\'spotclass\')">Spot Class</button>'
    if has_spot_class else ""
)

# ========== JAVASCRIPT ==========
js = f"""
const sampleData     = {json.dumps(sample_data)};
const colorMap       = {json.dumps(color_map)};
const spotClassColors = {json.dumps(SPOT_CLASS_COLORS)};
const sampleMeta     = {json.dumps(sample_label)};
const umap_trace_count = {umap_trace_count};
const has_umap       = {'true' if has_umap else 'false'};
const n_cols         = {n_cols};

let currentLeft  = '{default_left}';
let currentRight = '{default_right}';
let currentMode  = 'celltype';
let umapTracesData = null;   // deep-copied on first render

const gd = document.getElementById('main_fig');

// Capture UMAP traces once after initial render, then load defaults
gd.on('plotly_afterplot', function initView() {{
    gd.removeListener('plotly_afterplot', initView);
    umapTracesData = JSON.parse(JSON.stringify(gd.data));
    updateView();
    attachHoverEvents();
}});

// ---- Build spatial traces from raw data in JS ----
function buildSpatialTraces(sampleId, col) {{
    if (!sampleId || !sampleData[sampleId]) return [];
    const d = sampleData[sampleId];
    const n = d.x.length;
    const xax = col === 1 ? 'x'  : ('x' + col);
    const yax = col === 1 ? 'y'  : ('y' + col);
    const traces = [];

    if (currentMode === 'celltype') {{
        // Group by cell type
        const groups = {{}};
        for (let i = 0; i < n; i++) {{
            const ct = d.ct[i];
            if (!groups[ct]) groups[ct] = {{x: [], y: []}};
            groups[ct].x.push(d.x[i]);
            groups[ct].y.push(d.y[i]);
        }}
        for (const [ct, pts] of Object.entries(groups)) {{
            traces.push({{
                type: 'scattergl',
                x: pts.x, y: pts.y,
                mode: 'markers',
                marker: {{size: 4, color: colorMap[ct] || '#888', opacity: 0.8}},
                name: ct,
                legendgroup: ct,
                showlegend: false,
                meta: ct,
                xaxis: xax, yaxis: yax,
                hovertemplate: sampleId + ' | ' + ct + '<extra></extra>',
            }});
        }}
    }} else if (currentMode === 'spotclass' && d.sc) {{
        // Group by spot class
        const groups = {{}};
        for (let i = 0; i < n; i++) {{
            const sc = d.sc[i] || 'unknown';
            if (!groups[sc]) groups[sc] = {{x: [], y: []}};
            groups[sc].x.push(d.x[i]);
            groups[sc].y.push(d.y[i]);
        }}
        for (const [sc, pts] of Object.entries(groups)) {{
            traces.push({{
                type: 'scattergl',
                x: pts.x, y: pts.y,
                mode: 'markers',
                marker: {{size: 4, color: spotClassColors[sc] || '#888', opacity: 0.8}},
                name: sc,
                legendgroup: 'sc__' + sc,
                showlegend: col === 1,   // show legend entry only from left panel
                meta: 'sc__' + sc,
                xaxis: xax, yaxis: yax,
                hovertemplate: sampleId + ' | ' + sc + '<extra></extra>',
            }});
        }}
    }}
    return traces;
}}

// ---- Main update: rebuild all non-UMAP traces ----
function updateView() {{
    if (!umapTracesData) return;

    const leftTraces  = buildSpatialTraces(currentLeft,  1);
    const rightTraces = buildSpatialTraces(currentRight, 2);

    // In spotclass mode hide UMAP (no matching legend); in celltype keep it
    const umap = (currentMode === 'celltype') ? umapTracesData : [];
    const allTraces = [...umap, ...leftTraces, ...rightTraces];

    // Preserve axis layout (zoom, ranges)
    Plotly.react(gd, allTraces, gd.layout);

    // Update panel titles
    const leftLabel  = sampleMeta[currentLeft]  || currentLeft;
    const rightLabel = sampleMeta[currentRight] || currentRight;
    const titleUpdates = {{'annotations[0].text': leftLabel, 'annotations[1].text': rightLabel}};
    if (has_umap) titleUpdates['annotations[2].text'] = 'UMAP (neurons)';
    Plotly.relayout(gd, titleUpdates);
}}

function setMode(mode) {{
    currentMode = mode;
    document.querySelectorAll('.mode-btn').forEach(b => b.classList.remove('active'));
    const btn = document.getElementById('btn-' + mode);
    if (btn) btn.classList.add('active');
    updateView();
}}

document.getElementById('selectLeft').addEventListener('change',  function() {{
    currentLeft = this.value;
    updateView();
}});
document.getElementById('selectRight').addEventListener('change', function() {{
    currentRight = this.value;
    updateView();
}});

// ---- Hover-to-highlight (cell type mode only) ----
function attachHoverEvents() {{
    gd.on('plotly_hover', function(eventData) {{
        if (currentMode !== 'celltype') return;
        const hovered = eventData.points[0].data.meta;
        if (!hovered || hovered.startsWith('sc__')) return;

        const match = [], other = [];
        gd.data.forEach(function(trace, i) {{
            if (trace.meta === hovered) match.push(i);
            else other.push(i);
        }});
        if (match.length) Plotly.restyle(gd, {{'marker.opacity': 1.0}},  match);
        if (other.length) Plotly.restyle(gd, {{'marker.opacity': 0.04}}, other);
    }});

    gd.on('plotly_unhover', function() {{
        if (currentMode !== 'celltype') return;
        Plotly.restyle(gd, {{'marker.opacity': 0.8}},
            gd.data.map((_, i) => i));
    }});
}}
"""

# ========== RENDER & ASSEMBLE HTML ==========
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
  <title>Slide-seq Spatial Viewer</title>
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
      cursor: pointer; min-width: 260px;
    }}
    select:hover, select:focus {{ border-color: #888; outline: none; }}
    .btn-row {{ display: flex; gap: 5px; }}
    .mode-btn {{
      background: #fff; color: #444; border: 1px solid #ccc;
      border-radius: 4px; padding: 5px 11px; font-size: 11px; cursor: pointer;
      transition: background 0.1s;
    }}
    .mode-btn:hover {{ background: #f0f0f0; border-color: #999; }}
    .mode-btn.active {{
      background: #222; color: #fff; border-color: #222; font-weight: 600;
    }}
    .hint {{ font-size: 10px; color: #999; margin-top: 4px; }}
    #main_fig {{ width: 100%; }}
  </style>
</head>
<body>
  <h2>Slide-seq Spatial Viewer &mdash; Neurons</h2>
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
      <span class="ctrl-label">Color by</span>
      <div class="btn-row">
        <button class="mode-btn active" id="btn-celltype" onclick="setMode('celltype')">Cell Type</button>
        {spotclass_btn}
      </div>
      <span class="hint">Hover a dot to highlight that cell type across all panels</span>
    </div>
  </div>

  {plotly_div_html}

  <script>
    {js}
  </script>
</body>
</html>
"""

with open(args.output, "w", encoding="utf-8") as f:
    f.write(html)

print(f"\nSaved: {args.output}")
print(f"  File size: {os.path.getsize(args.output) / 1e6:.1f} MB")
