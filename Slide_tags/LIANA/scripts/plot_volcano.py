import pandas as pd
import numpy as np
import json
from pathlib import Path

# Paths (relative to scripts/ directory, same pattern as CCC_liana.py)
base_dir = Path.cwd().parents[0]
dea_path = base_dir / "out" / "first_try" / "dea_result.csv"
out_dir = base_dir / "out" / "figures" / "volcano_plots"
out_dir.mkdir(parents=True, exist_ok=True)

# Load DEA results
df = pd.read_csv(dea_path, index_col=0)
df = df.reset_index().rename(columns={"index": "gene"})

epsilon = 1e-300
df["-log10(padj)"] = -np.log10(df["padj"].clip(lower=epsilon))

cell_types = sorted(df["subclass_name"].unique().tolist())
print(f"Building combined HTML for {len(cell_types)} cell types...")

# Build per-cell-type data dict for embedding in HTML
# NaN padj -> 1.0 so JSON serialization and JS comparisons are safe
plot_data = {}
for ct in cell_types:
    sub = df[df["subclass_name"] == ct].copy()
    sub["padj"] = sub["padj"].fillna(1.0)
    plot_data[ct] = {
        "gene":           sub["gene"].tolist(),
        "x":              sub["log2FoldChange"].round(4).tolist(),
        "y":              sub["-log10(padj)"].round(4).tolist(),
        "padj":           sub["padj"].tolist(),
        "baseMean":       sub["baseMean"].round(2).tolist(),
    }

data_json    = json.dumps(plot_data)
ct_list_json = json.dumps(cell_types)

html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Volcano Plots — CORT vs OIL</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
    * {{ box-sizing: border-box; }}
    body {{
      font-family: Arial, sans-serif;
      margin: 0;
      padding: 16px 24px;
      background: #f4f6f8;
      color: #222;
    }}
    h1 {{ font-size: 1.3em; margin: 0 0 14px; }}
    .controls {{
      display: flex;
      flex-wrap: wrap;
      gap: 24px;
      align-items: center;
      background: white;
      border: 1px solid #ddd;
      border-radius: 6px;
      padding: 12px 16px;
      margin-bottom: 14px;
    }}
    .ctrl-group {{ display: flex; align-items: center; gap: 8px; }}
    .ctrl-group label.lbl {{ font-weight: bold; white-space: nowrap; }}
    select {{
      font-size: 0.95em;
      padding: 4px 8px;
      border: 1px solid #bbb;
      border-radius: 4px;
      background: white;
      cursor: pointer;
    }}
    .radio-row {{ display: flex; gap: 14px; }}
    .radio-row label {{
      display: flex;
      align-items: center;
      gap: 4px;
      cursor: pointer;
      font-size: 0.95em;
    }}
    #sig-count {{
      font-size: 0.9em;
      color: #555;
      margin-left: auto;
    }}
    #volcano-plot {{
      background: white;
      border: 1px solid #ddd;
      border-radius: 6px;
      width: 100%;
      height: 620px;
    }}
  </style>
</head>
<body>
  <h1>Volcano Plots — CORT vs OIL (DESeq2)</h1>
  <div class="controls">
    <div class="ctrl-group">
      <label class="lbl" for="ct-select">Cell type:</label>
      <select id="ct-select"></select>
    </div>
    <div class="ctrl-group">
      <span class="lbl">FDR threshold:</span>
      <div class="radio-row">
        <label><input type="radio" name="fdr" value="0.05" checked> 0.05</label>
        <label><input type="radio" name="fdr" value="0.1"> 0.10</label>
        <label><input type="radio" name="fdr" value="0.5"> 0.50</label>
      </div>
    </div>
    <div class="ctrl-group">
      <span class="lbl">|log2FC| threshold:</span>
      <div class="radio-row">
        <label><input type="radio" name="lfc" value="0.1"> 0.1</label>
        <label><input type="radio" name="lfc" value="1.0" checked> 1.0</label>
      </div>
    </div>
    <span id="sig-count"></span>
  </div>
  <div id="volcano-plot"></div>

  <script>
    const allData   = {data_json};
    const cellTypes = {ct_list_json};

    // Populate cell type dropdown
    const sel = document.getElementById('ct-select');
    cellTypes.forEach(ct => {{
      const opt = document.createElement('option');
      opt.value = ct;
      opt.textContent = ct;
      sel.appendChild(opt);
    }});

    function getSelectedFdr() {{
      return parseFloat(document.querySelector('input[name=fdr]:checked').value);
    }}

    function getSelectedLfc() {{
      return parseFloat(document.querySelector('input[name=lfc]:checked').value);
    }}

    function drawPlot() {{
      const ct     = sel.value;
      const fdr    = getSelectedFdr();
      const lfc    = getSelectedLfc();
      const d      = allData[ct];
      const hlineY = -Math.log10(fdr);

      const sigX = [], sigY = [], sigText = [];
      const nsX  = [], nsY  = [], nsText  = [];

      for (let i = 0; i < d.gene.length; i++) {{
        const isSig = (d.padj[i] < fdr) && (Math.abs(d.x[i]) > lfc);
        const tip = `<b>${{d.gene[i]}}</b><br>`
                  + `log2FC: ${{d.x[i].toFixed(3)}}<br>`
                  + `padj: ${{d.padj[i].toExponential(2)}}<br>`
                  + `baseMean: ${{d.baseMean[i].toFixed(1)}}`;
        if (isSig) {{
          sigX.push(d.x[i]); sigY.push(d.y[i]); sigText.push(tip);
        }} else {{
          nsX.push(d.x[i]); nsY.push(d.y[i]); nsText.push(tip);
        }}
      }}

      document.getElementById('sig-count').textContent =
        `${{sigX.length}} significant / ${{d.gene.length}} genes`;

      const traces = [
        {{
          x: nsX, y: nsY,
          mode: 'markers', type: 'scatter',
          name: 'NS',
          marker: {{ color: 'lightgrey', size: 5, opacity: 0.6 }},
          text: nsText,
          hovertemplate: '%{{text}}<extra></extra>'
        }},
        {{
          x: sigX, y: sigY,
          mode: 'markers', type: 'scatter',
          name: `Significant`,
          marker: {{ color: 'crimson', size: 6, opacity: 0.85 }},
          text: sigText,
          hovertemplate: '%{{text}}<extra></extra>'
        }}
      ];

      const layout = {{
        title: {{ text: ct, font: {{ size: 15 }} }},
        xaxis: {{ title: 'log2 Fold Change (CORT / OIL)', zeroline: true, zerolinecolor: '#aaa' }},
        yaxis: {{ title: '-log10(padj)' }},
        shapes: [
          {{
            type: 'line', xref: 'x', yref: 'paper',
            x0: -lfc, x1: -lfc, y0: 0, y1: 1,
            line: {{ dash: 'dash', color: '#444', width: 1 }}
          }},
          {{
            type: 'line', xref: 'x', yref: 'paper',
            x0: lfc, x1: lfc, y0: 0, y1: 1,
            line: {{ dash: 'dash', color: '#444', width: 1 }}
          }},
          {{
            type: 'line', xref: 'paper', yref: 'y',
            x0: 0, x1: 1, y0: hlineY, y1: hlineY,
            line: {{ dash: 'dash', color: '#444', width: 1 }}
          }}
        ],
        annotations: [{{
          xref: 'paper', yref: 'y',
          x: 1.01, y: hlineY,
          text: `FDR ${{fdr}}`,
          showarrow: false,
          font: {{ size: 11, color: '#444' }},
          xanchor: 'left'
        }}],
        legend: {{ orientation: 'v', x: 1.08, y: 1 }},
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
        margin: {{ t: 50, b: 60, l: 70, r: 110 }}
      }};

      Plotly.react('volcano-plot', traces, layout);
    }}

    sel.addEventListener('change', drawPlot);
    document.querySelectorAll('input[name=fdr]').forEach(r =>
      r.addEventListener('change', drawPlot)
    );
    document.querySelectorAll('input[name=lfc]').forEach(r =>
      r.addEventListener('change', drawPlot)
    );

    drawPlot();
  </script>
</body>
</html>
"""

out_path = out_dir / "volcano_plots_all.html"
out_path.write_text(html_content)
print(f"Done. Single combined HTML saved to:\n  {out_path}")
