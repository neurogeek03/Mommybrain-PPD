"""
HTML builder for the spatial gene expression interactive viewer.
Produces a fully standalone HTML file — no server required.
"""

import json
import plotly.io as pio


def build_html(fig, expr_data, celltype_data, genes, celltypes, sample_names,
               output_path, point_size=2, global_min=0.0, global_max=1.0,
               celltype_color_map=None):
    """
    Build and write a standalone HTML file.

    Parameters
    ----------
    fig             : plotly Figure — used for initial layout/render
    expr_data       : {gene: {sample: [float list]}} — normalized expression values
    celltype_data   : {sample: [str list]} — cell type per cell, per sample
    genes           : list of gene symbol strings
    celltypes       : sorted list of all cell type strings
    sample_names    : list of 3 sample name strings (order matches figure traces)
    output_path     : str — path to write the HTML file
    point_size      : int
    global_min      : float — global expression min for colorscale
    global_max      : float — global expression max for colorscale
    celltype_color_map : {celltype: "#hex"} — categorical color assignments
    """
    if celltype_color_map is None:
        celltype_color_map = {}

    # Round expression to 4 dp to limit file size
    expr_rounded = {
        gene: {
            sample: [round(v, 4) for v in vals]
            for sample, vals in sd.items()
        }
        for gene, sd in expr_data.items()
    }

    # Serialize JS data blobs
    genes_js      = json.dumps(genes)
    samples_js    = json.dumps(sample_names)
    expr_js       = json.dumps(expr_rounded)
    celltype_js   = json.dumps(celltype_data)
    celltypes_js  = json.dumps(celltypes)
    ct_colors_js  = json.dumps(celltype_color_map)
    fig_json      = pio.to_json(fig)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Spatial Gene Expression Viewer</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ background: #1a1a1a; color: #e0e0e0; font-family: 'Segoe UI', Helvetica, sans-serif; }}
    #controls {{
      display: flex;
      align-items: center;
      gap: 24px;
      padding: 12px 20px;
      background: #2a2a2a;
      border-bottom: 1px solid #444;
      flex-wrap: wrap;
    }}
    .ctrl-group {{ display: flex; align-items: center; gap: 8px; }}
    .ctrl-label {{ font-size: 12px; color: #aaa; white-space: nowrap; text-transform: uppercase; letter-spacing: 0.05em; }}
    select {{
      background: #333;
      color: #e0e0e0;
      border: 1px solid #555;
      border-radius: 4px;
      padding: 5px 10px;
      font-size: 13px;
      min-width: 220px;
      cursor: pointer;
    }}
    select:focus {{ outline: none; border-color: #7ec8e3; }}
    .radio-group {{ display: flex; gap: 14px; }}
    .radio-group label {{ color: #e0e0e0; font-size: 13px; cursor: pointer; display: flex; align-items: center; gap: 4px; }}
    input[type="radio"] {{ cursor: pointer; accent-color: #7ec8e3; }}
    #plotly-div {{ width: 100%; }}
  </style>
</head>
<body>
  <div id="controls">
    <div class="ctrl-group">
      <span class="ctrl-label">Gene</span>
      <select id="gene-select"></select>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">Color by</span>
      <div class="radio-group">
        <label><input type="radio" name="colormode" value="expression" checked> Expression</label>
        <label><input type="radio" name="colormode" value="celltype"> Cell Type</label>
      </div>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">Filter cell type</span>
      <select id="celltype-select"></select>
    </div>
  </div>
  <div id="plotly-div"></div>

  <script>
    // === Embedded data ===
    var GENES             = {genes_js};
    var SAMPLES           = {samples_js};
    var EXPR_DATA         = {expr_js};
    var CELLTYPE_DATA     = {celltype_js};
    var CELLTYPES         = {celltypes_js};
    var CELLTYPE_COLOR_MAP = {ct_colors_js};
    var GLOBAL_MIN        = {global_min};
    var GLOBAL_MAX        = {global_max};
    var POINT_SIZE        = {point_size};

    // === Populate gene dropdown ===
    var geneSelect = document.getElementById('gene-select');
    GENES.forEach(function(g) {{
      var opt = document.createElement('option');
      opt.value = g;
      opt.textContent = g;
      geneSelect.appendChild(opt);
    }});

    // === Populate cell type dropdown ===
    var ctSelect = document.getElementById('celltype-select');
    var allOpt = document.createElement('option');
    allOpt.value = 'All';
    allOpt.textContent = 'All cell types';
    ctSelect.appendChild(allOpt);
    CELLTYPES.forEach(function(ct) {{
      var opt = document.createElement('option');
      opt.value = ct;
      opt.textContent = ct;
      ctSelect.appendChild(opt);
    }});

    // === Initial render from embedded figure JSON ===
    var figData = {fig_json};
    Plotly.newPlot('plotly-div', figData.data, figData.layout, {{responsive: true}});

    // === Core update function ===
    // Non-selected cells are set to null so Plotly renders them as transparent.
    // This avoids re-sending x/y coordinates and keeps the file compact.
    function updatePlot() {{
      var gene      = geneSelect.value;
      var colorMode = document.querySelector('input[name="colormode"]:checked').value;
      var filterCT  = ctSelect.value;

      for (var i = 0; i < SAMPLES.length; i++) {{
        var sample  = SAMPLES[i];
        var ctArr   = CELLTYPE_DATA[sample];
        var exprArr = EXPR_DATA[gene][sample];
        var isLast  = (i === SAMPLES.length - 1);
        var colors;

        if (colorMode === 'expression') {{
          colors = exprArr.map(function(val, idx) {{
            return (filterCT === 'All' || ctArr[idx] === filterCT) ? val : null;
          }});
        }} else {{
          colors = ctArr.map(function(ct, idx) {{
            return (filterCT === 'All' || ct === filterCT)
              ? (CELLTYPE_COLOR_MAP[ct] || '#888888')
              : null;
          }});
        }}

        Plotly.restyle('plotly-div', {{
          'marker.color':      [colors],
          'marker.colorscale': [colorMode === 'expression' ? 'Viridis' : null],
          'marker.cmin':       [colorMode === 'expression' ? GLOBAL_MIN : null],
          'marker.cmax':       [colorMode === 'expression' ? GLOBAL_MAX : null],
          'marker.showscale':  [isLast && colorMode === 'expression'],
        }}, [i]);
      }}
    }}

    // === Event listeners ===
    geneSelect.addEventListener('change', updatePlot);
    ctSelect.addEventListener('change', updatePlot);
    document.querySelectorAll('input[name="colormode"]').forEach(function(radio) {{
      radio.addEventListener('change', updatePlot);
    }});
  </script>
</body>
</html>"""

    with open(output_path, 'w') as f:
        f.write(html)

    print(f"Saved: {output_path}")