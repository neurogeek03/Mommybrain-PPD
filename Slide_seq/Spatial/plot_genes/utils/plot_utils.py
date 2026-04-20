"""
Plotly figure construction for the spatial gene expression Dash viewer.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def build_figure(sample_names, spatial_data, color_data, hover_celltypes,
                 colorscale="Viridis", cmin=None, cmax=None,
                 show_colorbar=True, point_size=2):
    """
    Build a 1-row spatial scatter figure.

    Parameters
    ----------
    sample_names    : list of str — one per subplot
    spatial_data    : {sample: {"x": [...], "y": [...]}}
    color_data      : {sample: list} — floats for expression, hex strings for cell type,
                      None elements are rendered transparent by Plotly
    hover_celltypes : {sample: [str]} — cell type label per cell (for hover text)
    colorscale      : str — Plotly colorscale name (used only when colors are numeric)
    cmin, cmax      : float or None — colorscale range
    show_colorbar   : bool — show colorbar on the last subplot
    point_size      : int

    Returns
    -------
    plotly.graph_objects.Figure
    """
    n = len(sample_names)
    fig = make_subplots(
        rows=1,
        cols=n,
        subplot_titles=sample_names,
        horizontal_spacing=0.05,
    )

    for i, sample in enumerate(sample_names):
        is_last = (i == n - 1)
        colors = color_data[sample]
        celltypes = hover_celltypes[sample]

        # Determine if colors are numeric (expression) or categorical (cell type)
        numeric_colors = isinstance(colors[0], (int, float)) if colors else True

        marker = dict(
            color=colors,
            size=point_size,
            colorscale=colorscale if numeric_colors else None,
            cmin=cmin if numeric_colors else None,
            cmax=cmax if numeric_colors else None,
            showscale=is_last and show_colorbar and numeric_colors,
            colorbar=dict(x=1.02, title="log1p") if (is_last and numeric_colors) else None,
        )

        fig.add_trace(
            go.Scattergl(
                x=spatial_data[sample]["x"],
                y=spatial_data[sample]["y"],
                mode="markers",
                marker=marker,
                text=celltypes,
                hovertemplate="<b>%{text}</b><extra></extra>",
                showlegend=False,
                name=sample,
            ),
            row=1, col=i + 1,
        )

    fig.update_layout(
        plot_bgcolor="black",
        paper_bgcolor="#1a1a1a",
        font=dict(color="white"),
        height=580,
        margin=dict(t=50, b=10, l=10, r=90),
    )

    for i in range(1, n + 1):
        suffix = "" if i == 1 else str(i)
        fig.update_layout(**{
            f"xaxis{suffix}": dict(showticklabels=False, showgrid=False, zeroline=False,
                                   scaleanchor=f"y{suffix}", scaleratio=1),
            f"yaxis{suffix}": dict(showticklabels=False, showgrid=False, zeroline=False),
        })

    return fig
