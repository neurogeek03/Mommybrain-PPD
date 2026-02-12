
import pandas as pd
import os
import plotly.graph_objects as go
from pathlib import Path
import numpy as np

def generate_heatmap(dataframe, filename_suffix, output_base_path):
    """
    Generates and saves an interactive heatmap and a static PNG for the given dataframe.
    Returns the generated heatmap dataframe of RAW counts.
    """
    
    print(f"\n--- Generating heatmap for '{filename_suffix}' group ---")
    
    # Filter out rejects
    plot_df = dataframe[dataframe['spot_class'] != 'reject'].copy()
    print(f"Using {len(plot_df)} rows after filtering rejects.")

    if not plot_df.empty and 'first_type' in plot_df.columns and 'second_type' in plot_df.columns:
        plot_df.dropna(subset=['first_type', 'second_type'], inplace=True)
        
        if plot_df.empty:
            print("No data left after dropping NaNs. Skipping plot generation.")
            return None

        print(f"Creating heatmap from {len(plot_df)} data points.")
        links_df = plot_df.groupby(['first_type', 'second_type']).size().reset_index(name='value')
        
        print(f"Found {len(links_df)} unique connections.")
        print("Top 5 connections (raw counts):")
        print(links_df.nlargest(5, 'value'))

        # Create the raw count matrix to be returned for the diff plot
        raw_heatmap_df = links_df.pivot(index='first_type', columns='second_type', values='value').fillna(0)

        # --- Log Transform for Plotting ---
        links_df['log_value'] = np.log1p(links_df['value'])
        log_heatmap_df = links_df.pivot(index='first_type', columns='second_type', values='log_value').fillna(0)
        
        # Create dataframe for hover text (raw counts)
        hover_df = links_df.pivot(index='first_type', columns='second_type', values='value').fillna(0)
        hover_df = hover_df.reindex(index=log_heatmap_df.index, columns=log_heatmap_df.columns, fill_value=0)
        text_df = hover_df.applymap(lambda count: f'{count:.0f}')

        print(f"Heatmap dimensions: {log_heatmap_df.shape[0]} first_types x {log_heatmap_df.shape[1]} second_types.")
        
        # Create and save the figure
        fig = go.Figure(data=go.Heatmap(
                           z=log_heatmap_df.values,
                           x=log_heatmap_df.columns,
                           y=log_heatmap_df.index,
                           colorscale='Viridis',
                           hoverongaps=False,
                           text=text_df.values, # Use text property
                           hovertemplate='<b>First Type</b>: %{y}<br>' +
                                         '<b>Second Type</b>: %{x}<br>' +
                                         '<b>Raw Count</b>: %{text}<extra></extra>')) # Use %{text}

        fig.update_layout(
            title=f'Log-Transformed Heatmap of Connections ({filename_suffix})',
            xaxis_title="Second Type", yaxis_title="First Type",
            xaxis=dict(tickangle=-45))

        html_plot_file = output_base_path / f"connections_heatmap_log_{filename_suffix}.html"
        fig.write_html(html_plot_file, include_plotlyjs=True)
        print(f"Successfully created log heatmap and saved to {html_plot_file}")
        
        try:
            png_plot_file = output_base_path / f"connections_heatmap_log_{filename_suffix}.png"
            fig.write_image(png_plot_file, width=2400, height=1600, scale=1)
            print(f"Successfully created log PNG heatmap and saved to {png_plot_file}")
        except (ValueError, RuntimeError) as e:
            print(f"\n--- Could not create PNG for {filename_suffix} ---\n{e}\n---")

        return raw_heatmap_df # Return the raw counts for diff calculation
    else:
        print(f"Could not create plot for '{filename_suffix}'. No data or required columns missing.")
        return None

def merge_metadata():
    # Define paths
    project_path = Path.cwd().parents[0]
    output_base = project_path / 'out'
    output_base.mkdir(exist_ok=True, parents=True)
    csv_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/OUT/FINAL_RCTD_newgenelist/merged_metadata_csvs"
    metadata_file = "/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR/slide_seq_metadata.csv"
    output_file = output_base / "merged_slide_seq_metadata.csv"

    # Read and process data
    df_list = [pd.read_csv(os.path.join(csv_dir, f)).assign(sample=f.split('_')[0]) for f in os.listdir(csv_dir) if f.endswith(".csv")]
    merged_df = pd.concat(df_list, ignore_index=True)
    metadata_df = pd.read_csv(metadata_file)
    final_df = pd.merge(merged_df, metadata_df, on='sample', how='left').rename(columns={'treatment': 'condition'})
    final_df.to_csv(output_file, index=False)
    
    print(f"Successfully merged {len(df_list)} files and saved to {output_file}")
    print("\nValue counts of 'condition' column:\n", final_df['condition'].value_counts())

    # --- Data Splitting and Plotting ---
    final_df['condition'] = final_df['condition'].astype(str)
    cort_mask = final_df['condition'].str.contains('CORT', na=False)
    df_cort = final_df[cort_mask]
    df_other = final_df[~cort_mask]

    print(f"\nSplitting data into two groups:\n1. 'CORT' group: {len(df_cort)} rows\n2. 'OTHER' group: {len(df_other)} rows")

    # Generate a heatmap for each group and get the data back
    heatmap_cort = generate_heatmap(df_cort, 'CORT', output_base)
    heatmap_other = generate_heatmap(df_other, 'OTHER', output_base)

    # --- Differential Heatmap Logic ---
    if heatmap_cort is not None and heatmap_other is not None:
        print("\n--- Generating differential heatmap (OTHER - CORT) ---")

        common_rows = heatmap_cort.index.intersection(heatmap_other.index)
        common_cols = heatmap_cort.columns.intersection(heatmap_other.columns)
        print(f"Found {len(common_rows)} common first_types and {len(common_cols)} common second_types for diff plot.")

        if len(common_rows) > 0 and len(common_cols) > 0:
            cort_reindexed = heatmap_cort.reindex(index=common_rows, columns=common_cols, fill_value=0)
            other_reindexed = heatmap_other.reindex(index=common_rows, columns=common_cols, fill_value=0)
            diff_df = other_reindexed - cort_reindexed
            
            print("Top 5 largest increases in 'OTHER' group:\n", diff_df.stack().nlargest(5))
            print("\nTop 5 largest decreases in 'OTHER' group (i.e., increases in 'CORT'):\n", diff_df.stack().nsmallest(5))

            fig_diff = go.Figure(data=go.Heatmap(
                                   z=diff_df.values, x=diff_df.columns, y=diff_df.index,
                                   colorscale='RdBu', zmid=0))

            fig_diff.update_layout(
                title='Differential Heatmap (OTHER - CORT)',
                xaxis_title="Second Type", yaxis_title="First Type",
                xaxis=dict(tickangle=-45))
            
            diff_html_file = output_base / "connections_heatmap_DIFF.html"
            fig_diff.write_html(diff_html_file, include_plotlyjs=True)
            print(f"\nSuccessfully created differential heatmap and saved to {diff_html_file}")
            
            try:
                diff_png_file = output_base / "connections_heatmap_DIFF.png"
                fig_diff.write_image(diff_png_file, width=2400, height=1600, scale=1)
                print(f"Successfully created differential PNG and saved to {diff_png_file}")
            except (ValueError, RuntimeError) as e:
                print(f"\n--- Could not create PNG for DIFF ---\n{e}\n---")
        else:
            print("No common connections found between CORT and OTHER groups. Skipping differential plot.")

if __name__ == "__main__":
    merge_metadata()
