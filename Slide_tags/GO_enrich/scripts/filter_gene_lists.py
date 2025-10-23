import pandas as pd
from pathlib import Path 

project_dir = Path.cwd().parents[0]
out_dir = project_dir / 'out' / 'gene_lists'
out_dir.mkdir(parents=True, exist_ok=True)
data_folder = project_dir / 'EdgeR'
data_dir = data_folder / 'out' / 'edger_lrt'
print(data_dir)

for tags_file in data_dir.glob('*.tsv'):
    print(tags_file)
    celltype = tags_file.stem.split("_edgeR")[0]
    print(f'Found celltype {celltype}...')
    tags_df = pd.read_csv(tags_file, sep="\t", index_col=0)
    print(tags_df.head())
    filtered_UP = tags_df[
                    (tags_df['logFC'] > 0.1) &
                    (tags_df['FDR'] < 0.1)
                ]
    genes_up = filtered_UP.index
    out_path_up = out_dir / f'new_UP_{celltype}_genes_tags.csv'
    genes_up.to_series().to_csv(out_path_up, index=False, header=False)

    filtered_DOWN = tags_df[
                    (tags_df['logFC'] < 0.1) &
                    (tags_df['FDR'] < 0.1)
                ]
    genes_down = filtered_DOWN.index
    out_path_down = out_dir / f'new_DOWN_{celltype}_genes_tags.csv'
    genes_down.to_series().to_csv(out_path_down, index=False, header=False)