from multiprocessing import parent_process
from re import T
import pandas as pd 
from pathlib import Path 
import tqdm
import re 

project_dir = Path.cwd().parents[0]

seq_dir = project_dir / 'seq'
tags_dir = project_dir/ 'tags' / 'CORT_vs_OIL'
out_path = project_dir / 'out'
out_path.mkdir(parents=True, exist_ok=True)


for tags_file in tags_dir.glob("*.tsv"):
    celltype = tags_file.stem.split("_edgeR")[0]
    print(f'Found celltype {celltype}, looking for matches in Slide-seq data...')
    
    tags_filepath = tags_file
    tags_df = pd.read_csv(tags_filepath, sep="\t", index_col=None, skiprows=1)
    tags_df.columns = ['gene','logFC','logCPM','LR','PValue','FDR']
    #print(tags_df.head())
    
    for dir in seq_dir.iterdir():
        for seq_file in dir.glob("*.tsv"):
            if re.match(celltype, seq_file.stem):
                print('Found a match!')
                seq_filepath = dir / seq_file
                seq_df = pd.read_csv(seq_filepath, sep="\t", index_col=None, skiprows=1)
                seq_df.columns = ['gene','logFC','logCPM','LR','PValue','FDR']
                print('Merging...')
                merged = seq_df.merge(tags_df, on = 'gene')
                print(merged.head())
                
                print('filtering merged...')
                filtered = merged[
                    (merged['logFC_x'].abs() > 0.1) &
                    (merged['FDR_x'] < 0.05)
                ]
                print (filtered.head())
                
                comparison_dir = out_path / dir.name
                comparison_dir.mkdir(parents=True, exist_ok=True)

                if not filtered.empty:
                    filtered_path = comparison_dir / f'filtered{celltype}_tags_seq_merged.tsv'
                    filtered.to_csv(filtered_path, sep="\t") 
                    print("Filtered results saved.")
                else:
                    print("No rows passed the filter. Nothing saved.")
                