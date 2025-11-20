# functions used for CAST 
###########################################################################
# Preview first 5 rows of the experession matrix: 
import numpy as np

def preview_adata_X(adata, n=5):
    # Extract first n rows
    X = adata.X[:n]
    if hasattr(X, "toarray"):
        X = X.toarray()
    
    # Print the preview
    print(f"=== Preview of adata.X (first {n} rows) ===")
    print(X)
    
    arr = X.toarray() if hasattr(X, "toarray") else X
    # Flatten for diagnostics
    flat = arr.ravel()
    nonzero = flat[flat != 0]

    print("\n=== Quick diagnostics ===")

    if len(nonzero) == 0:
        print("No nonzero values detected in these rows.")
    else:
        print(f"Nonzero entries in preview: {len(nonzero)}")
        print("First few nonzero values:", nonzero[:20])

        integer_like = np.allclose(nonzero, np.round(nonzero))
        print("Values appear integer-like:", bool(integer_like))

        if integer_like:
            print("→ Data looks like raw counts.")
        else:
            print("→ Data does not look like raw counts.")

        print(f"Range of nonzero values: min={nonzero.min()}, max={nonzero.max()}")

    # Check if all values are integer-like
    if np.all(np.floor(X) == X):
        print("Values appear integer-like: True → Data is raw")
    else:
        print("Values appear integer-like: False → Data may be normalized or processed")

###########################################################################
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def plot_brains(
    csv_path,
    group_col,
    point_size,
    x_col="x",
    y_col="y",
    color_col="cluster_alias",
    outdir=Path("out/plots")
    ):
    # Load the data
    df = pd.read_csv(csv_path)

    # Create output directory
    outdir.mkdir(exist_ok=True, parents=True)

    # Group the data
    grouped = df.groupby(group_col)

    # Loop through each group
    for group_value, group_df in grouped:
        plt.figure(figsize=(6, 6))

        # Scatter x vs y, colored by cluster_alias
        scatter = plt.scatter(
            group_df[x_col],
            group_df[y_col],
            c=group_df[color_col].astype("category").cat.codes,
            s=point_size,
        )

        plt.title(f"{group_col}: {group_value}")
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.gca().set_aspect("equal", "box")

        # Save figure
        filename = f"{group_col}_{group_value}.png"
        out_path=Path(outdir/filename)
        plt.savefig(out_path, dpi=300)
        plt.close()

    print(f"Saved plots to {outdir}")