import pandas as pd
import matplotlib.pyplot as plt
import io

all_cells_path = '/scratch/mfafouti/Mommybrain/IMN_Trajectory/out/subclass_replicate_representation_all_cells.csv'
data = pd.read_csv(all_cells_path)

# subset_path = ''
data.set_index("subclass_name", inplace=True)

# 3. Create the plot
# We only plot CORT and OIL columns
ax = data[["CORT", "OIL"]].plot(kind='bar', figsize=(30, 6), width=0.8)

# 4. Formatting
plt.title("Comparison of CORT and OIL levels by Subclass", fontsize=14)
plt.xlabel("Subclass Name", fontsize=12)
plt.ylabel("Value", fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend(title="Condition")

plt.tight_layout()

plt.savefig('/scratch/mfafouti/Mommybrain/IMN_Trajectory/out/counts/all_cells_counts.png')