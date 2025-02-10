import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys, re

def plot_correlation_heatmap(file_path, outpath, groups):
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    df_log2 = np.log2(df + 0.01)

    for group_name, patterns in groups.items():
        # Filter columns based on regex patterns
        selected_cols = []
        for pattern in patterns:
            selected_cols.extend([col for col in df_log2.columns if re.search(pattern, col)])
        
        if not selected_cols: #check if any sample selected
            print(f"No samples found for group: {group_name}")
            continue

        group_df = df_log2[selected_cols]
        correlation_matrix = group_df.corr(method='pearson')
        correlation_matrix = correlation_matrix.fillna(0)
        
        if correlation_matrix.empty: #check if the correlation_matrix is empty
            print(f"No correlation can be calculated for group: {group_name}. May be only one sample in this group.")
            continue

        n_sample = len(selected_cols)
        plt.figure(figsize=(0.5*n_sample, 0.4*n_sample))
        sns.heatmap(correlation_matrix, cmap="coolwarm", vmin=0.8, vmax=1)
        plt.title(f"Correlation Heatmap: {group_name}")
        plt.tight_layout()
        plt.savefig(outpath  + "/corr.heatmap." + f"{group_name}".replace(" ","_") + ".png")
        plt.show()


if __name__ == "__main__":
    filepath = sys.argv[1]
    outpath = sys.argv[2]
    groups = {
        "All vs All": [".*"],
        "Head": [r".*h[123]$"], 
        "Body": [r".*b[123]$"]
    }
    try:
        plot_correlation_heatmap(filepath, outpath, groups)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
    except Exception as e: # Catch any other potential exceptions
        print(f"An unexpected error occurred: {e}")


