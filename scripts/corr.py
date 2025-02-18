import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys, re

def main():
    filepath = snakemake.input
    outfile = snakemake.output
    target_groups = snakemake.params.target_groups
    sample_mapping = snakemake.params.sample_mapping
    valid_samples = [s for s in counts_df.columns if sample_mapping.get(s.split(".")[0]) in target_groups]

    df = pd.read_csv(filepath, sep='\t', index_col=0)
    df_log2 = np.log2(df + 0.01)

    group_df = df_log2[valid_samples]
    correlation_matrix = group_df.corr(method='pearson')
    correlation_matrix = correlation_matrix.fillna(0)

    if correlation_matrix.empty: #check if the correlation_matrix is empty
        print(f"No correlation can be calculated for group: {group_name}. May be only one sample in this group.")
        continue

    n_sample = len(samples)
    plt.figure(figsize=(0.4*n_sample+1, 0.4*n_sample))
    sns.heatmap(correlation_matrix, cmap="coolwarm", vmin=0.8, vmax=1)
    #plt.title(f"Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(outfile)
    plt.show()

if __name__ == "__main__":
    main()
