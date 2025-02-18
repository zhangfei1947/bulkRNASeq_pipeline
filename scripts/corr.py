import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys, re

def main():
    filepath = snakemake.input
    print(filepath)
    outfile = snakemake.output
    target_groups = snakemake.params.target_groups
    sample_mapping = snakemake.params.sample_mapping

    df = pd.read_csv(filepath, sep='\t', index_col=0)
    df_log2 = np.log2(df + 0.01)
    valid_samples = [s for s in df_log2.columns if sample_mapping.get(s.split(".")[0]) in target_groups]

    group_df = df_log2[valid_samples]
    correlation_matrix = group_df.corr(method='pearson')
    correlation_matrix = correlation_matrix.fillna(0)

    n_sample = len(samples)
    plt.figure(figsize=(0.4*n_sample+1, 0.4*n_sample))
    sns.heatmap(correlation_matrix, cmap="coolwarm", vmin=0.8, vmax=1)
    #plt.title(f"Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(outfile)
    plt.show()

if __name__ == "__main__":
    main()
