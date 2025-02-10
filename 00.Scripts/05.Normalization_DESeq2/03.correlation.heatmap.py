import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
import random
import sys

def create_heatmap_and_dotplots(filepath):
    try:
        df = pd.read_csv(filepath, sep='\t', index_col=0)
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return
    except pd.errors.ParserError:
        print(f"Error: Could not parse the file at {filepath}. Check the separator.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during file reading: {e}")
        return

    replicate_groups = {}
    for col in df.columns:
        base_name = col.rsplit(".", 1)[0]
        if base_name not in replicate_groups:
            replicate_groups[base_name] = []
        replicate_groups[base_name].append(col)

    df_means = pd.DataFrame(index=df.index)
    for base_name, replicates in replicate_groups.items():
        df_means[base_name] = df[replicates].mean(axis=1)

    df_log2_means = np.log2(df_means + 1)

    # 1. Heatmap with clustering (clustering samples/columns)
    plt.figure(figsize=(10, 8))
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_log2_means.T) # Transpose for column clustering
    clustering = AgglomerativeClustering(linkage='ward')
    clustering.fit(scaled_data)
    sample_order = np.argsort(clustering.labels_)
    ordered_data = scaled_data[sample_order].T # Transpose back for heatmap

    sns.heatmap(ordered_data, cmap="viridis", yticklabels=df_log2_means.index, xticklabels=df_log2_means.columns[sample_order], cbar_kws={'label': 'Z-score (Log2 Mean Read Count)'})
    plt.title("Heatmap of Log2 Mean Read Counts with Hierarchical Clustering (Samples)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig("heatmap_means_samples.png")
    plt.show()

    # 2. Correlation dot plots (with random pairs)
    df_log2 = np.log2(df + 1)
    all_samples = df_log2.columns.tolist()

    for base_name, replicates in replicate_groups.items():
        n_replicates = len(replicates)
        if n_replicates >= 2:
            n_plots = n_replicates * (n_replicates - 1) // 2 + 5 # Add 5 random pairs
            fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 5))
            if n_plots == 1:
                axes = [axes]
            k = 0
            # Replicate pairs
            for i in range(n_replicates):
                for j in range(i + 1, n_replicates):
                    x = df_log2[replicates[i]]
                    y = df_log2[replicates[j]]
                    correlation, p_value = pearsonr(x, y)
                    sns.scatterplot(x=x, y=y, ax=axes[k])
                    axes[k].text(0.05, 0.95, f"r: {correlation:.2f}", transform=axes[k].transAxes)
                    axes[k].set_xlabel(replicates[i])
                    axes[k].set_ylabel(replicates[j])
                    axes[k].set_aspect('equal', adjustable='box')
                    k += 1

            # Random pairs
            random_pairs = []
            while len(random_pairs) < 5:
                sample1 = random.choice(all_samples)
                sample2 = random.choice(all_samples)
                if sample1 != sample2 and (sample1, sample2) not in random_pairs and (sample2, sample1) not in random_pairs:
                    random_pairs.append((sample1, sample2))
            
            for sample1, sample2 in random_pairs:
                x = df_log2[sample1]
                y = df_log2[sample2]
                correlation, p_value = pearsonr(x, y)
                sns.scatterplot(x=x, y=y, ax=axes[k])
                axes[k].text(0.05, 0.95, f"r: {correlation:.2f}", transform=axes[k].transAxes)
                axes[k].set_xlabel(sample1)
                axes[k].set_ylabel(sample2)
                axes[k].set_aspect('equal', adjustable='box')
                k += 1

            plt.suptitle(f"Correlation Dot Plots for {base_name} and Random Pairs")
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            plt.savefig(f"{base_name}_correlation_random.png")
            plt.show()
        else:
            print(f"Less than 2 replicates detected for {base_name}. Skipping dot plot generation.")


# Example usage:
filepath = sys.argv[1]
create_heatmap_and_dotplots(filepath)