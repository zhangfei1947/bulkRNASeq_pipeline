import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import sys

def main():
    infile = snakemake.input

    df = pd.read_csv(infile, sep='\t', index_col=0)
    df_log2 = np.log2(df + 0.01)

    # Scaling the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_log2.T) # Transpose for PCA on samples

    # Performing PCA
    pca = PCA(n_components=3)  # Reduce to 3 principal components
    pca_result = pca.fit_transform(scaled_data)

    # Creating a DataFrame for PCA results
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2', 'PC3'])
    pca_df['Sample'] = df.columns

    color_mapping = snakemake.params.color_mapping
    sample_colors = snakemake.params.sample_colors

    # Extracting metadata for coloring
    pca_df['Tissue'] = pca_df['Sample'].str.split("_").str[2].str[0]
    pca_df['Day'] = pca_df['Sample'].str.split("_").str[0]
    pca_df['Treatment'] = pca_df['Sample'].str.split("_").str[1]

    # Plotting
    n = 5
    fig, axes = plt.subplots(1, n, figsize=(n*(4.2), 4))

    sns.scatterplot(x='PC1', y='PC2', hue='Tissue', data=pca_df, ax=axes[0], s=100)
    axes[0].set_title('PCA Colored by Tissue')

    sns.scatterplot(x='PC1', y='PC2', hue='Day', data=pca_df, ax=axes[1], s=100)
    axes[1].set_title('PCA Colored by Day')

    sns.scatterplot(x='PC1', y='PC2', hue='Treatment', data=pca_df, ax=axes[2], s=100)
    axes[2].set_title('PCA Colored by Treatment')

    sns.scatterplot(x='PC2', y='PC3', hue='Treatment', data=pca_df, ax=axes[3], s=100)
    axes[3].set_title('PCA Colored by Treatment')

    sns.scatterplot(x='PC1', y='PC3', hue='Treatment', data=pca_df, ax=axes[4], s=100)
    axes[4].set_title('PCA Colored by Treatment')

    for ax in axes[:3]:
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)') #add explained variance ratio
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)') #add explained variance ratio
        ax.set_aspect('auto', adjustable='box') # Make subplots square

'''
    axes[3].set_xlabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)') #add explained variance ratio
    axes[3].set_ylabel(f'PC3 ({pca.explained_variance_ratio_[2]*100:.2f}%)') #add explained variance ratio
    axes[3].set_aspect('auto', adjustable='box') # Make subplots square

    axes[4].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)') #add explained variance ratio
    axes[4].set_ylabel(f'PC3 ({pca.explained_variance_ratio_[2]*100:.2f}%)') #add explained variance ratio
    axes[4].set_aspect('auto', adjustable='box') # Make subplots square
'''

    plt.tight_layout()
    plt.savefig(outpath+"/pca_plots.pdf")
    plt.savefig(outpath+"/pca_plots.png")
    plt.show()

if __name__ == "__main__":
    main()

