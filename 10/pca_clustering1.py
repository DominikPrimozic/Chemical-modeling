# -*- coding: utf-8 -*-
"""
Created on Thu May 29 11:15:41 2025

@author: domin
"""

import pandas as pd

# Load the CSV file
df = pd.read_csv('combined_pubchem_descriptors.csv')  # Replace 'your_file.csv' with your file path

first_cols = ['CID', 'SMILES']
other_cols = [col for col in df.columns if col not in first_cols]
df = df[first_cols + other_cols]

df = df[(df['BoilingPoint'] <= 1250) & (df['MeltingPoint'] <= 1250)]

columns_to_keep = [
    'CID', 'SMILES',
    'MeltingPoint', 'BoilingPoint',  
    'MolWt', 'MolLogP', 'TPSA',
    'fr_benzene', 'fr_pyridine', 'fr_furan',
    'LabuteASA', 'BalabanJ', 'BertzCT',
    'Chi0', 'Chi1', 'Chi2n', 'Chi3v', 'Chi4v',
    'Kappa1', 'Kappa2', 'Kappa3',
    'MaxPartialCharge', 'MinPartialCharge'
]

# Keep only the selected columns
df = df[columns_to_keep]

df = df[
    (df['TPSA'] <= 120) &
    (df['Kappa1'] <= 30) &
    (df['Kappa2'] >= -20) &
    (df['Kappa3'] <= 100) &
    (df['Chi4v'] <= 8) &
    (df['Chi3v'] <= 10) &
    (df['Chi1'] <= 20) &
    (df['Chi2n'] <= 18) &
    (df['Chi0'] <= 25) &
    (df['LabuteASA'] <= 250) &
    (df['BertzCT'] <= 1250)
]

import matplotlib.pyplot as plt
import seaborn as sns

df_boiling = df.dropna(subset=['BoilingPoint']).copy()
df_melting = df.dropna(subset=['MeltingPoint']).copy()

features = [col for col in df.columns if col not in ['BoilingPoint', 'MeltingPoint', 'CID', 'SMILES']]

# Plot for BoilingPoint
plt.figure(figsize=(18, len(features)//3 * 4))
for i, feature in enumerate(features):
    plt.subplot(len(features)//3 + 1, 3, i+1)
    sns.scatterplot(x=df_boiling[feature], y=df_boiling['BoilingPoint'], alpha=0.6)
    plt.xlabel(feature)
    plt.ylabel('BoilingPoint')
plt.suptitle('Feature Influence on Boiling Point')
plt.tight_layout()
plt.show()

# Plot for MeltingPoint
plt.figure(figsize=(18, len(features)//3 * 4))
for i, feature in enumerate(features):
    plt.subplot(len(features)//3 + 1, 3, i+1)
    sns.scatterplot(x=df_melting[feature], y=df_melting['MeltingPoint'], alpha=0.6)
    plt.xlabel(feature)
    plt.ylabel('MeltingPoint')
plt.suptitle('Feature Influence on Melting Point')
plt.tight_layout()
plt.show()


corr_boiling = df_boiling.drop(columns=['CID', 'SMILES']).corr()
corr_melting = df_melting.drop(columns=['CID', 'SMILES']).corr()

# Extract correlations with target
target_corr_boiling = corr_boiling['BoilingPoint'].drop('BoilingPoint')
target_corr_melting = corr_melting['MeltingPoint'].drop('MeltingPoint')

# Plot heatmaps separately
plt.figure(figsize=(6, len(target_corr_boiling)*0.4))
sns.heatmap(target_corr_boiling.to_frame(), annot=True, cmap='coolwarm', cbar=False)
plt.title('Feature Correlation with BoilingPoint')
plt.show()

plt.figure(figsize=(6, len(target_corr_melting)*0.4))
sns.heatmap(target_corr_melting.to_frame(), annot=True, cmap='coolwarm', cbar=False)
plt.title('Feature Correlation with MeltingPoint')
plt.show()


from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np

def scale_and_pca(df_target, target_name):
    features = df_target.drop(columns=['CID', 'SMILES', 'BoilingPoint', 'MeltingPoint']).copy()
    features = features.dropna()
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(features)
    
    pca = PCA()
    pca.fit(X_scaled)
    
    plt.figure(figsize=(8,5))
    plt.plot(np.cumsum(pca.explained_variance_ratio_), marker='o')
    plt.xlabel('Number of PCA Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.title(f'PCA Explained Variance for {target_name}')
    plt.grid(True)
    plt.show()
    
    n_components = np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.95) + 1
    
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X_scaled)
    
    pca_df = pd.DataFrame(X_pca, columns=[f'PC{i+1}' for i in range(n_components)], index=features.index)
    
    # Align by index labels, not positional
    pca_df[target_name] = df_target.loc[features.index, target_name].values
    pca_df['CID'] = df_target.loc[features.index, 'CID'].values
    pca_df['SMILES'] = df_target.loc[features.index, 'SMILES'].values
    
    return pca_df, pca, scaler

pca_boiling_df, pca_boiling, scaler_boiling = scale_and_pca(df_boiling, 'BoilingPoint')
pca_melting_df, pca_melting, scaler_melting = scale_and_pca(df_melting, 'MeltingPoint')


from sklearn.cluster import DBSCAN

def apply_dbscan(pca_df, n_components):
    X = pca_df[[f'PC{i+1}' for i in range(n_components)]].values
    dbscan = DBSCAN(eps=0.8, min_samples=10)
    clusters = dbscan.fit_predict(X)
    pca_df['Cluster'] = clusters
    print(pca_df['Cluster'].value_counts())
    return pca_df

pca_boiling_df = apply_dbscan(pca_boiling_df, pca_boiling.n_components_)
pca_melting_df = apply_dbscan(pca_melting_df, pca_melting.n_components_)

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_boiling_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=False)
plt.title('DBSCAN Clusters on PCA Components')
plt.show()

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_melting_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=False)
plt.title('DBSCAN Clusters on PCA Components')
plt.show()


plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
sns.boxplot(x='Cluster', y='BoilingPoint', data=pca_boiling_df)
plt.title('Melting Point distribution by Cluster')

plt.subplot(1,2,2)
sns.boxplot(x='Cluster', y='MeltingPoint', data=pca_melting_df)
plt.title('Boiling Point distribution by Cluster')

plt.tight_layout()
plt.show()

from sklearn.cluster import KMeans

def apply_kmeans(pca_df, n_components, n_clusters=4):
    # Extract the PCA component columns
    component_cols = [f'PC{i+1}' for i in range(n_components)]
    X = pca_df[component_cols].values

    # Fit KMeans
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X)

    # Add cluster labels to DataFrame
    pca_df = pca_df.copy()
    pca_df['Cluster'] = clusters
    
    return pca_df



pca_boiling_df = apply_kmeans(pca_boiling_df, pca_boiling.n_components_, n_clusters=4)
pca_melting_df = apply_kmeans(pca_melting_df, pca_melting.n_components_, n_clusters=4)


plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_boiling_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=False)
plt.title('KMeans Clusters on Boiling Point PCA')
plt.show()

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_melting_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=False)
plt.title('KMeans Clusters on Melting Point PCA')
plt.show()


from sklearn.cluster import AgglomerativeClustering

def apply_agglomerative(pca_df, n_components, n_clusters=4):
    # Extract the PCA component columns
    component_cols = [f'PC{i+1}' for i in range(n_components)]
    X = pca_df[component_cols].values

    # Fit Agglomerative Clustering
    agg = AgglomerativeClustering(n_clusters=n_clusters)
    clusters = agg.fit_predict(X)

    # Add cluster labels to DataFrame
    pca_df = pca_df.copy()
    pca_df['Cluster'] = clusters
    
    return pca_df


pca_boiling_df = apply_agglomerative(pca_boiling_df, pca_boiling.n_components_, n_clusters=4)
pca_melting_df = apply_agglomerative(pca_melting_df, pca_melting.n_components_, n_clusters=4)


plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_boiling_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=False)
plt.title('Agglomerative Clusters on Boiling Point PCA')
plt.show()

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_melting_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=False)
plt.title('Agglomerative Clusters on Melting Point PCA')
plt.show()


from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

def plot_dendrogram(pca_df, n_components, method='ward', title='Dendrogram'):
    from scipy.cluster.hierarchy import linkage, dendrogram

    component_cols = [f'PC{i+1}' for i in range(n_components)]
    X = pca_df[component_cols].values

    # Compute linkage matrix
    Z = linkage(X, method=method)

    # Plot dendrogram
    plt.figure(figsize=(12, 6))
    dendrogram(Z, truncate_mode='lastp', p=30, leaf_rotation=90, leaf_font_size=10, show_contracted=True)
    plt.title(title)
    plt.xlabel('Sample Index or (Cluster Size)')
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.show()




plot_dendrogram(pca_boiling_df, pca_boiling.n_components_, title='Boiling Point PCA - Dendrogram')


plot_dendrogram(pca_melting_df, pca_melting.n_components_, title='Melting Point PCA - Dendrogram')
