# -*- coding: utf-8 -*-
"""
Created on Thu May 29 12:18:44 2025

@author: domin
"""

##combining 
# my + boiling
#my + melting

#water solubility clustering too
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

df = pd.read_csv('boiling.csv')
df = df.iloc[2:].reset_index(drop=True)

for col in df.columns:
    if col != 'SMILES':
        df[col] = pd.to_numeric(df[col], errors='coerce')

# Drop all rows with any NaNs
df = df.dropna().reset_index(drop=True)


##add TPS with RDKIt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Define a function to calculate TPSA from a SMILES string
def calc_tpsa(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return rdMolDescriptors.CalcTPSA(mol)
    else:
        return None
# Define functions to calculate H-bond donors and acceptors
def calc_hbd(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return rdMolDescriptors.CalcNumHBD(mol)
    else:
        return None

def calc_hba(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return rdMolDescriptors.CalcNumHBA(mol)
    else:
        return None

# Apply these functions to the dataframe
#df['HBD'] = df['SMILES'].apply(calc_hbd)
#df['HBA'] = df['SMILES'].apply(calc_hba)

# Apply the function to the 'SMILES' column and create a new 'TPSA' column
df['TPSA'] = df['SMILES'].apply(calc_tpsa)

# Optionally, drop rows where TPSA calculation failed (None values)
df = df.dropna(subset=['TPSA']).reset_index(drop=True)
##########################################################################################

# Define feature columns
features = [col for col in df.columns if col not in ['boiling_point', 'SMILES']]


# Plot for BoilingPoint
plt.figure(figsize=(18, len(features)//3 * 4))
for i, feature in enumerate(features):
    plt.subplot(len(features)//3 + 1, 3, i+1)
    sns.scatterplot(x=df[feature], y=df['boiling_point'], alpha=0.6)
    plt.xlabel(feature)
    plt.ylabel('Boiling Point')

    # Limit the number of ticks to 4 on both axes
    x_vals = df[feature].dropna()
    y_vals = df['boiling_point'].dropna()
    
    if not x_vals.empty:
        x_ticks = np.linspace(x_vals.min(), x_vals.max(), 4)
        plt.xticks(x_ticks.round(2))
    if not y_vals.empty:
        y_ticks = np.linspace(y_vals.min(), y_vals.max(), 4)
        plt.yticks(y_ticks.round(2))

plt.suptitle('Feature Influence on Boiling Point')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np

def scale_and_pca(df_target, target_name):
    features = df_target.drop(columns=['SMILES']).copy()
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
    pca_df['SMILES'] = df_target.loc[features.index, 'SMILES'].values
    
    return pca_df, pca, scaler

pca_boiling_df, pca_boiling, scaler_boiling = scale_and_pca(df, 'boiling_point')



from sklearn.cluster import DBSCAN

def apply_dbscan(pca_df, n_components):
    X = pca_df[[f'PC{i+1}' for i in range(n_components)]].values
    dbscan = DBSCAN(eps=1, min_samples=10)
    clusters = dbscan.fit_predict(X)
    pca_df['Cluster'] = clusters
    print(pca_df['Cluster'].value_counts())
    return pca_df

pca_boiling_df = apply_dbscan(pca_boiling_df, pca_boiling.n_components_)


import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_boiling_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=True)
plt.title('DBSCAN Clusters on PCA Components')
plt.show()

df = df.merge(pca_boiling_df[['SMILES', 'Cluster']], on='SMILES', how='left')


import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

features_to_check = ['MW', 'SLogP', 'boiling_point', "TPSA", "nRot", "SMR", 'nHBAcc', 'nHBDon', 'Zagreb2']
n_features = len(features_to_check)

# Calculate layout: 3 columns is a good visual balance
n_cols = 3
n_rows = (n_features + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*6, n_rows*5))
axes = axes.flatten()

anova_results = []

for i, feature in enumerate(features_to_check):
    ax = axes[i]
    sns.boxplot(x='Cluster', y=feature, data=df, ax=ax)
    ax.set_title(f'{feature} by Cluster', fontsize=12)
    ax.set_xlabel('Cluster')
    ax.set_ylabel(feature)

    # Perform ANOVA test (exclude cluster -1 if present)
    grouped = [df[df['Cluster'] == c][feature].dropna() for c in df['Cluster'].unique() if c != -1]
    f_stat, p_val = stats.f_oneway(*grouped)
    anova_results.append((feature, f_stat, p_val))

    # Annotate p-value on the plot
    ax.text(0.5, 0.95, f"p = {p_val:.4f}", ha='center', va='top', transform=ax.transAxes, fontsize=10)

# Hide unused subplots if any
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.suptitle('Feature Distributions by Cluster with ANOVA p-values', fontsize=16, y=1.02)
plt.show()

# Print ANOVA results below the plots
for feature, f_stat, p_val in anova_results:
    print(f"ANOVA for {feature}: F = {f_stat:.2f}, p = {p_val:.4f}")
    
from sklearn.cluster import KMeans
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def plot_elbow(pca_df, n_components, max_clusters=10):
    component_cols = [f'PC{i+1}' for i in range(n_components)]
    X = pca_df[component_cols].values

    inertias = []
    cluster_range = range(1, max_clusters + 1)

    for k in cluster_range:
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(X)
        inertias.append(kmeans.inertia_)

    plt.figure(figsize=(8, 5))
    plt.plot(cluster_range, inertias, marker='o')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Inertia (Within-Cluster SSE)')
    plt.title('Elbow Method for Optimal K')
    plt.grid(True)
    plt.show()
    
def apply_kmeans(pca_df, n_components, n_clusters=7):
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

plot_elbow(pca_boiling_df, n_components=pca_boiling.n_components_, max_clusters=10)


pca_boiling_df = apply_kmeans(pca_boiling_df, pca_boiling.n_components_, n_clusters=7)
df = df.drop(columns='Cluster', errors='ignore')
df = df.merge(pca_boiling_df[['SMILES', 'Cluster']], on='SMILES', how='left')

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_boiling_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=True)
plt.title('KMeans Clusters on Boiling Point PCA')
plt.show()
features_to_check = ['MW', 'SLogP', 'boiling_point', "TPSA", "nRot", "SMR", 'nHBAcc', 'nHBDon', 'Zagreb2']
n_features = len(features_to_check)

# Calculate layout: 3 columns is a good visual balance
n_cols = 3
n_rows = (n_features + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*6, n_rows*5))
axes = axes.flatten()

anova_results = []

for i, feature in enumerate(features_to_check):
    ax = axes[i]
    sns.boxplot(x='Cluster', y=feature, data=df, ax=ax)
    grouped = [df[df['Cluster'] == c][feature].dropna() for c in df['Cluster'].unique() if c != -1]
    

    ax.set_title(f'{feature} by Cluster', fontsize=12)
    ax.set_xlabel('Cluster')
    ax.set_ylabel(feature)

    # Perform ANOVA test (exclude cluster -1 if present)

    f_stat, p_val = stats.f_oneway(*grouped)
    anova_results.append((feature, f_stat, p_val))

    # Annotate p-value on the plot
    ax.text(0.5, 0.95, f"p = {p_val:.4f}", ha='center', va='top', transform=ax.transAxes, fontsize=10)

# Hide unused subplots if any
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.suptitle('Feature Distributions by Cluster with ANOVA p-values', fontsize=16, y=1.02)
plt.show()


from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage

def plot_dendrogram(pca_df, n_components):
    component_cols = [f'PC{i+1}' for i in range(n_components)]
    X = pca_df[component_cols].values

    linked = linkage(X, method='ward')

    plt.figure(figsize=(10, 6))
    dendrogram(linked, truncate_mode='level', p=5)
    plt.title('Hierarchical Clustering Dendrogram (Truncated)')
    plt.xlabel('Data Points')
    plt.ylabel('Distance')
    plt.show()
    
def apply_agglomerative(pca_df, n_components, n_clusters=7):
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

plot_dendrogram(pca_boiling_df, n_components=pca_boiling.n_components_)

pca_boiling_df = apply_agglomerative(pca_boiling_df, pca_boiling.n_components_, n_clusters=7)
df = df.drop(columns='Cluster', errors='ignore')
df = df.merge(pca_boiling_df[['SMILES', 'Cluster']], on='SMILES', how='left')

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_boiling_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend=True)
plt.title('Agglomerative Clusters on Boiling Point PCA')
plt.show()



features_to_check = ['MW', 'SLogP', 'boiling_point', "TPSA", "nRot", "SMR", 'nHBAcc', 'nHBDon', 'Zagreb2']
n_features = len(features_to_check)

# Calculate layout: 3 columns is a good visual balance
n_cols = 3
n_rows = (n_features + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*6, n_rows*5))
axes = axes.flatten()

anova_results = []

for i, feature in enumerate(features_to_check):
    ax = axes[i]
    sns.boxplot(x='Cluster', y=feature, data=df, ax=ax)
    grouped = [df[df['Cluster'] == c][feature].dropna() for c in df['Cluster'].unique() if c != -1]


    ax.set_title(f'{feature} by Cluster', fontsize=12)
    ax.set_xlabel('Cluster')
    ax.set_ylabel(feature)

    # Perform ANOVA test (exclude cluster -1 if present)
    

    f_stat, p_val = stats.f_oneway(*grouped)
    anova_results.append((feature, f_stat, p_val))

    # Annotate p-value on the plot
    ax.text(0.5, 0.95, f"p = {p_val:.4f}", ha='center', va='top', transform=ax.transAxes, fontsize=10)

# Hide unused subplots if any
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.suptitle('Feature Distributions by Cluster with ANOVA p-values', fontsize=16, y=1.02)
plt.show()


#############################################################################################
##Euclidian and Mahalabis distance
###########################################################################################

from scipy.spatial.distance import euclidean

def find_outliers_by_distance(df, feature_cols, id_col='SMILES', label_col='boiling_point', top_n=10, plot=False):
    # Step 1: Standardize the features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df[feature_cols])

    # Step 2: Compute mean vector
    mean_vector = X_scaled.mean(axis=0)

    # Step 3: Compute Euclidean distance from mean
    distances = np.linalg.norm(X_scaled - mean_vector, axis=1)
    df['distance_from_mean'] = distances

    # Step 4: Sort by distance
    df_sorted = df.sort_values(by='distance_from_mean', ascending=False).reset_index(drop=True)

    # Step 5: Optional visualization using PCA
    if plot:
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_scaled)
        
        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], 
                              c=distances, cmap='viridis', s=30, edgecolor='k', alpha=0.7)
        plt.colorbar(scatter, label='Distance from Mean')
        plt.title('PCA Projection Colored by Distance from Mean')
        plt.xlabel('PCA 1')
        plt.ylabel('PCA 2')

        # Highlight top N outliers
        for i in range(top_n):
            idx = df_sorted.index[i]
            plt.annotate(df.loc[idx, id_col], 
                         (X_pca[idx, 0], X_pca[idx, 1]), 
                         fontsize=8, color='red')

        plt.tight_layout()
        plt.show()

    return df_sorted[['distance_from_mean', id_col, label_col] + feature_cols]


outlier_df = find_outliers_by_distance(df, features, id_col='SMILES', label_col='boiling_point', top_n=5)
print(outlier_df.head())

import matplotlib.pyplot as plt

def plot_distance_histogram(df):
    plt.figure(figsize=(8, 5))
    plt.hist(df['distance_from_mean'], bins=30, color='skyblue', edgecolor='k')
    plt.xlabel('Distance from Mean')
    plt.ylabel('Number of Points')
    plt.title('Distribution of Distances from Mean')
    plt.show()
    
def plot_distance_vs_boiling(df):
    plt.figure(figsize=(10, 5))
    plt.scatter(df['boiling_point'], df['distance_from_mean'], color='mediumseagreen', edgecolor='k', alpha=0.7)
    plt.xlabel('Boiling Point')
    plt.ylabel('Distance from Mean')
    plt.title('Distance from Mean vs. Boiling Point')
    plt.show()
    
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df[features])
mean_vector = X_scaled.mean(axis=0)
df['distance_from_mean'] = np.linalg.norm(X_scaled - mean_vector, axis=1)

# Then plot
plot_distance_histogram(df)
plot_distance_vs_boiling(df)





from scipy.spatial.distance import mahalanobis

def compute_mahalanobis_distances(X):
    mean_vec = X.mean(axis=0)
    cov_matrix = np.cov(X, rowvar=False)
    inv_cov_matrix = np.linalg.inv(cov_matrix)
    
    # Compute Mahalanobis distance for each row
    m_distances = np.array([mahalanobis(x, mean_vec, inv_cov_matrix) for x in X])
    return m_distances

scaler = StandardScaler()
X_scaled = scaler.fit_transform(df[features])

# Compute Euclidean distances (you already do this)
mean_vector = X_scaled.mean(axis=0)
df['euclidean_distance'] = np.linalg.norm(X_scaled - mean_vector, axis=1)

# Compute Mahalanobis distances
df['mahalanobis_distance'] = compute_mahalanobis_distances(X_scaled)

def plot_mahalanobis_histogram(df):
    plt.figure(figsize=(8, 5))
    plt.hist(df['mahalanobis_distance'], bins=30, color='orange', edgecolor='k')
    plt.xlabel('Mahalanobis Distance from Mean')
    plt.ylabel('Number of Points')
    plt.title('Distribution of Mahalanobis Distances from Mean')
    plt.show()


plot_distance_histogram(df)          # Euclidean histogram
plot_mahalanobis_histogram(df)       # Mahalanobis histogram

plt.figure(figsize=(10, 5))
plt.scatter(df['boiling_point'], df['mahalanobis_distance'], 
            color='purple', edgecolor='k', alpha=0.7)
plt.xlabel('Boiling Point')
plt.ylabel('Mahalanobis Distance from Mean')
plt.title('Mahalanobis Distance vs. Boiling Point')
plt.show()

threshold_percentile = 0.95  # top 5%
threshold_euc = np.percentile(df['euclidean_distance'], threshold_percentile * 100)
threshold_mah = np.percentile(df['mahalanobis_distance'], threshold_percentile * 100)

euc_outliers = df[df['euclidean_distance'] >= threshold_euc]
mah_outliers = df[df['mahalanobis_distance'] >= threshold_mah]


print(euc_outliers[['SMILES', 'boiling_point', 'TPSA'] + features])
print(mah_outliers[['SMILES', 'boiling_point', 'TPSA'] + features])


from rdkit.Chem import Crippen

# Add logP
df['logP'] = df['SMILES'].apply(lambda smi: Crippen.MolLogP(Chem.MolFromSmiles(smi)))
top_mahalanobis = df.sort_values(by='mahalanobis_distance', ascending=False).head(int(0.05 * len(df)))
top_euclidean = df.sort_values(by='euclidean_distance', ascending=False).head(int(0.05 * len(df)))

print(top_mahalanobis[['SMILES', 'boiling_point', 'TPSA', 'logP']])
from rdkit.Chem import Lipinski

df['HBD'] = df['SMILES'].apply(lambda s: Lipinski.NumHDonors(Chem.MolFromSmiles(s)))
df['HBA'] = df['SMILES'].apply(lambda s: Lipinski.NumHAcceptors(Chem.MolFromSmiles(s)))
df['RotBonds'] = df['SMILES'].apply(lambda s: Lipinski.NumRotatableBonds(Chem.MolFromSmiles(s)))
cols = ['SMILES', 'boiling_point', 'TPSA', 'logP', 'HBD', 'HBA', 'RotBonds', 'mahalanobis_distance']
print(top_mahalanobis[cols].sort_values(by='mahalanobis_distance', ascending=False))

