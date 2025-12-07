import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen, Lipinski
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial.distance import euclidean, mahalanobis

# Load and clean data
df = pd.read_csv('logS.csv')
df = df.iloc[2:].reset_index(drop=True)

for col in df.columns:
    if col != 'SMILES':
        df[col] = pd.to_numeric(df[col], errors='coerce')

df = df.dropna().reset_index(drop=True)

# Descriptor calculations
def calc_tpsa(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcTPSA(mol) if mol else None

df['TPSA'] = df['SMILES'].apply(calc_tpsa)
df = df.dropna(subset=['TPSA']).reset_index(drop=True)

# Calculate logP
df['logP'] = df['SMILES'].apply(lambda smi: Crippen.MolLogP(Chem.MolFromSmiles(smi)))

# Calculate H-bond donors, acceptors, and rotatable bonds
df['HBD'] = df['SMILES'].apply(lambda s: Lipinski.NumHDonors(Chem.MolFromSmiles(s)))
df['HBA'] = df['SMILES'].apply(lambda s: Lipinski.NumHAcceptors(Chem.MolFromSmiles(s)))
df['RotBonds'] = df['SMILES'].apply(lambda s: Lipinski.NumRotatableBonds(Chem.MolFromSmiles(s)))

# Define features
features = ['TPSA', 'logP', 'HBD', 'HBA', 'RotBonds']

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df[features])
mean_vector = X_scaled.mean(axis=0)

# Euclidean distance
df['euclidean_distance'] = np.linalg.norm(X_scaled - mean_vector, axis=1)

# Mahalanobis distance
def compute_mahalanobis_distances(X):
    mean_vec = X.mean(axis=0)
    cov_matrix = np.cov(X, rowvar=False)
    inv_cov_matrix = np.linalg.inv(cov_matrix)
    return np.array([mahalanobis(x, mean_vec, inv_cov_matrix) for x in X])

df['mahalanobis_distance'] = compute_mahalanobis_distances(X_scaled)

# Visualizations
def plot_distance_histogram(df):
    plt.figure(figsize=(8, 5))
    plt.hist(df['euclidean_distance'] , bins=30, color='skyblue', edgecolor='k')
    plt.xlabel('Distance from Mean')
    plt.ylabel('Number of Points')
    plt.title('Euclidean Distance Histogram')
    plt.show()

def plot_mahalanobis_histogram(df):
    plt.figure(figsize=(8, 5))
    plt.hist(df['mahalanobis_distance'], bins=30, color='orange', edgecolor='k')
    plt.xlabel('Mahalanobis Distance from Mean')
    plt.ylabel('Number of Points')
    plt.title('Mahalanobis Distance Histogram')
    plt.show()

def plot_distance_vs_logS(df):
    plt.figure(figsize=(10, 5))
    plt.scatter(df['logS'], df['euclidean_distance'], color='green', edgecolor='k', alpha=0.7)
    plt.xlabel('logS Point')
    plt.ylabel('Euclidean Distance')
    plt.title('Euclidean Distance vs. logS Point')
    plt.show()

def plot_mahalanobis_vs_logS(df):
    plt.figure(figsize=(10, 5))
    plt.scatter(df['logS'], df['mahalanobis_distance'], color='purple', edgecolor='k', alpha=0.7)
    plt.xlabel('logS ')
    plt.ylabel('Mahalanobis Distance')
    plt.title('Mahalanobis Distance vs. logS Point')
    plt.show()

# Plot distances
plot_distance_histogram(df)
plot_distance_vs_logS(df)
plot_mahalanobis_vs_logS(df)
plot_mahalanobis_histogram(df)

# Identify outliers (top 5%)
top_percent = 0.05
top_n = int(len(df) * top_percent)

top_mahalanobis = df.sort_values(by='mahalanobis_distance', ascending=False).head(top_n)
top_euclidean = df.sort_values(by='euclidean_distance', ascending=False).head(top_n)

# Output descriptors of top outliers
cols = ['SMILES', 'logS', 'TPSA', 'logP', 'HBD', 'HBA', 'RotBonds',
        'mahalanobis_distance', 'euclidean_distance']

print("\nTop Outliers by Mahalanobis Distance:")
print(top_mahalanobis[cols].sort_values(by='mahalanobis_distance', ascending=False))

print("\nTop Outliers by Euclidean Distance:")
print(top_euclidean[cols].sort_values(by='euclidean_distance', ascending=False))

# Label outliers
df['is_mahalanobis_outlier'] = df.index.isin(top_mahalanobis.index)
df['is_euclidean_outlier'] = df.index.isin(top_euclidean.index)

# Function to plot boxplots
def plot_box_comparison(df, feature, outlier_col, title_suffix):
    plt.figure(figsize=(6, 5))
    sns.boxplot(x=df[outlier_col], y=df[feature], palette='Set2')
    plt.xlabel('Outlier' if 'outlier' in outlier_col else outlier_col)
    plt.ylabel(feature)
    plt.title(f'{feature} vs. Outlier Status ({title_suffix})')
    plt.xticks([0, 1], ['No', 'Yes'])
    plt.tight_layout()
    plt.show()

# Plot for each descriptor
descriptor_cols = ['logP', 'TPSA', 'HBD', 'HBA', 'RotBonds']

for feature in descriptor_cols:
    plot_box_comparison(df, feature, 'is_mahalanobis_outlier', 'Mahalanobis')
    plot_box_comparison(df, feature, 'is_euclidean_outlier', 'Euclidean')

# Summary statistics
def summarize_feature_differences(df, feature_cols, outlier_col):
    print(f"\nSummary Statistics for {outlier_col}:")
    for col in feature_cols:
        group_stats = df.groupby(outlier_col)[col].agg(['mean', 'median'])
        print(f"\n{col}:")
        print(group_stats)

summarize_feature_differences(df, descriptor_cols, 'is_mahalanobis_outlier')
summarize_feature_differences(df, descriptor_cols, 'is_euclidean_outlier')

