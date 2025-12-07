# -*- coding: utf-8 -*-
"""
Created on Thu May 29 09:13:55 2025

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

import pandas as pd
import numpy as np





import matplotlib.pyplot as plt
import seaborn as sns

# Set up the plotting grid
features = [col for col in df.columns if col not in ['BoilingPoint', 'MeltingPoint', 'CID', 'SMILES']]
n_cols = 3
n_rows = int(np.ceil(len(features) / n_cols))

plt.figure(figsize=(18, n_rows * 4))

for i, feature in enumerate(features):
    plt.subplot(n_rows, n_cols, i + 1)
    sns.scatterplot(x=df[feature], y=df['BoilingPoint'], alpha=0.6, label='BoilingPoint')
    sns.scatterplot(x=df[feature], y=df['MeltingPoint'], alpha=0.6, label='MeltingPoint')
    plt.xlabel(feature)
    plt.ylabel('Temperature')
    plt.legend()
    plt.tight_layout()

plt.suptitle('Feature Influence on Boiling and Melting Points', y=1.02, fontsize=16)
plt.show()

# Compute correlation matrix
corr = df.drop(columns=['CID', 'SMILES']).corr()

# Extract correlations with BoilingPoint and MeltingPoint
target_corr = corr[['BoilingPoint', 'MeltingPoint']].drop(['BoilingPoint', 'MeltingPoint'])

# Plot heatmap
plt.figure(figsize=(8, len(target_corr) * 0.5))
sns.heatmap(target_corr.sort_values('BoilingPoint', ascending=False), annot=True, cmap='coolwarm')
plt.title('Feature Correlation with BoilingPoint and MeltingPoint')
plt.show()


#######################################################################################

# Drop non-feature columns
feature_data = df.drop(columns=['CID', 'SMILES', 'BoilingPoint', 'MeltingPoint'])

# Compute correlation matrix
feature_corr = feature_data.corr().abs()  # absolute value to catch both + and - correlations

# Create a mask to ignore the diagonal and duplicates (upper triangle)
mask = np.triu(np.ones_like(feature_corr, dtype=bool))

# Use mask to filter only the upper triangle
upper = feature_corr.where(mask)

# Find feature pairs with correlation > 0.9 (can adjust threshold)
high_corr_pairs = [
    (column, row, upper.loc[row, column])
    for column in upper.columns
    for row in upper.index
    if pd.notnull(upper.loc[row, column]) and upper.loc[row, column] > 0.9
]

# Display high correlation pairs
for feat1, feat2, corr_val in high_corr_pairs:
    print(f"{feat1} and {feat2} have correlation of {corr_val:.2f}")
    
# Compute correlation matrix
feature_corr = feature_data.corr()

# Plot the heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(feature_corr, 
            annot=True, 
            fmt=".2f", 
            cmap='coolwarm', 
            square=True, 
            cbar_kws={"shrink": .8},
            linewidths=0.5)

plt.title("Feature-to-Feature Correlation Heatmap (Collinearity)", fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()

########################################################################################3
    
    

# Shuffle the DataFrame
df_shuffled = df.sample(frac=1, random_state=42).reset_index(drop=True)

# Calculate the split index
split_index = int(0.8 * len(df_shuffled))


# Split into training and testing sets
train_df = df_shuffled.iloc[:split_index]
test_df = df_shuffled.iloc[split_index:]


#no need for dropping with PCA
#drop pyridine for melting point 
#can drop all chis, kappa 1, labuteasa for high mw correlation


##############################################################
from sklearn.preprocessing import StandardScaler

# Select feature columns (exclude target and ID columns)
features = train_df.drop(columns=['CID', 'SMILES', 'BoilingPoint', 'MeltingPoint']).dropna()
train_df = train_df.loc[features.index]

# Initialize the scaler
scaler = StandardScaler()

# Fit and transform the features
scaled_features = scaler.fit_transform(features)

# Convert back to DataFrame for readability
scaled_df = pd.DataFrame(scaled_features, columns=features.columns)

# Optionally, re-add target columns and IDs
scaled_df['BoilingPoint'] = train_df['BoilingPoint'].values
scaled_df['MeltingPoint'] = train_df['MeltingPoint'].values
scaled_df['CID'] = train_df['CID'].values
scaled_df['SMILES'] = train_df['SMILES'].values

# Reorder columns if desired
scaled_df = scaled_df[['CID', 'SMILES'] + list(features.columns) + ['BoilingPoint', 'MeltingPoint']]

##########################################################
from sklearn.decomposition import PCA
feature_cols = features.columns  # same as before
X_scaled = scaled_df[feature_cols].values

# Initialize PCA - can specify number of components or keep all
pca = PCA()

# Fit PCA on scaled features
pca.fit(X_scaled)

# Explained variance ratio
explained_variance = pca.explained_variance_ratio_

# Plot cumulative explained variance to decide how many components to keep
plt.figure(figsize=(8,5))
plt.plot(np.cumsum(explained_variance), marker='o')
plt.xlabel('Number of PCA Components')
plt.ylabel('Cumulative Explained Variance')
plt.grid(True)
plt.show()

# Choose number of components that explain e.g. 95% variance
n_components = np.argmax(np.cumsum(explained_variance) >= 0.95) + 1

# Apply PCA with chosen number of components
pca = PCA(n_components=8)
X_pca = pca.fit_transform(X_scaled)

# Convert to DataFrame for further processing
pca_columns = [f'PC{i+1}' for i in range(8)]
pca_df = pd.DataFrame(X_pca, columns=pca_columns)

# Reattach target columns and IDs
pca_df = pd.DataFrame(X_pca, columns=pca_columns)

# Add target and ID columns from train_df, which matches X_pca rows
pca_df['CID'] = train_df['CID'].values
pca_df['SMILES'] = train_df['SMILES'].values
pca_df['BoilingPoint'] = train_df['BoilingPoint'].values
pca_df['MeltingPoint'] = train_df['MeltingPoint'].values

########################################################
from sklearn.cluster import DBSCAN

# Select PCA components for clustering
X = pca_df[[f'PC{i+1}' for i in range(8)]].values

# Initialize and fit DBSCAN (adjust eps and min_samples as needed)
dbscan = DBSCAN(eps=0.9, min_samples=5)
clusters = dbscan.fit_predict(X)

# Add cluster labels to the DataFrame
pca_df['Cluster'] = clusters

# Optional: print cluster counts
print(pca_df['Cluster'].value_counts())

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', legend='full')
plt.title('DBSCAN Clusters on PCA Components')
plt.show()

cluster_summary = pca_df.groupby('Cluster')[['BoilingPoint', 'MeltingPoint']].mean()
print(cluster_summary)



loadings = pd.DataFrame(pca.components_.T, 
                        columns=[f'PC{i+1}' for i in range(pca.n_components_)], 
                        index=feature_cols)

print(loadings)


print(loadings['PC1'].sort_values(ascending=False))
