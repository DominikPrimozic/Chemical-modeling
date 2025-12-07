import pandas as pd
import pubchempy as pcp
from tqdm import tqdm

# Load descriptors and take first 100 for testing
rdkit_df = pd.read_csv('pubchem_descriptors.csv')
rdkit_df = rdkit_df.head(50).copy()
rdkit_df['CID'] = rdkit_df['CID'].astype(int)
cids = rdkit_df['CID'].tolist()

# Function to get experimental values
def fetch_experimental_props(cid):
    try:
        c = pcp.Compound.from_cid(cid)
        d = c.to_dict()
        return {
            'CID': cid,
            'MeltingPoint': d.get('melting_point'),
            'BoilingPoint': d.get('boiling_point'),
            'WaterSolubility': d.get('solubility')
        }
    except Exception as e:
        print(f"Failed for CID {cid}: {e}")
        return {'CID': cid, 'MeltingPoint': None, 'BoilingPoint': None, 'WaterSolubility': None}

# Fetch for 100 molecules
experimental_data = [fetch_experimental_props(cid) for cid in tqdm(cids)]

# Merge results
experimental_df = pd.DataFrame(experimental_data)
experimental_df['CID'] = experimental_df['CID'].astype(int)

# Merge and select relevant columns
full_df = pd.merge(rdkit_df, experimental_df, on='CID', how='left')

columns_to_keep = [
    'CID', 'SMILES',
    'MeltingPoint', 'BoilingPoint', 'WaterSolubility',
    'MolWt', 'MolLogP', 'TPSA',
    'fr_benzene', 'fr_pyridine', 'fr_furan',
    'LabuteASA', 'BalabanJ', 'BertzCT',
    'Chi0', 'Chi1', 'Chi2n', 'Chi3v', 'Chi4v',
    'Kappa1', 'Kappa2', 'Kappa3',
    'MaxPartialCharge', 'MinPartialCharge'
]

result_df = full_df[columns_to_keep]


result_df.to_csv('combined_dataset_test.csv', index=False)
print(result_df.head())
