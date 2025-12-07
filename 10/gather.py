# -*- coding: utf-8 -*-
"""
Created on Wed May 28 16:15:48 2025

@author: domin
"""

import gzip

"""
input_gz = 'CID-PMID.gz'
output_txt = 'CID-IUPAC.txt'

with gzip.open(input_gz, 'rt', encoding='utf-8') as f_in, open(output_txt, 'w', encoding='utf-8') as f_out:
    for line in f_in:
        f_out.write(line)
print("Decompression done!")

cid_list_file = 'cid_list2.txt'

with open(output_txt, 'r', encoding='utf-8') as fin, open(cid_list_file, 'w', encoding='utf-8') as fout:
    for line in fin:
        parts = line.strip().split('\t')
        if parts and parts[0].isdigit():
            fout.write(parts[0] + '\n')
print("CID list extracted!")
"""
import random
cid_list_file = 'cid_list2.txt'

with open(cid_list_file, 'r') as f:
    cids = [int(line.strip()) for line in f if line.strip().isdigit()]
    
cids=cids[20000:]
#print(f"Total CIDs loaded: {len(cids)}")

import gzip
import random
import requests
import pubchempy as pcp
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from tqdm import tqdm
import re
import time


# Set up descriptor calculator
descriptor_names = [desc[0] for desc in Descriptors._descList]
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

# Helper to parse boiling/melting from PubChem PUG View JSON
def extract_value(section):
    info_list = section.get('Information', [])
    for info in info_list:
        val = info.get('Value', {})
        markup = val.get('StringWithMarkup', [])
        if markup:
            text = ''.join([m['String'] for m in markup])
            match = re.search(r'-?\d+\.?\d*', text)
            if match:
                return float(match.group())
    return None

def fetch_boiling_melting(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/"
    try:
        res = requests.get(url, timeout=10)
        res.raise_for_status()
        data = res.json()
    except Exception:
        return None, None

    melting = boiling = None
    try:
        sections = data['Record']['Section']
        for section in sections:
            if section.get('TOCHeading') == 'Chemical and Physical Properties':
                for sub in section.get('Section', []):
                    if sub.get('TOCHeading') == 'Experimental Properties':
                        for prop in sub.get('Section', []):
                            heading = prop.get('TOCHeading', '').lower()
                            if 'melting point' in heading:
                                melting = extract_value(prop)
                            elif 'boiling point' in heading:
                                boiling = extract_value(prop)
    except Exception:
        return None, None

    return boiling, melting

# Sample and collect data
random.shuffle(cids)
results = []
max_needed = 2000  # Adjust how many full-property compounds you want
checked = 0
fail = 0

for cid in tqdm(cids, desc="Filtering compounds"):
    if len(results) >= max_needed:
        break

    try:
        compound = pcp.Compound.from_cid(cid)
        smiles = compound.isomeric_smiles
        xlogp = compound.xlogp
        if not smiles or xlogp is None:
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        boiling, melting = fetch_boiling_melting(cid)
        if boiling is None or melting is None:
            continue

        descriptors = calculator.CalcDescriptors(mol)
        results.append([boiling, melting, xlogp, cid, smiles] + list(descriptors))
    except Exception:
        fail += 1
    checked += 1

    time.sleep(0.1)  # Respect rate limit

print(f"✔ Found {len(results)} usable compounds after checking {checked}. Failed: {fail}")

# Save
columns = ['BoilingPoint', 'MeltingPoint', 'LogP', 'CID', 'SMILES'] + descriptor_names
df = pd.DataFrame(results, columns=columns)
df.to_csv("filtered_pubchem_descriptors22.csv", index=False)
print("✅ CSV saved: filtered_pubchem_descriptors22.csv")
