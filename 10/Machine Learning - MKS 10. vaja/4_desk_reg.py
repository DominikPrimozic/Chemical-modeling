# Uvoz knjiznic in modulov
import pandas as pd           
import numpy as np            
from rdkit import Chem, RDLogger  
from rdkit.Chem import Descriptors 
from sklearn.linear_model import LinearRegression   
from sklearn.metrics import mean_squared_error, r2_score   
import matplotlib.pyplot as plt      

# 1) utisa (ne prikaze) RDKit parser napake 
RDLogger.DisableLog('rdApp.*')     

def load_and_clean(path):       # funkcija za branje .csv datoteke
    """Preberi CSV, pretvori melting_point v numeric in odstrani neveljavne vrstice.""" 
    df = pd.read_csv(path)      
    df['melting_point'] = pd.to_numeric(df['boiling_point'], errors='coerce')
    return df.dropna(subset=['melting_point']).reset_index(drop=True)  # prestevilcenje vrstic, v primeru brisanja

    # funkcija compute_descriptors kot argument rabi pandas DataFrame, pricakuje SMILES stolpec
def compute_descriptors(df):
    """Izracunaj osnovne RDKit descriptorje iz SMILES."""  
    # za vsak kljucnik ustvarimo NumPy polje dolzine len(df), 1 element na molekulo, napolnjeno z niclami
    descs = {                               
        'MolWt':    np.zeros(len(df)),      # vektor nicel (prikladna inicializacija) dolzine enake stevilu vrstic v df
        'LogP':     np.zeros(len(df)),      # s tem ze alociramo velikost glede na df
        'TPSA':     np.zeros(len(df)),
        'NumRot':   np.zeros(len(df)),
    }
    valid = []          # vpise dejanskih deskriptorjev (seznam True/Falsev)
    # zanka prek vseh SMILES, vrne zaporedje parov (indeks, vrednost)
    for i, smi in enumerate(df.SMILES):
        mol = Chem.MolFromSmiles(smi)  # iz SMILES proba RDKit narest molekulo (graf atomov in vezi)
        if mol is None:                
            valid.append(False)        # v seznam valid vnese False in vrstice ne bomo uporabili
            continue                   # preskoci na naslednjo vrstico
        valid.append(True)             # nastavi valid na True, ce je SMILES okej
        descs['MolWt'][i]  = Descriptors.MolWt(mol)      # za vsak veljaven mol izracunamo 4 deskriptorje
        descs['LogP'][i]   = Descriptors.MolLogP(mol)    
        descs['TPSA'][i]   = Descriptors.TPSA(mol)       
        descs['NumRot'][i] = Descriptors.NumRotatableBonds(mol)  
    # po koncu zanke, descs vsebuje za vsak i tocne deskriptorje, invalidni so ostali nicle
    df_desc = pd.DataFrame(descs)    
    df_desc['valid'] = valid     # dodan nov stolpec valid, kjer je za vsak i True ali False
    df_all = pd.concat([df.reset_index(drop=True), df_desc], axis=1)
    #df_all_valid filtrira vrstice, samo True ostanejo
    #drop(columns) izbrise stolpec valid, k ga ne nucamo vec
    #reset_index nanovo indeksira vrstice od 0 navzdol, da ni vrzeli
    return df_all[df_all.valid].drop(columns=['valid']).reset_index(drop=True)

# --- 1. NALOZI IN POCISTI NEUPORABNO ---
df_train = load_and_clean('train.csv')  # df_train uceni primeri (80% podatkov)
df_test  = load_and_clean('test.csv')   # df_test testni primeri (20% podatkov)

# --- 2. IZRACUN DESKRIPTORJEV ---
# iterira cez SMILES, naredi molekulo z RDKit
# ob veljavnem SMILES izracuna deskriptorje, vrne razsirjen DataFrame
# poleg originalnih stolpcev se stiri nove numericne stolpce kot vhod za regresijo
df_train = compute_descriptors(df_train)
df_test  = compute_descriptors(df_test)

# --- 3. PRIPRAVI X in y ZA REGRESIJO ---
# pripravis vhodne (X) in tarcno (y) spremenljivko/e
# selektivna izbira stolpcev, v X zgolj deskriptorji brez SMILES in Ttalisc
# tarca y za ucno mnozico pa privzamemo talisca
feature_cols = ['MolWt','LogP','TPSA','NumRot']
X_train, y_train = df_train[feature_cols], df_train['melting_point']
X_test,  y_test  = df_test[feature_cols],  df_test['melting_point']

# --- 4. Fit model ---
# Xtrain in ytrain damo v model, da ga naucimo parametrov beta
model = LinearRegression().fit(X_train, y_train)
# Xtest nato uporabimo v model.predict da napovemo talisca in jih primerjamo z testnimi (y_test)
y_pred = model.predict(X_test)       
residuals = y_test - y_pred

# --- 5. EVALVACIJA KAKOVOSTI MODELA ---
# Izracun MSE med pravimi in napovedanimi vrednostmi
# koren MSE, enaka enota kot talisce (°C)
mse  = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
r2   = r2_score(y_test, y_pred)        # koliko variance v podatkih model razlozi

print("Uporabljeni descriptorji:", feature_cols)
print(f"Stevilo ucnih primerov: {len(y_train)}, testnih: {len(y_test)}") 
print(f"RMSE: {rmse:.4f}")             
print(f"R^2  : {r2:.4f}")

# --- 6. IZPIS NEKAJ PRVIH PRIMEROV ZA VPOGLED V NAPOVED ---
print("\nPrvi primeri (Actual vs Predicted):")    
print(pd.DataFrame({                              # nastavi nov DataFrame s 4 stolpci
    'SMILES':   df_test['SMILES'],
    'Actual':   y_test,
    'Predicted': np.round(y_pred, 2),             # zaokrozitev na 2 decimalki
    'Residual': np.round(residuals, 2)
}).head(10).to_string(index=False))               # izbere le prvih 10 vrstic in predstavi kot besedilo

# --- 7. RAZTRESENI GRAF: dejanske vs. napovedane vrednosti ---
plt.figure(figsize=(6,6))      
plt.scatter(y_test, y_pred, alpha=0.6)  # os x=pravo talisce os y=napovedano
# mn = min(y_test.min(), y_pred.min()) , mx = max(y_test.max(), y_pred.max())
# izracunas min in max vrednost izmed vseh napovedanih in dejanskih, da mas cel razpon
mn, mx = min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())
plt.plot([mn, mx], [mn, mx], 'k--', lw=1)  
plt.xlabel('Dejansko melting_point')        
plt.ylabel('Napovedano melting_point')
plt.title('Dejansko vs. Napovedano')
plt.tight_layout()                  
plt.show()                          

# --- 8. HISTOGRAM PREOSTANKOV ---
plt.figure(figsize=(6,4))
plt.hist(residuals, bins=20, edgecolor='k', alpha=0.7) 
plt.xlabel('Residuali (Actual − Predicted)')
plt.ylabel('Stevilo molekul')
plt.title('Histogram napak')
plt.tight_layout()
plt.show()
