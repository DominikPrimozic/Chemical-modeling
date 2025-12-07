
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt

# 1) Nalozi in pocisti CSV (samo talisce kot tarca)
def load_and_clean(path, target):
    df = pd.read_csv(path)
    df[target] = pd.to_numeric(df[target], errors='coerce')
    return df.dropna(subset=[target]).reset_index(drop=True)

# nastavi tarcno spremenljivko
target = 'melting_point'

# 2) Nalozi ucni in testni sklop
df_train = load_and_clean('train.csv', target)
df_test  = load_and_clean('test.csv',  target)

# 3) Pripravi X (feature_cols) in y (tarca)
feature_cols = [
    'nHBAcc', 'nHBDon', 'nRot', 'RotRatio',
    'SLogP', 'SMR', 'MW', 'AMW',
    'WPath', 'WPol', 'Zagreb1', 'Zagreb2',
    'mZagreb1', 'mZagreb2'
]

X_train, y_train = df_train[feature_cols], df_train[target]
X_test,  y_test  = df_test[feature_cols],  df_test[target]

# 4) Fit model in napoved
model  = LinearRegression().fit(X_train, y_train)
y_pred = model.predict(X_test)
residuals = y_test - y_pred

# 5) Evalvacija
mse  = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
r2   = r2_score(y_test, y_pred)

print("Uporabljeni deskriptorji:", feature_cols)
print(f"Stevilo ucnih primerov: {len(y_train)}, testnih: {len(y_test)}")
print(f"RMSE: {rmse:.4f}")
print(f"R^2  : {r2:.4f}")

# 6) Prikaz prvih 10 primerov
print("\nPrvi primeri (Actual vs. Predicted):")
print(pd.DataFrame({
    'SMILES':    df_test['SMILES'],
    'Actual':    y_test,
    'Predicted': np.round(y_pred, 2),
    'Residual':  np.round(residuals, 2)
}).head(10).to_string(index=False))

# 7) Raztreseni graf: dejansko vs. napovedano
plt.figure(figsize=(6,6))
plt.scatter(y_test, y_pred, alpha=0.6)
mn, mx = min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())
plt.plot([mn, mx], [mn, mx], 'k--', lw=1)
plt.xlabel(f'Dejansko {target}')
plt.ylabel(f'Napovedano {target}')
plt.title(f'Dejansko vs. Napovedano ({target})')
plt.tight_layout()
plt.show()

# 8) Histogram ostankov
plt.figure(figsize=(6,4))
plt.hist(residuals, bins=20, edgecolor='k', alpha=0.7)
plt.xlabel('Residuali (Actual − Predicted)')
plt.ylabel('Stevilo molekul')
plt.title(f'Histogram napak ({target})')
plt.tight_layout()
plt.show()
