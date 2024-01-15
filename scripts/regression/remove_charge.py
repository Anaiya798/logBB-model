import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from rdkit.Chem import SaltRemover, MolFromSmiles, Draw, GetFormalCharge, MolToSmiles
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdFingerprintGenerator

un = rdMolStandardize.Uncharger() # remove extra charge


filename = 'final_dataset_logbb.csv'
df = pd.read_csv(filename)

print(df.shape)

df['Agglomeration'] = df['SMILES'].apply(lambda x: '.' in x)
df['SMILES_clear'] = df['SMILES'].apply(lambda x: max(x.split('.'), key=len))
df['ROMol'] = df['SMILES_clear'].apply(lambda x: MolFromSmiles(x))
print(df.shape)
df = df[~df["ROMol"].isna()]
print(df.shape)

df['FORMAL_CHARGE_ch'] = df['SMILES_clear'].apply(lambda x: GetFormalCharge(MolFromSmiles(x)))
df['SMILES_uncharge'] = df['ROMol'].apply(lambda x: MolToSmiles(un.uncharge(x), kekuleSmiles=True))
df['FORMAL_CHARGE_unch'] = df['SMILES_uncharge'].apply(lambda x: GetFormalCharge(MolFromSmiles(x)))

df.drop(labels=['SMILES_clear', 'ROMol', 'FORMAL_CHARGE_ch'], inplace=True, axis=1)

print(df['FORMAL_CHARGE_unch'].describe())

df.reset_index(drop=True, inplace=True)

df.to_csv(f'{filename[:-4]}_mod.csv', index=False)

