import pandas as pd
import numpy as np

from functions import fill_categorical, fill_numerical
from scipy.stats.contingency import odds_ratio

# Load dataset
data = pd.read_csv('../data/join.csv')

# Initialize list to store imputed columns
cat_cols = [col for col in data.columns if data[col].dtype == 'object']
num_cols = [col for col in data.columns if data[col].dtype != 'object']
cat_data = data[cat_cols]
num_data = data[num_cols]
exclude = list(num_data.isna().sum().sort_values()[::-1][:30].index)
num_data = num_data.drop(columns=exclude)

c_data = fill_categorical(cat_data)
n_data = fill_numerical(num_data)

final = pd.concat([n_data, c_data], axis=1)
final.to_csv('data/join_filled.csv', index=False)