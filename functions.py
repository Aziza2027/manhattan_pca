import pandas as pd
import numpy as np

from sklearn.impute import KNNImputer
from sklearn.preprocessing import OrdinalEncoder

def fill(arr):
    imputer = KNNImputer(n_neighbors=10)
    print(arr.shape)
    d = imputer.fit_transform(arr)
    print(d.shape)
    return d


def fill_numerical(Data):
    filled = fill(Data.copy())
    # print(len(Data.columns))
    # print(filled.shape)
    # print(Data.copy().shape)
    return pd.DataFrame(filled, columns=Data.columns)


def fill_categorical(Cat_data):
    # encoder = OrdinalEncoder(handle_unknown='use_encoded_value', unknown_value=-1)
    encoder = OrdinalEncoder()

    encoded_data = encoder.fit_transform(Cat_data.values)

    encoded_data[encoded_data == -1] = np.nan
    data = fill(encoded_data)

    encoded_data = np.round(data)

    # Inverse transform the encoded data to get the original data
    df_decoded = pd.DataFrame(encoder.inverse_transform(encoded_data), columns=Cat_data.columns)

    return df_decoded