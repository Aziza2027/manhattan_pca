import pandas as pd
import numpy as np

from sklearn.impute import KNNImputer
from sklearn.preprocessing import OrdinalEncoder

def fill_numerical(Data):
    imputer = KNNImputer(n_neighbors=5)
    d = imputer.fit_transform(Data)
    return pd.DataFrame(d, columns=Data.columns)


def fill_categorical(Cat_data):
    print('here')
    # encoder = OrdinalEncoder(handle_unknown='use_encoded_value', unknown_value=-1)
    encoder = OrdinalEncoder()
    print('here2')
    print(Cat_data)

    encoded_data = encoder.fit_transform(Cat_data.values)
    print('here3')


    encoded_data[encoded_data == -1] = np.nan
    print(encoded_data)
    data = fill_numerical(encoded_data)
    print('here4')

    encoded_data = np.round(data.values)
    print(encoded_data)

    # Inverse transform the encoded data to get the original data
    df_decoded = pd.DataFrame(encoder.inverse_transform(encoded_data), columns=Cat_data.columns)
    print(df_decoded)

    return df_decoded, encoded_data