import pandas as pd
import numpy as np

from scipy.stats import chi2_contingency
from scipy.stats.contingency import odds_ratio

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

def calculate_OR(case_group, control_group, wild, heterozygous, mutant,model, confidence_level = 0.95):
    # print('case')
    # print(sum(case_group == wild), wild)
    # print(sum(case_group == heterozygous), heterozygous)
    # print(sum(case_group == mutant), mutant)
    # print('con')
    # print(sum(control_group == wild), wild)
    # print(sum(control_group == heterozygous), heterozygous)
    # print(sum(control_group == mutant), mutant)

    def get_stats(cross_tabs):
        
        stats = []
        for cross_tab in cross_tabs:

            cross_tab[cross_tab==0] = 1
            res = odds_ratio(cross_tab)
            OR = round(res.statistic, 4)
            CI = res.confidence_interval(confidence_level=confidence_level)
            lower_ci, upper_ci = round(CI.low, 4), round(CI.high, 4)

            stats.extend([OR, (lower_ci, upper_ci)])
        
        return stats

    def get_crosstabs():
        c_tabs = []
        if model == 'dominant':
            r1 = [sum(control_group == heterozygous) + sum(control_group == mutant), sum(control_group == wild)]
            r2 = [sum(case_group == heterozygous) + sum(case_group == mutant), sum(case_group == wild)]
            c_tab = np.array([r1, r2])
            c_tab2 = np.array([r2, r1])
            c_tabs.extend([c_tab, c_tab2])

        elif model == 'recessive':
            r1 = [sum(case_group == mutant), sum(case_group == wild) + sum(case_group == heterozygous)]
            r2 = [sum(control_group == mutant), sum(control_group == wild) + sum(control_group == heterozygous)]

            c_tab = np.array([r1, r2])
            c_tab2 = np.array([r2, r1])

            c_tabs.extend([c_tab2, c_tab])
        
        elif model == 'codominant':
            r1 = [sum(case_group == mutant), sum(case_group == wild) + sum(case_group == heterozygous)]
            r2 = [sum(control_group == mutant), sum(control_group == wild) + sum(control_group == heterozygous)]
            r3 = [sum(case_group == wild), sum(case_group == heterozygous) + sum(case_group == mutant)]
            r4 = [sum(control_group == wild), sum(control_group == heterozygous) + sum(control_group == mutant)]
            r5 = [sum(case_group == heterozygous), sum(case_group == wild) + sum(case_group == mutant)]
            r6 = [sum(control_group == heterozygous), sum(control_group == wild) + sum(control_group == mutant)]

            c_tab = np.array([r1, r2])
            c_tab2 = np.array([r3, r4])
            c_tab3 = np.array([r5, r6])

            c_tabs.extend([c_tab2, c_tab3, c_tab])
        else: # allele
            cases = pd.Series(list(''.join(case_group)))
            conls = pd.Series(list(''.join(control_group)))
            
            r1 = [sum(cases == wild[0]), sum(cases == mutant[0])]
            r2 = [sum(conls == wild[0]),sum(conls == mutant[0])]

            c_tab = np.array([r1, r2])
            c_tab2 = np.array([r2, r1])

            c_tabs.extend([c_tab2, c_tab])

        return c_tabs 
    
    
    cross_tabs = get_crosstabs()
    results = get_stats(cross_tabs)

    return results

def calculate_P(case, gen):
    
    c_tab = pd.crosstab(case,gen)
    p_gen = chi2_contingency(c_tab).pvalue

    col1 = pd.concat([case,case]).reset_index(drop=True)
    col2 = pd.concat([gen.str[0],gen.str[1]]).reset_index(drop=True)
    c_tab2 = pd.crosstab(col1,col2).iloc[[1,0]]
    p_all = chi2_contingency(c_tab2).pvalue

    return p_all,p_gen
