import pandas as pd

path = '/home/bio/Downloads/GWAS_analysis-20230810T033758Z-001/GWAS_analysis/variants/'

DATA_PATH = path + 'SNP_with_id.vcf'

DATA_PATH_FILLED = path + 'filled.csv'

def read_vcf(data_path):
    df = pd.read_csv(DATA_PATH, delimiter=' ')
    df = df[['#[1]CHROM','[2]POS', '[3]ID']]
    df.columns = ['chr','pos', 'rs']
    chr = df.chr.str.split('.').str[0].str.split('_').str[1].astype(int)
    # print(chr)

    df.chr = chr
    # df[df.chr<24]

    return df

df = read_vcf(DATA_PATH)

df.to_csv('./code/data/rs_info.csv', index=False)

df = pd.read_csv(DATA_PATH_FILLED)

df['Status'] = df.iloc[:,0].apply(lambda x: int('k' not in x))

df.to_excel('./code/data/data_Ali.xlsx', engine='xlsxwriter', index=False)

print(df.Status.value_counts())
