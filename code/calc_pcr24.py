import pandas as pd

from functions import CalcPcr

output_file = "./code/data/confidence_calc24.txt"

# data_p = '/home/bio/projects/manhattan_pca/code/data/for Aziza.xlsx'
# p_val = None

data_p = '/home/bio/projects/manhattan_pca/code/data/data_Ali.xlsx'
p_val = '/home/bio/projects/manhattan_pca/code/data/R_stat_p.csv'


df = pd.read_excel(data_p)

if p_val:
    tmp = pd.read_csv(p_val)
    rss = tmp[tmp.codominant<=0.05].iloc[:,0].values
else:
    rss = df.columns[1:]


driv = CalcPcr()
results = driv.calc_stats(df[['Status'] + list(rss)])

with open(output_file, 'w') as f:
    f.write(results)

