from selenium import webdriver
from selenium.webdriver.common.keys import Keys

import pandas as pd

import sys

output_file = "./data/confidence_calc24.txt"

# data_p = '/home/bio/projects/manhattan_pca/code/data/for Aziza.xlsx'
# p_val = None

data_p = '/home/bio/projects/manhattan_pca/code/data/data_Ali.xlsx'
p_val = '/home/bio/projects/manhattan_pca/code/data/R_stat_p.csv'


df = pd.read_excel(data_p)

# Choose the appropriate web driver (e.g., Chrome, Firefox)
driver = webdriver.Chrome()  # You need to have the Chrome driver installed

# Navigate to the website
website_url = "https://calc.pcr24.ru/index.php"
driver.get(website_url)


if p_val:
    tmp = pd.read_csv(p_val)
    rss = tmp[tmp.codominant<=0.05].iloc[:,0].values
else:
    rss = df.columns[1:]

def get_res(values):
    try:
        def put(k,v):
            # Find and input value A
            key = f"input[name='{k}']"
            input_a = driver.find_element("css selector", key)
            input_a.clear()  # Clear any existing value
            input_a.send_keys(str(v))

        for k, v in values.items():
            put(k,v)

        calculate_button = driver.find_element("css selector", "input[name='submit']")
        calculate_button.click()

        driver.implicitly_wait(3)  # Wait for up to 10 seconds

        # Get the page source after calculations
        page_source = driver.page_source
        dfs = pd.read_html(page_source)
        df_ = []
        for i, name in zip([1,2,4,5], ['Allele', 'Codominant', 'Dominant', 'Recessive']):
            DF = dfs[i].set_index(dfs[i].columns[0])
            DF.columns = DF.iloc[0]
            DF = DF[1:]
            df_.append(name + '\n-----------------\n' +DF.to_string(col_space=10))

        return '\n\n'.join(df_)

    except Exception as e:
        print("An error occurred:", str(e))



def get_stats(d):

    if len(set(list(d.index[0]))) > 1:
        my_list = list(d.index)
        my_list[0], my_list[1] = my_list[1], my_list[0]
        d = d.loc[my_list, :]

    a = d.index.str.cat()
    a,b = a[0], a[-1]

    try:
        sg3 = d.iloc[2,0]
        cg3 = d.iloc[2,1]
    except:
        sg3, cg3 = 0, 0
        
    if a == b:
        print(d)
        b = d.index.str.cat()[-2]
        print('ERROR', a, b, sg3, cg3)

        # sys.exit()

    values = {
        'a1':a, 
        'a2': b, 
        'sg1': d.iloc[0,0], 
        'sg2': d.iloc[1,0], 
        'sg3': sg3, 
        'cg1': d.iloc[0,1], 
        'cg2': d.iloc[1,1], 
        'cg3': cg3}


    res = get_res(values)
    return res

results = ''
for rs in rss:
    
    res = df.pivot_table(index=rs, columns='Status', aggfunc='size', fill_value=0)
    res = rs + '\n' + get_stats(res) + '\n\n============================================================\n\n'
    results += res


# results
# Close the browser
driver.quit()
with open(output_file, 'w') as f:
    f.write(results)

