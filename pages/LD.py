import subprocess
import pandas as pd

import streamlit as st 
from PIL import Image

import base64
import gzip


df = pd.read_csv('in/sites.txt', delimiter='\t', header=None)
df = df[df[0].str.startswith('N')]
df['chr'] = df[0].str.split('_').str[-1].str.split('.').str[0].astype(int)
df.columns = ['chr', 'pos', 'Chr']
df['com'] = df.chr + df.pos.astype(str)

ids = pd.read_csv('in/SNP_annotated_only_with_id.vcf', delimiter='\t', skiprows=767, usecols=range(3))
ids.columns = ['chr', 'pos', 'id_']
ids['com'] = ids.chr + ids.pos.astype(str)

mer = pd.merge(df, ids, on =['chr', 'pos'], how='inner')
mer[['chr', 'pos', 'id_']].to_csv('./in/in_new.snp', sep='\t', index=False)

def get_ld_block(file_='D_num'):
    d = pd.read_csv('./in/in_new.snp', sep='\t')
    f=gzip.open(f'./out_wbs/{file_}.site.gz','r')
    site=str(f.read())[2:-3]
    df2 = pd.DataFrame([[x.split('\\t')[0], int(x.split('\\t')[1])] for x in site.split('\\n')], columns=['chr', 'pos'])
    rs = pd.merge(df2,d, on=['chr', 'pos'], how='left').id_.to_list()
    f=gzip.open(f'./out_wbs/{file_}.TriangleV.gz','r')
    txt = str(f.read())
    data = [x.split('\\t')for x in txt.split('\\n')[1:-1]]

    new = []
    new2 = []
    for idx, x in enumerate(data):
        add = []
        for j in new:
            add.append(j[idx-1])
        new2.append(add + [1] + x)
        add.extend(x)
        new.append(add)
    else:
        add = []
        for j in new:
            add.append(j[idx])
        new2.append(add + [1])
    print(new2)
    print(rs)

    df = pd.DataFrame(new2, columns=rs, index=rs)

    return df, file_


def calc_LD(start, end, chr, proc):
    print('start')
    chr = mer.chr[mer.Chr==chr].unique()[0]
    command_ = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/D_num   -Region  {chr}:{start}:{end} -OutPng -SeleVar 1 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng -ShowNum'
    command = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/D   -Region  {chr}:{start}:{end} -OutPng -SeleVar 1 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng'
    command2_ = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/R_num   -Region  {chr}:{start}:{end} -OutPng -SeleVar 2 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng -ShowNum'
    command2 = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/R   -Region  {chr}:{start}:{end} -OutPng -SeleVar 2 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng'
    subprocess.run('rm out_wbs/*', shell=True)

    ims = []
    df = []

    if proc[0]:
        res = subprocess.run(command, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res.stdout[-100:]):
            # im = Image.open('out_wbs/D.png')
            ims.append('out_wbs/D.svg')
    if proc[1]:
        res1 = subprocess.run(command_, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res1.stdout[-100:]):
            df.append(get_ld_block('D_num'))
            ims.append('out_wbs/D_num.svg')
    if proc[2]:
        res1 = subprocess.run(command2, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res1.stdout[-100:]):
            ims.append('out_wbs/R.svg')
    if proc[3]:
        res1 = subprocess.run(command2_, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res1.stdout[-100:]):
            df.append(get_ld_block('R_num'))
            ims.append('out_wbs/R_num.svg')

    return ims, df


options = mer.Chr.unique()

st.markdown("### LD for all chromosomes: \nhttps://ld-block.netlify.app/")

# Display the dropdown and store the col2selected value
chr = st.selectbox("Select chromosome:", options)

col1, col2 = st.columns(2)

with col1:
    start = st.number_input("Start:", value=0)
    D_prime = st.checkbox('D-prime', value=True)
    D_prime_ = st.checkbox('D-prime with munbers')
with col2:
    end = st.number_input("End:", value=30000000)
    R2 = st.checkbox('R^2')
    R2_ = st.checkbox('R^2 with munbers')

_,_, col, _,_ = st.columns(5)

with col:
    st.write("")  # Spacer for alignment
    if st.button("Calculate LD"):
        with st.spinner("Calculating..."):
            ims, dfs = calc_LD(start, end, chr, [D_prime, D_prime_, R2, R2_])

def render_svg(svg):
    """Renders the given svg string."""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    # html = r'<img src="data:image/svg+xml;base64,%s" width="1000" height="1000" />' % b64
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    print('=======================================')
    st.write(html, unsafe_allow_html=True)        

def display(file):
    f = open(file,"r")
    lines = f.readlines()
    line_string=''.join(lines)

    render_svg(line_string)

@st.cache_data
def convert_df_to_excel(df_tmp):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df_tmp.to_csv().encode('utf-8')


try:
    for im in ims:
        print(im)
        display(im)
        
    res = mer[((mer.Chr==chr)&(mer.pos>=start)&(mer.pos<=end))][['chr', 'pos', 'id_']].copy()
    
    
    for df_ in dfs:
        st.markdown('\n\n\n')
        st.markdown(f'### {df_[1]}')
        st.dataframe(df_[0])
        
        st.download_button(
            label="Download data",
            data=convert_df_to_excel(df_[0]),
            file_name=f'{df_[1][0]}_data.csv',
            mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
        )
                
    
    st.markdown('### Data:')
    if res.shape[0]:
        st.dataframe(res)


except Exception as e:
    print(e)
    pass
