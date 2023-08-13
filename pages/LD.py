import subprocess
import pandas as pd

import streamlit as st 
from PIL import Image


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

def calc_LD(start, end, chr, proc):
    print('start')
    chr = mer.chr[mer.Chr==chr].unique()[0]
    command_ = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/D_num   -Region  {chr}:{start}:{end} -OutPng -SeleVar 1 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng -ShowNum'
    command = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/D   -Region  {chr}:{start}:{end} -OutPng -SeleVar 1 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng'
    command2_ = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/R_num   -Region  {chr}:{start}:{end} -OutPng -SeleVar 2 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng -ShowNum'
    command2 = f'./LDBlockShow -InVCF in/SNP_annotated_only_with_id.vcf -OutPut  out_wbs/R   -Region  {chr}:{start}:{end} -OutPng -SeleVar 2 -NoShowLDist 1000000000 -SpeSNPName in/in_new.snp -OutPng'
    subprocess.run('rm out_wbs/*', shell=True)

    ims = []

    if proc[0]:
        res = subprocess.run(command, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res.stdout[-100:]):
            im = Image.open('out_wbs/D.png')
            ims.append(im)
    if proc[1]:
        res1 = subprocess.run(command_, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res1.stdout[-100:]):
            im = Image.open('out_wbs/D_num.png')
            ims.append(im)
    if proc[2]:
        res1 = subprocess.run(command2, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res1.stdout[-100:]):
            im = Image.open('out_wbs/R.png')
            ims.append(im)
    if proc[3]:
        res1 = subprocess.run(command2_, shell=True,stdout = subprocess.PIPE)
        if 'done' in str(res1.stdout[-100:]):
            im = Image.open('out_wbs/R_num.png')
            ims.append(im)

    return ims

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
            ims = calc_LD(start, end, chr, [D_prime, D_prime_, R2, R2_])
        

try:
    for im in ims:
        st.image(im)
    res = mer[((mer.Chr==chr)&(mer.pos>=start)&(mer.pos<=end))][['chr', 'pos', 'id_']].copy()
    st.markdown('### Data:')
    if res.shape[0]:
        st.dataframe(res)
except:
    pass
