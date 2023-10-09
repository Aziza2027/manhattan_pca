
import streamlit as st
import pandas as pd


st.set_page_config(page_title = "This is a Multipage WebApp")
st.title("Results")
# stm.sidebar.success("Select Any Page from here")


r_stat_p = './code/data/R_stat_p.csv'
txt_file_path = './code/data/confidence.txt'
txt_file_path2 = './code/data/confidence_calc24.txt'
indel_txt_file_path = './code/data/confidence_INDEL_r.txt'
indel_txt_file_path2 = './code/data/confidence_INDEL_calcpcr.txt'
excel_file_path = './code/data/data_Ali.xlsx'

data = pd.read_csv("./code/data/original.csv")
removed = pd.read_csv("./code/data/filled_dropped.csv")
indels = pd.read_csv("./code/data/significant_GENOTYPES.csv")

# --------------------------------------------------------------------

st.markdown('## P values')
my_large_df = pd.read_csv(r_stat_p) #This is your 'my_large_df'

st.dataframe(my_large_df)

@st.cache_data
def convert_df_to_csv(df):
  # IMPORTANT: Cache the conversion to prevent computation on every rerun
  return df.to_csv().encode('utf-8')

st.download_button(
  label="Download data as CSV",
  data=convert_df_to_csv(my_large_df),
  file_name='p_val.csv',
  mime='text/csv',
)

# --------------------------------------------------------------------

st.markdown('## Confidence interval and odds ratios')
st.markdown('### SNPs')
st.markdown('#### R')

# Read the content of the text file
with open(txt_file_path, 'r') as txt_file:
    txt_content = txt_file.read()

st.text_area('', value=txt_content, height=400)

@st.cache_data
def convert_txt_to_bytes(txt_content):
    return txt_content.encode('utf-8')

st.download_button(
    label="Download as TXT",
    data=convert_txt_to_bytes(txt_content),
    file_name='confidence_interval.txt',
    mime='text/plain',
)

st.markdown('#### Calc Pcr 24')

# Read the content of the text file
with open(txt_file_path2, 'r') as txt_file2:
    txt_content2 = txt_file2.read()

st.text_area('', value=txt_content2, height=400)

@st.cache_data
def convert_txt_to_bytes2(txt_content2):
    return txt_content2.encode('utf-8')

st.download_button(
    label="Download as TXT",
    data=convert_txt_to_bytes2(txt_content2),
    file_name='confidence_interval_calc24.txt',
    mime='text/plain',
)

# --------------------------------------------------------------------

st.markdown('### INDELS')
st.markdown('#### R')

# Read the content of the text file
with open(indel_txt_file_path, 'r') as txt_file:
    txt_content = txt_file.read()

st.text_area('', value=txt_content, height=400)

@st.cache_data
def convert_txt_to_bytes(txt_content):
    return txt_content.encode('utf-8')

st.download_button(
    label="Download as TXT",
    data=convert_txt_to_bytes(txt_content),
    file_name='confidence_interval_indel.txt',
    mime='text/plain',
)

st.markdown('#### Calc Pcr 24')

# Read the content of the text file
with open(indel_txt_file_path2, 'r') as txt_file2:
    txt_content2 = txt_file2.read()

st.text_area('', value=txt_content2, height=400)

@st.cache_data
def convert_txt_to_bytes2(txt_content2):
    return txt_content2.encode('utf-8')

st.download_button(
    label="Download as TXT",
    data=convert_txt_to_bytes2(txt_content2),
    file_name='confidence_interval_calc24_indel.txt',
    mime='text/plain',
)

st.markdown('## Significant INDELS')

st.dataframe(indels)

@st.cache_data
def convert_df_to_excel(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return indels.to_csv().encode('utf-8')

st.download_button(
    label="Download data",
    data=convert_df_to_excel(removed),
    file_name='significant_indels.csv',
    mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
)



# --------------------------------------------------------------------

st.markdown('## Genotypes(Imputed)')

# Read the content of the Excel file
my_large_df = pd.read_excel(excel_file_path)

# Display the content of the Excel file
st.dataframe(my_large_df.iloc[:20, :20])

@st.cache_data
def convert_df_to_excel(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')


st.download_button(
    label="Download full data",
    data=convert_df_to_excel(my_large_df),
    file_name='genotype_impotuted.csv',
    mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
)


#######################################


st.markdown('## Genotypes')
st.markdown('#### With missing values, all SNPs and INDELs')

@st.cache_data
def convert_df_to_excel(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return data.to_csv().encode('utf-8')

st.download_button(
    label="Download original data",
    data=convert_df_to_excel(data),
    file_name='original.csv',
    mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
)

#######################################


st.markdown('## Removed genotypes: INDELS, SNPs with 3+ alleles')

st.dataframe(removed)

@st.cache_data
def convert_df_to_excel(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return removed.to_csv().encode('utf-8')

st.download_button(
    label="Download data",
    data=convert_df_to_excel(removed),
    file_name='removed.csv',
    mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
)


#######################################



