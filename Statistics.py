
import streamlit as st
import pandas as pd


st.set_page_config(page_title = "This is a Multipage WebApp")
st.title("This is the Home Page Geeks.")
# stm.sidebar.success("Select Any Page from here")

# --------------------------------------------------------------------

st.markdown('## P values')
my_large_df = pd.read_csv('./code/data/R_stat_p.csv') #This is your 'my_large_df'

st.dataframe(my_large_df)

@st.cache_data
def convert_df_to_csv(df):
  # IMPORTANT: Cache the conversion to prevent computation on every rerun
  return df.to_csv().encode('utf-8')

st.download_button(
  label="Download data as CSV",
  data=convert_df_to_csv(my_large_df),
  file_name='large_df.csv',
  mime='text/csv',
)

# --------------------------------------------------------------------

st.markdown('## Confidence interval and odds ratios')

txt_file_path = './code/data/confidence.txt'

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
    file_name='text_file.txt',
    mime='text/plain',
)

# --------------------------------------------------------------------

st.markdown('## Genotypes(Imputed)')

excel_file_path = './code/data/data_Ali.xlsx'

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
    file_name='CI.csv',
    mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
)


#######################################


st.markdown('## Genotypes(Original, with missing values)')

data = pd.read_csv("./code/data/original.csv")


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


st.markdown('## Removed genotypes(3+ alleles)')

data = pd.read_csv("./code/data/filled_dropped.csv")

st.dataframe(data)

@st.cache_data
def convert_df_to_excel(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return data.to_csv().encode('utf-8')

st.download_button(
    label="Download data",
    data=convert_df_to_excel(data),
    file_name='removed.csv',
    mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
)
