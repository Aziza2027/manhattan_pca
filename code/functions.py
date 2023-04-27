import numpy as np
import pandas as pd

from scipy.stats import chi2_contingency
from scipy.stats.contingency import odds_ratio

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import seaborn as sns
import matplotlib.pyplot as plt

from adjustText import adjust_text

from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from sklearn.preprocessing import OrdinalEncoder

from ast import literal_eval


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


def preprocess_data(rs, chr, pos, p_val):
    df = pd.concat([rs, chr, pos, p_val], axis=1)
    df.columns = ['rsid', 'CHR', 'POS','p_val']
    df = df.sort_values(['CHR', 'POS'], ascending=[True, True])
    df['neg_p_val']  = - np.log10(df['p_val'])
    running_pos = 0
    cumulative_pos = []

    for _, group_df in df.groupby('CHR'):  
        cumulative_pos.append(group_df['POS'] + running_pos)
        running_pos += group_df['POS'].max()
        
    df['cumulative_pos'] = pd.concat(cumulative_pos)
    df.sort_values(by='CHR', inplace=True)

    return df


def get_chr(rsR, rs, chr, pos):
    test = pd.concat([rs, chr, pos], axis=1)
    test.columns = ['rs', 'chr', 'pos']
    merged = pd.merge(test, rsR)
    return merged.chr, merged.pos


def get_manhattan(rs, chr, pos, p_val, title):
    df = preprocess_data(rs, chr, pos, p_val)

    # Copy data
    my_data = df.copy()

    # Thresholds
    threshold = 0.05  # Red line
    display = 0.05  # p-value below this threshold will be displayed
    neg_l_t = -np.log10(threshold)
    neg_l_d = -np.log10(display)

    # Create figure
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # Plot scatter
    fig.add_trace(
        go.Scatter(
            x=my_data['cumulative_pos'], 
            y=my_data['neg_p_val'],
            mode='markers',
            marker=dict(size=6, color=my_data['CHR'], colorscale=px.colors.qualitative.Dark24),
            showlegend=False,
            hovertext=df.rsid.values,
            hoverinfo="text",
        ),
        secondary_y=False
    )

    # Set x-axis labels and tick positions
    fig.update_xaxes(
        tickmode='array',
        tickvals=my_data.groupby('CHR')['cumulative_pos'].median(),
        ticktext=my_data['CHR'].unique(),
        title_text='Chromosome'
    )

    # Set y-axis labels
    fig.update_yaxes(
        title_text='—Log10(p-value)',
        secondary_y=False
    )

    # Add horizontal line at threshold
    fig.add_shape(
        type="line",
        x0=my_data['cumulative_pos'].min(),
        y0=neg_l_t,
        x1=my_data['cumulative_pos'].max(),
        y1=neg_l_t,
        line=dict(color="red", width=1, dash='dash'),
        secondary_y=False,
    )

    # Add text label for threshold
    fig.add_annotation(
        x=-1,
        y=neg_l_t-0.11,
        text=f"p-value = {threshold}",
        showarrow=False,
        font=dict(size=10, family='italic')
    )

    # Add annotations for significant data points
    sig_points = my_data[my_data['neg_p_val'] > neg_l_d]
    # print(sig_points.rsid.values)
    annotations = []
    for i, row in sig_points.iterrows():
        annotations.append(
            dict(
                x=row['cumulative_pos'], 
                y=row['neg_p_val'], 
                text=row['rsid'], 
                # text="outlier",
                showarrow=False,
                yshift=9,
                xshift=28,
                font=dict(size=12, family='italic')
            )
        )
    layout = go.Layout(
        # title="Box plot of Amino acid deletions",    
        plot_bgcolor="#FFFFFF",
        barmode="stack",
        xaxis=dict(
            # domain=[0, 0.5],
            # title="substitutions",
            linecolor="#BCCCDC",
        ),
        yaxis=dict(
            # title="frequency",
            linecolor="#BCCCDC"
        )
        )

    fig.update_layout(
        title_text=title,
        # xaxis=dict(showgrid=True, zeroline=False),
        # yaxis=dict(showgrid=True, zeroline=False),
        # hovermode='closest',
        margin=dict(l=50, r=50, b=50, t=50),
        annotations=annotations,
        showlegend=False
    )
    fig.update_layout(layout)

    fig.show()


def get_manhattan2(rs, chr, pos, p_val):

    df = preprocess_data(rs, chr, pos, p_val)

    my_data = df.copy()

    threshold = 0.05 # red line  
    display = 0.05 # p-value below this threshold will be displayed
    neg_l_t = -np.log10(threshold)
    neg_l_d = -np.log10(display)

    g = sns.relplot(
        data = my_data,
        x = 'cumulative_pos',
        y = 'neg_p_val',
        aspect = 4,
        hue = 'CHR',
        # palette = ['grey', 'black'] * 11,
        palette = 'Set1',
        linewidth=0,
        s=20,
        legend=None
    )

    g.ax.set_xlabel('Chromosome')
    g.ax.set_ylabel('—Log10(p-value)')

    g.ax.set_xticks(my_data.groupby('CHR')['cumulative_pos'].mean())

    g.ax.set_xticklabels(df['CHR'].unique())
    g.ax.axhline(neg_l_t, ls='--', linewidth=1, color='red')
    g.ax.text(-1, neg_l_t-0.11, f'p-value = {threshold}', fontsize=10, va='bottom', fontstyle='italic')

    # g.fig.suptitle('GWAS plot showing association between SNPs on autosomes and speeding')

    annotations = my_data[my_data['neg_p_val'] > neg_l_d].apply(lambda p : g.ax.annotate(p['rsid'], (p['cumulative_pos'], p['neg_p_val'])), axis=1).to_list()

    adjust_text(annotations, arrowprops = {'arrowstyle' : '-', 'color' : 'blue'})
    plt.plot()
    # plt.savefig('../visualizations/fig.jpg')


def to_array(df):
    def convert_to_array(col):
        """
        Convert genotypic information to an array.
        """
        w, m = col.name[-3:].split('>')
        wild, heterozygous, mutant = w+w, w+m, m+m
    
        replace_values = {
            wild : '[0, 0]',
            heterozygous : '[0, 1]',
            mutant : '[1, 1]'
        }

        res = col.replace(replace_values)
        try:
            res = res.apply(lambda x: np.array(literal_eval(str(x))))
        except:
            print(res.unique(), col.name)
            
        return res
        
    result = df.apply(convert_to_array)

    return result

def factorize(df):
    df = df.apply(lambda x: pd.factorize(x)[0])
    return df


def PCA3d(data, groups, title='3D PCA'):
    X_reduced = PCA(n_components=3).fit_transform(data)
    # Create a dataframe with the reduced data
    df = pd.DataFrame(X_reduced, columns=['PC1', 'PC2', 'PC3'])
    df['Groups'] = groups
    # Create the scatter plot with Plotly
    fig = px.scatter_3d(df, x='PC1', y='PC2', z='PC3', color='Groups', 
                        symbol='Groups', opacity=0.9,
                        width=800, height=800)

    # Set the title and axis labels
    fig.update_layout(title=title, 
                      template='plotly_white', 
                    scene=dict(xaxis_title='1st eigenvector', yaxis_title='2nd eigenvector', 
                                zaxis_title='3rd eigenvector'))
    fig.update_traces(marker_size=6)
    fig.update_coloraxes(showscale=False)
    return fig

def PCA2d(data, groups, title = '2D PCA'):
    layout = go.Layout(
        plot_bgcolor="#FFFFFF",
        barmode="stack",
        xaxis=dict(
            linecolor="#BCCCDC",
        ),
        yaxis=dict(
            linecolor="#BCCCDC"
        )
    )


    X_reduced = PCA(n_components=2).fit_transform(data)
    # Create a dataframe with the reduced data
    df = pd.DataFrame(X_reduced, columns=['PC1', 'PC2'])
    df['Groups'] = groups
    # Create the scatter plot with Plotly
    fig = px.scatter(df[::-1], x='PC1', y='PC2', color='Groups', 
                        symbol='Groups', opacity=0.9,
                        width=700, height=500)

    # Set the title and axis labels
    fig.update_layout(title=title, 
                      template='plotly_white', 
                    #   scene=dict(xaxis_title='1st eigenvector', yaxis_title='2nd eigenvector', 
                    #              zaxis_title='3rd eigenvector')
                                )
    # fig.update_layout(template='plotly_white', )
    fig.update_traces(marker_size=8)
    fig.update_coloraxes(showscale=False)
    return fig

