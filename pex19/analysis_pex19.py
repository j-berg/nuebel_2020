import os
import pandas as pd
import numpy as np
import scipy
from statsmodels.stats.multitest import multipletests
import xpressplot as xp
import matplotlib
import matplotlib.pyplot as plt
from numpy import mean, std
from math import sqrt

def cohen_d(x,y):
    return (
        mean(x) - mean(y)) / sqrt((std(x, ddof=1) ** 2 + std(y, ddof=1) ** 2) / 2.0
    )

"""Import data
"""
data = pd.read_csv(
    os.path.join(os.getcwd(), "pex19", "pex_proteomics.txt"),
    sep='\t',
    index_col=0)
data = data.T
data.columns = [
    'pex19\u0394-1',
    'pex19\u0394-2',
    'pex19\u0394-3',
    'pex19\u0394-6',
    'pex19\u0394-7',
    'pex19\u0394msp1\u0394-1',
    'pex19\u0394msp1\u0394-2',
    'pex19\u0394msp1\u0394-3',
    'pex19\u0394msp1\u0394-4',
    'pex19\u0394msp1\u0394-5'
]

"""Build metadata
"""
meta = pd.DataFrame()
meta[0] = data.columns
meta[1] = [
    'pex19\u0394',
    'pex19\u0394',
    'pex19\u0394',
    'pex19\u0394',
    'pex19\u0394',
    'pex19\u0394msp1\u0394',
    'pex19\u0394msp1\u0394',
    'pex19\u0394msp1\u0394',
    'pex19\u0394msp1\u0394',
    'pex19\u0394msp1\u0394'
]

# PEX15 not measured in dataset
# ATG36 not measured in dataset
pex_list = [
    "PEX13",
    "PEX11",
    "PEX2",
    "PEX14",
    "PEX4",
    "PEX22",
    "PEX17",
    "PEX3",
    "PEX25",
    "PEX29",
    "PEX5",
    "PEX7",
    "PEX6",
    "PEX1",
    "POT1",
    "MDH3",
    "LYS1",
    "YJL185C",
    "SCS7",
    "PEX30",
    "PEX31"
]

"""Heatmap of PEX proteins
"""
sample_colors = {
    'pex19\u0394':'grey',
    'pex19\u0394msp1\u0394':'black'
}

data_scaled, data_labeled = xp.prep_data(data, meta, gene_scale=True)

xp.pca(
    data_scaled,
    meta,
    sample_colors,
    save_fig=os.path.join(os.getcwd(), "pex19", "pca_all.png"))

gene_info = pd.DataFrame()
gene_info[0] = pex_list
gene_info[1] = [
    'Importomer',
    'Importomer',
    'Importomer',
    'Importomer',
    'Importomer',
    'Importomer',
    'Importomer',
    'Importomer',
    'Importomer',
    'Importomer',
    'Cytosolic Cargo Receptors',
    'Cytosolic Cargo Receptors',
    'Membrane Anchored',
    'Membrane Anchored',
    'Matrix Proteins',
    'Matrix Proteins',
    'Matrix Proteins',
    'Matrix Proteins',
    'Matrix Proteins',
    'Other',
    'Other'
]
gene_colors = {
    'Importomer':'#a49284',
    'Cytosolic Cargo Receptors':'#f6efd0',
    'Membrane Anchored':'#41796b',
    'Matrix Proteins':'lightblue',
    'Other':'#d95f02'
}

xp.heatmap(
    data_scaled,
    meta,
    sample_palette = sample_colors,
    gene_info = gene_info,
    gene_palette = gene_colors,
    gene_list = pex_list,
    col_cluster = False,
    row_cluster = False,
    cbar_kws = {'label':'z-score'},
    figsize = (5,6)
)

# Format legends
f = lambda m,c: plt.plot([],[],marker='o', color=c, ls="none", markeredgewidth=0.5, markeredgecolor='black')[0]
handles = [f("s", list(gene_colors.values())[i]) for i in range(len(list(gene_colors.values())))]
first_legend = plt.legend(handles, list(gene_colors.keys()), bbox_to_anchor=(15.5, 1.2705), loc=2, borderaxespad=0., title='Gene Group')

# Add the legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)

g = lambda m,c: plt.plot([],[],marker='o', color=c, ls="none", markeredgewidth=0.5, markeredgecolor='black')[0]
handles_g = [f("s", list(sample_colors.values())[i]) for i in range(len(list(sample_colors.values())))]
plt.legend(handles_g, list(sample_colors.keys()), bbox_to_anchor=(15, 1.2705), loc=1, borderaxespad=0., title='Samples')

# Save and show figure
plt.savefig(
    os.path.join(os.getcwd(), "pex19", "heatmap.png"),
    dpi=1800,
    bbox_inches='tight'
)

"""Volcano Plot
"""
exp = [
    'pex19\u0394msp1\u0394-1',
    'pex19\u0394msp1\u0394-2',
    'pex19\u0394msp1\u0394-3',
    'pex19\u0394msp1\u0394-4',
    'pex19\u0394msp1\u0394-5'
]
cont = [
    'pex19\u0394-1',
    'pex19\u0394-2',
    'pex19\u0394-3',
    'pex19\u0394-6',
    'pex19\u0394-7'
]

data_volcano = data.copy()
data_volcano['log$_2$(Fold Change)'] = np.log2(
    data_volcano[exp].sum(axis=1) / data_volcano[cont].sum(axis=1))
data_volcano['p-val'] = scipy.stats.ttest_ind(
    data_volcano[exp], data_volcano[cont], axis=1)[1]
data_volcano['padj'] = multipletests(
    data_volcano['p-val'].values.tolist(),
    alpha=0.1,
    method='fdr_bh',
    is_sorted=False,
    returnsorted=False)[1]
data_volcano['-log$_1$$_0$(p-value)'] = -1 * np.log10(data_volcano['p-val'].values)

for i, r in data_volcano.iterrows():
    data_volcano.at[i, 'Cohen\'s $\mathit{d}$'] = abs(
        cohen_d(
            data[exp].loc[i].values,
            data[cont].loc[i].values
        )
    )

data.stack().plot.hist(bins=100)
data_volcano['p-val'].hist()
data_volcano['padj'].hist()

data_volcano.to_csv(
    os.path.join(os.getcwd(), "pex19", "data_corrected.txt"),
    sep='\t'
)

p_thresh = -1 * np.log10(0.05)

data_volcano_select = data_volcano.loc[pex_list + ['MSP1']]
data_volcano_label = data_volcano.loc[[
    'MSP1',
    'ERV1',
    'LYS1',
    'PEX13',
    'SKY1',
    'PEX11',
    'PEX2',
    'MDH3',
    'THI80',
    'YKL187C'
]][['log$_2$(Fold Change)', '-log$_1$$_0$(p-value)']]

ax = plt.gca()
data_volcano.plot.scatter(
    'log$_2$(Fold Change)', '-log$_1$$_0$(p-value)',
    c='black',
    s=data_volcano['Cohen\'s $\mathit{d}$'] * 10,
    alpha=0.3,
    grid=False,
    ax=ax)
data_volcano_select.plot.scatter(
    'log$_2$(Fold Change)', '-log$_1$$_0$(p-value)',
    c='red',
    s=data_volcano_select['Cohen\'s $\mathit{d}$'] * 10,
    alpha=1,
    grid=False,
    ax=ax)
l_box = matplotlib.patches.Rectangle((1,p_thresh), 1.5, 8.15, alpha=0.4, color="lightblue", zorder=-10)
r_box = matplotlib.patches.Rectangle((-4.5,p_thresh), 3.5, 8.15, alpha=0.4, color="lightblue", zorder=-10)
ax.set_facecolor('w')
ax.add_patch(l_box)
ax.add_patch(r_box)
ax.axhline(p_thresh, ls='--', color="blue")
ax.axvline(-1, ls='--', color="blue")
ax.axvline(1, ls='--', color="blue")
for index, row in data_volcano_label.iterrows():
    if index == 'PEX13':
        x_pos = -0.2
        y_pos = 0.18
    elif index == 'PEX2':
        x_pos = -0.65
        y_pos = -0.12
    elif index == "ERV1":
        x_pos = -0.65
        y_pos = 0.1
    elif index == "THI80":
        x_pos = -0.65
        y_pos = 0.1
    else:
        x_pos = 0.07
        y_pos = 0.1
    ax.text(
        row[0] + x_pos, row[1] + y_pos,
        str(index),
        horizontalalignment='left',
        size='medium',
        color='black',
        weight='semibold')
ax.text(
    -4.35, 5.85,
    "n = 11",
    horizontalalignment='left',
    size='medium',
    color='blue',
    weight='semibold')
ax.text(
    2, 5.85,
    "n = 8",
    horizontalalignment='left',
    size='medium',
    color='blue',
    weight='semibold')

# Save and show figure
plt.savefig(
    os.path.join(os.getcwd(), "pex19", "volcano_plot_update.png"),
    dpi=1800,
    bbox_inches='tight'
)

g1 = plt.scatter([],[], s=10, marker='o', color='black')
g2 = plt.scatter([],[], s=20, marker='o', color='black')
g3 = plt.scatter([],[], s=50, marker='o', color='black')
plt.legend(
    (g1, g2, g3),
    ['1', '2', '5'],
    bbox_to_anchor=(0, 0),
    loc=1, borderaxespad=1,
    title='Cohen\'s $\mathit{d}$')
plt.savefig(
    os.path.join(os.getcwd(), "pex19", "volcano_plot_legend.png"),
    dpi=1800,
    bbox_inches='tight'
)


data_volcano.loc[
    (data_volcano['log$_2$(Fold Change)'] < -1) &
    (data_volcano['p-val'] < 0.05)
]

data_volcano.loc[
    (data_volcano['log$_2$(Fold Change)'] > 1) &
    (data_volcano['p-val'] < 0.05)
]
