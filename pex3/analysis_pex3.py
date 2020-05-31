import pandas as pd
import numpy as np
import xpressplot as xp
import matplotlib
import matplotlib.pyplot as plt

"""Import data
"""
data = pd.read_csv(
    '/Users/jordan/Desktop/nuebel_2020/pex3/pex3_scaled.txt',
    sep='\t',
    index_col=0)
data.columns = [
    'pex19\u0394-1',
    'pex19\u0394-2',
    'pex19\u0394msp1\u0394-1',
    'pex19\u0394msp1\u0394-2',
    'pex3\u0394-1',
    'pex3\u0394-2',
    'pex3\u0394-3',
    'pex3\u0394msp1\u0394-1',
    'pex3\u0394msp1\u0394-2',
    'pex3\u0394msp1\u0394-3',
]

data = data[[
    'pex3\u0394-1',
    'pex3\u0394-2',
    'pex3\u0394-3',
    'pex3\u0394msp1\u0394-1',
    'pex3\u0394msp1\u0394-2',
    'pex3\u0394msp1\u0394-3',
]]

"""Build metadata
"""
meta = pd.DataFrame()
meta[0] = data.columns
meta[1] = [
    #'pex19\u0394',
    #'pex19\u0394',
    #'pex19\u0394msp1\u0394',
    #'pex19\u0394msp1\u0394',
    'pex3\u0394',
    'pex3\u0394',
    'pex3\u0394',
    'pex3\u0394msp1\u0394',
    'pex3\u0394msp1\u0394',
    'pex3\u0394msp1\u0394',
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
    'pex3\u0394':'darkgrey',
    'pex3\u0394msp1\u0394':'black'
}

data_scaled, data_labeled = xp.prep_data(data, meta, gene_scale=True)




xp.pca(
    data_scaled,
    meta,
    sample_colors,
    n_components=2,
    save_fig='/Users/jordan/Desktop/nuebel_2020/pex3/pca_all.png')

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
    figsize = (3,6)
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
    '/Users/jordan/Desktop/nuebel_2020/pex3/heatmap.png',
    dpi=1800,
    bbox_inches='tight'
)

"""Volcano Plot
"""
label_points = {
    'PEX3':[-3.85,	0.4],
    'PEX11':[1.45,	1.95],
    'PEX2':[1.55,	3.3],
    'PEX15':[1.65,	2.2],
    'MSP1':[-4,	    1.75],
}

xp.volcano(
    data[[
        'pex3\u0394-1',
        'pex3\u0394-2',
        'pex3\u0394-3',
        'pex3\u0394msp1\u0394-1',
        'pex3\u0394msp1\u0394-2',
        'pex3\u0394msp1\u0394-3'
    ]],
    meta.loc[4:],
    'pex3\u0394msp1',
    'pex3\u0394-',
    alpha=0.4,
    y_threshold=2,
    x_threshold=[-1,1],
    highlight_points=[pex_list + ['MSP1']],
    highlight_color=['red'],
    whitegrid=True,
    label_points=label_points,
    figsize=(5,4)
)

l_box = matplotlib.patches.Rectangle((-4.5,2), 3.5, 4.15, alpha=0.4, color="lightblue", zorder=-10)
r_box = matplotlib.patches.Rectangle((1,2), 1.75, 4.15, alpha=0.4, color="lightblue", zorder=-10)
ax = plt.gca()
ax.add_patch(l_box)
ax.add_patch(r_box)

# Save and show figure
plt.savefig(
    '/Users/jordan/Desktop/nuebel_2020/pex3/pex3/volcano_plot.png',
    dpi=1800,
    bbox_inches='tight'
)
