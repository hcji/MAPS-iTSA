# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 10:04:33 2022

@author: jihon
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
from core.core import *

'''
1. Generating Sensing Matrix
'''

# pool_matrix = generate_sensing_matrix(n_pools = 9, n_drugs = 15, n_replicates = 3)
# pool_matrix = pd.DataFrame(pool_matrix)
# pool_matrix.columns = ['Palbociclib', 'Panobinostat',
#                        'Raltitrexed', 'Methotrexate', 'Vemurafenib', 'Fimepinostat',
#                        'Olaparib', 'Bafetinib', 'SCIO-469', 'OTS964', 'SL-327', 'Abemaciclib',
#                        'CCT137690', 'Belumosudil', 'Parthenolide']

# Here we use the pre-generated sensing matrix, because the generation has randomness
pool_matrix = pd.read_excel('data/sensing_matrix_15drugs.xlsx')
pool_matrix[np.isnan(pool_matrix)] = 0


'''
2. Drug Target Analysis with LASSO
'''

# read protein table
protein_table = pd.read_csv('data/preprocessed_15drugs.csv')

# processing
scores, fold_changes = post_analysis(protein_table, pool_matrix, drug_num = 3)

# save results
scores.to_csv('results/independent_15drugs.csv')


'''
3. Plot Results
'''

# verification drugs
drug_targets = {'Palbociclib': ['CDK4', 'CDK6'], 
                'Panobinostat': ['HDAC1', 'HDAC2'],
                'Fimepinostat': ['HDAC1'],
                'Raltitrexed': ['TYMS'], 
                'Methotrexate': ['DHFR'], 
                'Vemurafenib': ['BRAF'],
                'SCIO-469': ['MAPK14'],
                'SL-327': ['MAP2K1', 'MAP2K2']}


for d, genes in drug_targets.items():
    plot_results(d, scores, fold_changes, genes, fc_thres = 1.1, score_thres=0.2, top_markers = 8)


# exploratory drugs
drugs = ['Bafetinib', 'Abemaciclib', 'OTS964', 'Olaparib', 'CCT137690', 'Belumosudil', 'Parthenolide']

for i, d in enumerate(drugs):
    plot_results(d, scores, fold_changes, fc_thres = 1.1, score_thres=0.15, top_markers = 8)


'''
4. Plot Boxplots
'''

drugs = ['Panobinostat', 'Fimepinostat', 'Bafetinib', 'Abemaciclib']
genes = ['HDAC1', 'HDAC1', 'BRAF', 'CSNK2A2']

for i in range(len(drugs)):
    plot_boxplot(genes[i], drugs[i], protein_table, pool_matrix)


'''
5. Plot Target Ranking
'''

plt_data = []
for drug, target in drug_targets.items():
    for t in target:
        w = np.where(scores['Gene Symbol'] == t)[0][0]
        if drug not in scores.columns:
            rank = np.nan
            mark = str('not included')
            plt_data.append([drug, t, rank, mark])
            continue
        score = scores[drug][w]
        rank = len(np.where(scores[drug] >= score)[0])
        if rank > 25:
            mark = 'rank > 25'
        elif rank > 10:
            mark = '10 > rank >= 25'
        elif rank > 5:
            mark = '5 > rank >= 10'
        elif rank > 1:
            mark = '1 > rank >= 5'
        else:
            mark = 'rank = 1'
        plt_data.append([drug, t, rank, mark])
plt_data = pd.DataFrame(plt_data)
plt_data.columns = ['drug', 'target', 'rank', 'mark']

flatui = ['red', 'pink', 'orange', 'purple', 'navy']
my_cmap = ListedColormap(sns.color_palette(flatui).as_hex())

plt.figure(dpi = 300, figsize = (12, 5))
sns.scatterplot(x='drug', y='target', hue='mark', data=plt_data, s = 400, palette = my_cmap.colors,
                hue_order = ['rank = 1', '1 > rank >= 5', '5 > rank >= 10', '10 > rank >= 25', 'not included'])
plt.xticks(fontsize = 20, rotation=20)
plt.yticks(fontsize = 20)
plt.xlabel('')
plt.ylabel('')
plt.legend(loc='upper right', bbox_to_anchor=(1.39, 1), fontsize = 20, markerscale = 2)




