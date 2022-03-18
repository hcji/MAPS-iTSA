# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 08:25:14 2022

@author: jihon
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

panob = pd.read_csv('data/t_test_fimep.csv')
fimep = pd.read_csv('data/t_test_panob.csv')

protein_table = pd.read_csv('data/preprocessed_15drugs.csv')


def add_gene_name(data):
    gene = []
    for i in data.index:
        acc = panob['Accession'][i]
        w = np.where(protein_table['Accession'] == acc)[0][0]
        g = protein_table['Gene Symbol'][w]
        gene.append(g)
    data['Gene Symbol'] = gene
    return data

panob = add_gene_name(panob)
fimep = add_gene_name(fimep)


def plot_ranking(data, gene):
    w = np.array([i for i in data.index if data.loc[i, 'Gene Symbol'] in gene])
    data = data.sort_values(by = '-logPval', ascending = False)
    data = data.reset_index(drop = True)
    
    plt.figure(dpi = 300)
    plt.scatter(1 + np.arange(len(data)), data.loc[:,'-logPval'], s = 25, color = '#C0C0C0')
    plt.scatter(1 + w, data.loc[w,'-logPval'], color='#FF0000', s = 25)
    
    texts = []
    for i, s in enumerate(gene):
        x = w[i]
        y = data.loc[x,'-logPval']
        texts.append(plt.text(x+200, y, s, color = 'red', fontsize = 14, weight = 'semibold'))
    adjust_text(texts, force_points=0.01, force_text=0.01)

    plt.xlabel('rank of t-test', fontsize = 18)
    plt.ylabel('-log10 p value', fontsize = 18)
    plt.xticks(fontsize = 15, rotation=20)
    plt.yticks(fontsize = 15)


plot_ranking(panob, 'HDAC1')
plot_ranking(fimep, 'HDAC1')

