# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 09:44:12 2022

@author: jihon
"""


import numpy as np
import pandas as pd

ptable_1 = pd.read_csv('data/protein_table_9drugs.csv')
ptable_2 = pd.read_csv('data/protein_table_15drugs.csv')


def add_gene_name(data):   
    desp = data['Description'].values
    gene = []
    for s in desp:
        
        g = s[s.find('GN=') + 3:].split(' ')[0]
        gene.append(g)
    data['Gene Symbol'] = gene
    return data


def preprocess_single(data):
    columns = ['Accession',
               '# PSMs',
               'Gene Symbol']
    val_cols = [s for s in data.columns if 'Abundances (Normalized):' in s]
    val_cols = val_cols[:len(val_cols)-2]
    columns += val_cols
    new_data = data[columns]
    new_columns = columns[:3] + ['Abundances_Pool_{}'.format(i) for i in range(9)]
    new_data.columns = new_columns
    return new_data
    

def combine_double(ptable_1, ptable_2):  
    index = list(set(ptable_1['Accession']) | set(ptable_2['Accession']))
    
    val_cols_1 = [s for s in ptable_1.columns if 'Abundances (Normalized):' in s]
    val_cols_2 = [s for s in ptable_2.columns if 'Abundances (Normalized):' in s]
    
    x1 = ptable_1[val_cols_1]
    x2 = ptable_2[val_cols_2]
    r1 = np.mean(x1.iloc[:,[-1,-2]], axis = 1)
    r2 = np.mean(x2.iloc[:,[-1,-2]], axis = 1)
    
    for i in range(x1.shape[0]):
        x1.iloc[i,:] /= r1[i]
    for i in range(x2.shape[0]):
        x2.iloc[i,:] /= r2[i]
    
    result = []
    for ind in index:
        w1 = np.where(ptable_1['Accession'] == ind)[0]
        if len(w1) > 0:
            w1 = w1[0]
            gene = ptable_1['Gene Symbol'][w1]
            psm1 = ptable_1['# PSMs'][w1]
            if (psm1 <= 5) and (r1[w1] <= 500):
                val1 = [np.nan] * 9
            else:
                val1 = list(x1.iloc[w1,:9])
        else:
            psm1 = 0
            val1 = [np.nan] * 9
        w2 = np.where(ptable_2['Accession'] == ind)[0]
        if len(w2) > 0:
            w2 = w2[0]
            gene = ptable_2['Gene Symbol'][w2]
            psm2 = ptable_2['# PSMs'][w2]
            if (psm2 <= 5) and (r2[w2] <= 500):
                val2 = [np.nan] * 9
            else:
                val2 = list(x2.iloc[w2,:9])
        else:
            psm2 = 0
            val2 = [np.nan] * 9
        psm = psm1 + psm2
        val = np.array(val1 + val2)
        val /= np.nanmax(val)
        item = [ind, gene, psm] + list(val)
        result.append(item)
    result = pd.DataFrame(result)
    result.columns = ['Accession', 'Gene Symbol', '# PSMs'] + ['Abundances_Pool_{}'.format(i) for i in range(len(val))]
    return result


ptable_1 = add_gene_name(ptable_1)
ptable_2 = add_gene_name(ptable_2)

new_ptable_1 = preprocess_single(ptable_1)
new_ptable_2 = preprocess_single(ptable_2)

new_ptable_1.to_csv('data/preprocessed_9drugs.csv', index = False)
new_ptable_2.to_csv('data/preprocessed_15drugs.csv', index = False)

combined_table = combine_double(ptable_1, ptable_2)
combined_table.to_csv('data/preprocessed_combined.csv', index = False)