# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 13:50:13 2021

@author: jihon
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from sko.GA import GA
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from scipy.stats import ttest_ind
from sklearn import linear_model
from tqdm import tqdm
from adjustText import adjust_text


def generate_sensing_matrix(n_pools = 9, n_drugs = 15, n_replicates = 3):
    '''
    Task: 
        Generate sensing matrix, which indicate how the drugs pool.
    Parameters:
        n_pools: int, the number of pools of the assay
        n_drugs: int, the number of drugs to be tested
        n_replicate: int, the number of replicates of the testing for each drug
    '''
    m, n, a = n_pools, n_drugs, n_replicates
    
    # optimize function, calculate the total correlations
    def fun(v):
        d1 = mat(v)
        l1 = np.dot(d1.T, d1) - a * np.eye(n)
        l1 = np.sum(l1 ** 2)
        l2 = np.sum(d1, axis = 1)
        l2 = np.sum((l2 - np.mean(l2)) ** 2)
        return l1 + l2
    
    # convert vector to sensing matrix
    def mat(v):
        d = v.reshape((m,n))
        d1 = np.array([np.argsort(d[:,i]) for i in range(d.shape[1])]).T
        d1[d1 < a] = 1
        d1[d1 >= a] = 0
        return d1    
    
    # use GA to optimize the vector
    ga = GA(func = fun, n_dim = m*n, max_iter = 500, lb = np.repeat(0, m*n), ub = np.repeat(1, m*n))
    x, y = ga.run()
    pool = mat(x)
    
    # plot history
    his = pd.DataFrame(ga.all_history_Y)
    fig, ax = plt.subplots(2, 1, figsize=(5,5), dpi=300)
    ax[0].plot(his.index, his.values, '.', color='red')
    his.min(axis=1).cummin().plot(kind='line')
    plt.xlabel('Index', fontsize = 12)
    plt.show()    
    return pool


def predict_targets(drugs, threshold = 0.3):
    '''
    Task: 
        Predict the potential targets for the drugs based on the structural similairy.
    Parameters:
        drugs, drug list.
    '''
    
    # load drugbank and remove items related to the testing drugs
    drugbank = pd.read_csv('drugbank/DrugBank_DTI.csv')
    k = np.array([n not in list(drugs['Name']) for n in drugbank['Name']])
    drugbank = drugbank[k]
    drugbank = drugbank.reset_index(drop=True)
    
    # calculat fingerprints of drugs in drugbank
    fps = []
    for smi in drugbank['SMILES']:
        try:
            m = Chem.MolFromSmiles(smi)
        except:
            m = None
        if m is None:
            fp = np.nan
        else:
            try:
                fp = AllChem.GetMorganFingerprintAsBitVect(m, 2)
            except:
                fp = np.nan
        fps.append(fp)
    
    # predict targets
    results = {}
    for i in drugs.index:
        drug, prot, score = [], [], []
        n, smi = drugs.loc[i, 'Name'], drugs.loc[i, 'SMILES']
        m = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2)
        for i, f in enumerate(fps):
            if type(f) is float:
                continue
            else:
                s = DataStructs.DiceSimilarity(fp, f)
                if s > threshold:
                    p = drugbank.loc[i, 'UniProt ID']
                    d = drugbank.loc[i, 'DrugBank ID']
                    if p not in prot:
                        drug.append(d)
                        prot.append(p)
                        score.append(s)
                    else:
                        j = prot.index(p)
                        if s > score[j]:
                            score[j] = s
                            drug[j] = d
        res = pd.DataFrame({'target': prot, 'reference': drug, 'similarity': score})
        res = res.sort_values('similarity')
        results[n] = res
    return results
    
    
def group_drugs(drugs, preds, pools):
    '''
    Task: 
        Divide testing drugs into different groups, each group for an assay
    Parameters:
        inp, the path of data file of drug list.
        preds, the predicted targets output by *predict_targets* function
        pools, the sensing matrix generated by *generate_sensing_matrix* function
    '''

    m, n, d = pools.shape[0], pools.shape[1], len(drugs)
     
    pred_tars = []
    for name in drugs['Name']:
        pred_tars.append(list(preds[name]['target']))

    def jaccard(a, b):
        unions = len(set(a).union(set(b)))
        intersections = len(set(a).intersection(set(b)))
        return intersections / (unions + 10 ** -6)

    def cal_tot_score(g):
        fp_sims = np.zeros((len(g), len(g)))
        for i, a in enumerate(g):
            for j, b in enumerate(g):
                tr1 = pred_tars[a]
                tr2 = pred_tars[b]
                fp_sims[i,j] = jaccard(tr1, tr2)
                fp_sims[j,i] = jaccard(tr1, tr2)
        fp_sims = fp_sims - np.eye(len(g))
        return np.sum(fp_sims ** 2)
    
    def optimizer(x):
        v = np.argsort(x)
        k = int(len(x) / n)
        s = 0
        for i in range(k):
            s += cal_tot_score(v[i*m : (i+1)*m])
        return s
     
    if len(drugs) % pools.shape[1] != 0:
        raise IOError('Drug number is not the integral multiple of the column number of sensing matrix')
    
    ub = np.ones(d)
    lb = np.zeros(d)
    ga = GA(func=optimizer, n_dim=d, size_pop=200, max_iter=100, lb=lb, ub=ub)
    best_v, minimum_similarity = ga.run()
    
    y_his = pd.DataFrame(ga.all_history_Y)
    fig, ax = plt.subplots(2, 1, dpi = 300)
    ax[0].plot(y_his.index, y_his.values, '.', color='red')
    y_his.min(axis=1).cummin().plot(kind='line')
    plt.xlabel('Index', fontsize = 12)
    plt.show()
    
    output = pd.DataFrame()
    v = np.argsort(best_v)
    k = int(len(v) / n)
    for i in range(k):
        output['Group {}'.format(i)] = list(drugs['Name'][v[i*n : (i+1)*n]])
    
    return output


def post_analysis(protein_table, pool_matrix, drug_num = 3):
    '''
    Task: 
        Identify drug-target interaction based on lasso algorithm
    Parameters:
        protein_table, the table of the quantitative proteomics.
        pool_matrix, the group matrix generated by *generate_sensing_matrix* function
        drug_num, which drug is used as reference
    '''
    
    X = np.nan_to_num(pool_matrix.values).astype(float)
    
    protein_table['Gene Symbol'] = protein_table['Gene Symbol'].astype(str)
    protein_table['# PSMs'] = [max( np.array(str(i).split(';')).astype(int)) for i in protein_table['# PSMs']]
    data = protein_table[protein_table['# PSMs'] > 2]
    data = data.reset_index(drop=True)
    
    k = np.array(['KRT' not in s for s in data['Gene Symbol'].values])
    data = data.iloc[k,:]
    data = data.reset_index(drop=True)
    
    cols = ['Abundances_Pool_{}'.format(i) for i in range(len(pool_matrix))]
    scores, fold_changes = [], []
    for i in tqdm(data.index):
        p = data.loc[i, 'Accession']
        g = data.loc[i, 'Gene Symbol']
        y = data.loc[i, cols].values.astype(float)
        kk = np.where(~np.isnan(y))[0]
        if len(kk) == 0:
            continue
        elif np.nanmax(y) == 0:
            continue
        else:
            y1 = y[kk]
            y1 = np.log2(y1 / np.min(y1))
        
        X1 = X[kk,:]
        '''
        for i in range(X1.shape[0]):
            sums = np.sum(X1[i,:])
            if sums == 0:
                sums = 1
            X1[i,:] /= sums
        '''
        _, _, coefs = linear_model.lars_path(X1, y1, positive=True, method="lasso")
        ind_1 = [np.where(X[:,i] == 1)[0] for i in range(X.shape[1])]
        ind_2 = [np.where(X[:,i] == 0)[0] for i in range(X.shape[1])]
        foldc = [np.nanmedian(y[ind_1[i]]) / np.nanmedian(y[ind_2[i]]) for i in range(X.shape[1])]
        
        xx = np.sum(np.abs(coefs.T), axis=1)
        xx /= xx[-1]
        
        '''
        x2 = ['Others', 'Others', 'Others', 'MTX', 'Others', 'MTX', 'Others', 'Others', 'MTX']
        y2 = np.log2(y)
        pltdata = pd.DataFrame({'pools': x2, 'log2_abundance': y2})
        plt.figure(dpi = 300)
        sns.boxplot(x="pools", y="log2_abundance", data=pltdata)
        plt.xticks(fontsize = 14)
        plt.yticks(fontsize = 14)
        plt.xlabel("group", fontsize=18)
        plt.ylabel("log2 abundance", fontsize=18)
        '''
        
        '''
        plt.figure(figsize = (4,6), dpi = 250)
        for j in range(coefs.shape[0]):
            if j == 1:
                plt.plot(xx, coefs[j,:], color='r', label='Panob')
            elif j == 5:
                plt.plot(xx, coefs[j,:], color='b', label='Fimep')
            else:
                plt.plot(xx, coefs[j,:], color='grey')
                
        plt.legend(fontsize = 20, loc='upper left')
        plt.ylim(0, 2.5)
        ymin, ymax = plt.ylim()
        # plt.vlines(xx, ymin, ymax, linestyle="dashed")
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        # plt.xlabel("|coef| / max|coef|", fontsize=14)
        # plt.ylabel("Coefficients", fontsize=14)
        # plt.axis("tight")
        plt.show()
        '''
        vip = coefs[:,drug_num] - np.min(coefs[np.nonzero(coefs[:,drug_num]),3])
        vip[vip < 0] = 0
        res_1 = [p, g] + list(vip)
        res_2 = [p, g] + list(foldc)
        scores.append(res_1)
        fold_changes.append(res_2)
    scores = pd.DataFrame(scores)
    scores.columns = ['Accession', 'Gene Symbol'] + list(pool_matrix.columns)
    fold_changes = pd.DataFrame(fold_changes)
    fold_changes.columns = ['Accession', 'Gene Symbol'] + list(pool_matrix.columns)
    return scores, fold_changes
    

def plot_results(drug, scores, fold_changes, true_targets=[], fc_thres = 1.05, score_thres=0.1, top_markers = 20):
    '''
    Task: 
        Scatter plot of the specific drug
    '''
    
    gene = scores['Gene Symbol'].values
    score = scores[drug].values
    fold_change = np.log2(fold_changes[drug].values)
    pltdata = pd.DataFrame({'Gene Symbol': gene, 'Score': score, 'Fold Change': fold_change})

    group = []
    for i in range(len(gene)):
        if (fold_change[i] < np.log2(fc_thres)) and (score[i] < score_thres):
            group.append('Not significant')
        elif (fold_change[i] >= np.log2(fc_thres)) and (score[i] < score_thres):
            group.append('Not significant')
        elif (fold_change[i] < np.log2(fc_thres)) and (score[i] >= score_thres):
            group.append('Not significant')
        else:
            group.append('Significant')
    pltdata['Group'] = group
    pltdata = pltdata.sort_values(by = 'Group', ascending=False)
    pltdata = pltdata.reset_index(drop = True)
    
    plt.figure(dpi = 250)
    # plt.title(drug, loc='left', fontsize = 20)
    flatui = ['#FF0000', '#C0C0C0']
    # scatter_cmap = ListedColormap(sns.color_palette(flatui).as_hex())
    
    sns.scatterplot(data=pltdata, x="Fold Change", y="Score", hue="Group",
                    palette = flatui,
                    hue_order=['Significant', 'Not significant'] )
    
    markers = pltdata[pltdata['Group'] == 'Significant']
    markers = markers.sort_values(by = 'Score', ascending=False)
    markers = markers.reset_index(drop = True)
    markers = markers.iloc[:min(len(markers), top_markers),:]
    texts = []
    for i in markers.index:
        x, y, s = markers.loc[i, 'Fold Change'], markers.loc[i, 'Score'], str(markers.loc[i, 'Gene Symbol']).split(';')[0]
        if s in true_targets:
            texts.append(plt.text(x, y, s, color = 'red', fontsize = 14, weight = 'semibold'))
        else:
            texts.append(plt.text(x, y, s, fontsize = 13))
    
    adjust_text(texts, force_points=0.2, force_text=0.2,
            expand_points=(1, 1), expand_text=(1, 1),
            arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    plt.axvline(x = np.log2(fc_thres), ls = '--', color = 'black', lw = 1)
    plt.axhline(y = score_thres, ls = '--', color = 'black', lw = 1)
    plt.legend(loc='upper left', fontsize = 12)
    plt.xticks(fontsize = 15, rotation = 20)
    plt.yticks(fontsize = 15)
    plt.xlabel('log2 FC', fontsize = 18)
    plt.ylabel('LASSO score', fontsize = 18)
    plt.show()
    pass


def plot_drugs(gene, scores, fold_changes, fc_thres = 1.05, score_thres=0.08, top_markers = 10):
    '''
    Task: 
        Scatter plot of the specific target
    '''
    
    try:
        i = np.where(scores['Gene Symbol'] == gene)[0][0]
    except:
        return None
    drug = scores.columns[2:]
    score = scores.iloc[i,2:]
    fold_change = fold_changes.iloc[i,2:]
    pltdata = pd.DataFrame({'Drug': drug, 'Score': score, 'Fold Change': fold_change})
    
    group = []
    for i in range(len(drug)):
        if (fold_change[i] < np.log2(fc_thres)) and (score[i] < score_thres):
            group.append('Others')
        elif (fold_change[i] >= np.log2(fc_thres)) and (score[i] < score_thres):
            group.append('Others')
        elif (fold_change[i] < np.log2(fc_thres)) and (score[i] >= score_thres):
            group.append('Others')
        else:
            group.append('Significant')
    pltdata['Group'] = group
    pltdata = pltdata.sort_values(by = 'Group', ascending=False)
    pltdata = pltdata.reset_index(drop = True)
    
    plt.figure(dpi = 250)
    sns.scatterplot(data=pltdata, x="Fold Change", y="Score", hue="Group")

    markers = pltdata[pltdata['Group'] == 'Significant']
    markers = markers.sort_values(by = 'Score', ascending=False)
    markers = markers.reset_index(drop = True)
    markers = markers.iloc[:min(len(markers), top_markers),:]
    texts = []
    for j in markers.index:
        x, y, s = markers.loc[j, 'Fold Change'], markers.loc[j, 'Score'], str(markers.loc[j, 'Drug']).split(';')[0]
        texts.append(plt.text(x, y, s, fontsize = 13))
    
    adjust_text(texts, force_points=0.2, force_text=0.2,
            expand_points=(1, 1), expand_text=(1, 1),
            arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    # plt.axvline(x = np.log2(fc_thres), ls = '--', color = 'black', lw = 1)
    plt.axhline(y = score_thres, ls = '--', color = 'black', lw = 1)
    plt.legend(loc='upper left', fontsize = 12)
    plt.xticks(fontsize = 15, rotation=45)
    plt.yticks(fontsize = 15)
    plt.xlabel('log2 FC', fontsize = 18)
    plt.ylabel('score', fontsize = 18)
    plt.show()


def plot_boxplot(gene, drug, protein_table, pool_matrix):
    protein_table['# PSMs'] = [max( np.array(str(i).split(';')).astype(int)) for i in  protein_table['# PSMs']]
    data = protein_table[protein_table['# PSMs'] > 4]
    data = data.reset_index(drop=True)
   
    cols = ['Abundances_Pool_{}'.format(i) for i in range(len(pool_matrix))]
    i = np.where(protein_table['Gene Symbol'] == gene)[0][0]
    j = np.where(pool_matrix[drug] == 1)[0]
    k = np.where(pool_matrix[drug] == 0)[0]
    y = data.loc[i, cols].values.astype(float)

    x2 = ['Others'] * len(y)
    y2 = np.log2(y)
    for jj in j:
        x2[jj] = drug
    pltdata = pd.DataFrame({'pools': x2, 'log2_abundance': y2})
    plt.figure(figsize=(6,5), dpi = 300)
    sns.boxplot(x="pools", y="log2_abundance", data=pltdata, order=[drug, 'Others'])
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 14)
    plt.xlabel("")
    plt.ylabel("log2 abundance", fontsize=16)
    print(ttest_ind(y2[j], y2[k]))


    
def merge_replicate(ptable_1, ptable_2):  
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
            val1 = list(x1.iloc[w1,:9])
            val1 = np.nan_to_num(val1, 0)
        else:
            psm1 = 0
            val1 = [0] * 9
        w2 = np.where(ptable_2['Accession'] == ind)[0]
        if len(w2) > 0:
            w2 = w2[0]
            gene = ptable_2['Gene Symbol'][w2]
            psm2 = ptable_2['# PSMs'][w2]
            val2 = list(x2.iloc[w2,:9])
            val2 = np.nan_to_num(val2, 0)
        else:
            psm2 = 0
            val2 = [0] * 9
        psm = psm1 + psm2
        val = np.array(val1) + np.array(val2)
        val /= np.nanmax(val)
        item = [ind, gene, psm] + list(val)
        result.append(item)
    result = pd.DataFrame(result)
    result.columns = ['Accession', 'Gene Symbol', '# PSMs'] + ['Abundances_Pool_{}'.format(i) for i in range(len(val))]
    return result



if __name__ == '__main__':
    
    pass