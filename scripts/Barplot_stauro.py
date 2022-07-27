# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:37:22 2022

@author: DELL
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kins = pd.read_csv('data/UniprotID_family-Kinase.txt', sep = '\t')
kins = kins['ID'].values
atpb = pd.read_csv('data/UniprotID_ATPbinding.txt', sep = '\t')
atpb = atpb['ID'].values
tars = np.union1d(kins, atpb)

x, y1, y2 = [], [], []
for f in os.listdir('results'):
    cel = f.split('_')[1]
    tab = pd.read_csv('results/{}'.format(f))
    sig = tab[tab['Staurosporine'] >= 0.15]['Accession'].values
    trues = np.intersect1d(sig, tars)
    falses = np.setdiff1d(sig, tars)
    x.append(cel)
    y1.append(len(trues))
    y2.append(len(falses))


df = pd.DataFrame(data={'kinase & ATP binding': y1, 'others': y2})
df.index = x

plt.rcParams["figure.dpi"] = 300
plt.rcParams["font.size"] = 15
ax = df.plot(kind='bar', stacked=True, figsize=(8, 6), rot=0, xlabel='Cell type', ylabel='Count')
for c in ax.containers:
    labels = np.round([v.get_height() if v.get_height() > 0 else '' for v in c]).astype(int)
    ax.bar_label(c, labels=labels, label_type='center')
