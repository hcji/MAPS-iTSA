# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 10:38:45 2022

@author: DELL
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from core.core import *


n_pools, n_replicates = 9, 3

max_rep = []

for n_drugs in tqdm(np.arange(5, 29, 2)):
    max_r = []
    for i in range(15):
        pool_matrix = generate_sensing_matrix(n_pools = n_pools, n_drugs = n_drugs, n_replicates = n_replicates)
        corr = np.dot(pool_matrix.T, pool_matrix) - n_replicates * np.eye(n_drugs)
        max_r.append(np.max(corr))
    max_rep.append(max_r)

max_rep = np.array(max_rep).T
max_num = np.mean(max_rep, axis = 1)

x = np.array(np.repeat(np.arange(5,29,2), 15))
y = max_rep.reshape((len(x),))

plt.figure(dpi = 300)
sns.violinplot(x = x, y = y)
plt.xlabel('number of drugs')
plt.ylabel('max overlapped pools')
