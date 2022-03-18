# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:14:26 2022

@author: jihon
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

method = ['TPP-TR', 'TPP-CCR', 'iTSA/PISA', 'Pooled-TPP (9*9)', 'Pooled-TPP (9*15)']
drugTMT = np.array([1/4, 1/2, 1, 9, 15])


plt.figure(dpi = 300, figsize = (14, 5))
sns.barplot(x=method, y=np.log2(drugTMT + 0.02), palette='pastel')
plt.xticks(fontsize = 18, rotation=20)
plt.yticks(fontsize = 18)
plt.xlabel('')
plt.ylabel('No. of drugs per TMT Expt.(Log2)', fontsize = 20)
plt.ylim(-4, 7)