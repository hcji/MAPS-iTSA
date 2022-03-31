# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 08:19:58 2022

@author: jihon
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sko.GA import GA
from seaborn import heatmap


m, n, a = 9, 15, 3

def fun(v):
    d1 = mat(v)
    l1 = np.dot(d1.T, d1) - a * np.eye(n)
    l1 = np.sum(l1 ** 2)
    l2 = np.sum(d1, axis = 1)
    l2 = np.sum((l2 - np.mean(l2)) ** 2)
    return l1 + l2


def mat(v):
    d = v.reshape((m,n))
    d1 = np.array([np.argsort(d[:,i]) for i in range(d.shape[1])]).T
    d1[d1 < a] = 1
    d1[d1 >= a] = 0
    return d1  


# randomly initial
v0 = np.random.random(m * n)
adapting_0 = np.reshape(v0, (m,n))
sensing_0 = mat(v0)
corr_0 = np.dot(sensing_0.T, sensing_0) - a * np.eye(n)

plt.figure(dpi = 300)
heatmap(adapting_0, cbar = False, cmap = 'coolwarm', linewidths = 0.5, linecolor = 'black')
plt.xticks([])
plt.yticks([])

plt.figure(dpi = 300)
heatmap(sensing_0, cbar = False, cmap = ['#FFFFFF', '#FF0000'], linewidths = 0.5, linecolor = 'black')
plt.xticks([])
plt.yticks([])

plt.figure(dpi = 300)
heatmap(-corr_0, cbar = False, cmap = 'BuGn_r', linewidths = 0.5, linecolor = 'black')
plt.xticks([])
plt.yticks([])


# optimizing with GA
ga = GA(func = fun, n_dim = m*n, max_iter = 500, lb = np.repeat(0, m*n), ub = np.repeat(1, m*n))
x, y = ga.run()

adapting = np.reshape(x, (m,n))
sensing = mat(x)
corr = np.dot(sensing.T, sensing) - a * np.eye(n)

plt.figure(dpi = 300)
heatmap(adapting, cbar = False, cmap = 'coolwarm', linewidths = 0.5, linecolor = 'black')
plt.xticks([])
plt.yticks([])

plt.figure(dpi = 300)
heatmap(sensing, cbar = False, cmap = ['#FFFFFF', '#FF0000'], linewidths = 0.5, linecolor = 'black')
plt.xticks([])
plt.yticks([])

plt.figure(dpi = 300)
heatmap(-corr, cbar = False, cmap = 'BuGn_r', linewidths = 0.5, linecolor = 'black')
plt.xticks([])
plt.yticks([])
