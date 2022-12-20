#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:25:30 2021

@author: daniel
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt


mas10_original = {'N': 185_053, 'grid': ['25x19x19', '50x38x38', '100x76x76', '200x152x152', '400x304x304'], 'time': np.array([66.82, 117.42, 329.32, 2112.38, 2139.41])}
mas8_original = {'N': 373_336, 'grid': ['25x19x19', '50x38x38', '100x76x76', '200x152x152', '400x304x304'], 'time': np.array([128.99, 221.07, 637.32, 4228.14, 4255.59])}
mas5_original = {'N': 1_464_542, 'grid': ['25x19x19', '50x38x38', '100x76x76'], 'time': np.array([508.37, 930.29, 4074.60])}

mas10_lowtol = {'N': 185_053, 'grid': ['25x19x19', '50x38x38', '100x76x76', '200x152x152', '400x304x304'], 'time': np.array([9.63, 37.06, 260.09, 1994.44, 1995.75])}
mas8_lowtol = {'N': 373_336, 'grid': ['25x19x19', '50x38x38', '100x76x76', '200x152x152', '400x304x304'], 'time': np.array([18.12, 69.89, 509.74, 4058.39, 4058.70])}
mas5_lowtol = {'N': 1_464_542, 'grid': ['25x19x19', '50x38x38', '100x76x76'], 'time': np.array([77.95, 309.49, 2044.94])}

mas10_sparse = {'time': np.array([1.25, 56.82, 57.98, 152.76, 154.12, 162.83, 164.62, 171.00, 173.65, 210.40]),
                'range': np.linspace(0.5, 5, 10)}
mas8_sparse = {'time': np.array([2.54, 120.93, 123.39, 271.54, 274.36, 285.69, 289.38, 293.80, 299.20, 329.70]),
               'range': np.linspace(0.5, 5, 10)}
mas5_sparse = {'time': np.array([9.84, 411.72, 420.96, 526.52, 537.42, 1726.64, 1741.31, 1757.70, 1779.51, 1825.22]),
               'range': np.linspace(0.5, 5, 10)}


data = {'mas10': (mas10_original, mas10_lowtol, mas10_sparse),
        'mas8': (mas8_original, mas8_lowtol, mas8_sparse),
        'mas5': (mas5_original, mas5_lowtol, mas5_sparse)}
colors = ['red', 'blue', 'green']

file = open("timing_data.pkl", "wb")
pickle.dump(data, file)
file.close()

# file = open("timing_data.pkl", "rb")
# output = pickle.load(file)

fig, ax = plt.subplots(figsize=(10, 8))
handles1 = []
for i, sample in enumerate(data.keys()):
#    handles1.append(ax.plot(np.arange(len(data[sample][0]['time']))+1, data[sample][0]['time']/60/24,
#                    color=colors[i], label=f'{sample}, {data[sample][0]["N"]}', marker='.')[0])
#    ax.plot(np.arange(len(data[sample][1]['time']))+1, data[sample][1]['time']/60/24, color=colors[i], marker='.', linestyle='dashed')
    ax.plot(data[sample][2]['range'], data[sample][2]['time']/data[sample][0]["N"], color=colors[i], marker='.', linestyle='dotted')
   
handles2 = []    
handles2.append(ax.plot(0, 0, color='k', marker='.', label='Original')[0])
handles2.append(ax.plot(0, 0, color='k', marker='.', linestyle='dashed', label='low tol')[0])
handles2.append(ax.plot(0, 0, color='k', marker='.', linestyle='dotted', label='sparse')[0])
#first_legend = ax.legend(handles=handles1, fontsize=16, loc='upper left')
#ax.add_artist(first_legend)
ax.legend(handles=handles2, fontsize=16, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.1))

ax.set_xlim(0.5, 5.5)
ax.set_xlabel('Gridstep', fontsize=20)
ax.set_xticks(ticks=[1, 2, 3, 4, 5])
ax.set_xticklabels([grid.replace('x', '\n') for grid in mas10_original['grid']], fontsize=16)
ax.set_ylabel('seconds/N', fontsize=20)
#plt.ylabel('Minutes/N$_{wds}$', fontsize=20)
#plt.yticks(ticks=np.arange(0, 80, 12), fontsize=18)
ax.grid(b=1, axis='y')

plt.show()

