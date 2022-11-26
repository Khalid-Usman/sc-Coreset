#/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 08:47:56 2020

@author: mac
"""

import numpy as np
from hausdorff import hausdorff_distance
from geosketch import gs
import matplotlib.pyplot as plt

#data = pd.read_csv('/Users/mac/Desktop/pbmc3k/pbmc3k_pca50.txt', sep=',', header=None)
#X = data.values()
#X = np.loadtxt('/Users/mac/Desktop/pbmc3k/pbmc3k_pca50.txt', delimiter=',') # pbmc3k
#X = np.loadtxt('/home/khalid/Datasets/pbmc68k_fbpca50.txt', delimiter=' ') # pbmc68k
#X = np.loadtxt('/Users/mac/Desktop/pbmc3k/pbmc3k_pca50.txt', delimiter=',') # pbmc3k
#X = np.loadtxt('/Users/mac/Desktop/zeisel/zeisel_fbpca10.txt', delimiter=',') # zeisel
X = np.loadtxt('/home/khalid/Datasets/neurons_fbpca100.txt', delimiter=' ') # mca

coresetSizes = [1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
#coresetSizes = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

def printError(X, Y):
    #print(f"Hausdorff distance test: {hausdorff_distance(X, Y, distance='manhattan')}")
    #print(f"Hausdorff distance test: {hausdorff_distance(X, Y, distance='euclidean')}")
    #print(f"Hausdorff distance test: {hausdorff_distance(X, Y, distance='chebyshev')}")
    #print(f"Hausdorff distance test: {hausdorff_distance(X, Y, distance='cosine')}")
    return hausdorff_distance(X, Y, distance='cosine')

for i in coresetSizes:
    
    #coresetSize = (i+1)*100
    coresetSize = i
    print(coresetSize)
    '''
    #### Geo Sketch  ####
    smallest = 1
    largest = 0
    total = 0
    average = 0
    for i in range(10):
        sketch_index = gs(X, coresetSize, replace=False)
        Y = X[sketch_index]
        dist = printError(X,Y)
        total += dist
        if dist < smallest:
            smallest = dist
        if dist > largest:
            largest = dist
    average = total/10
    print(average)
    print(smallest)
    print(largest)
    
    ###   Uniform  ######
    smallest = 1
    largest = 0
    total = 0
    average = 0
    for i in range(10):
        idx = np.random.randint(X.shape[0], size=coresetSize)
        Y = X[idx]
        dist = printError(X,Y)
        total += dist
        if dist < smallest:
            smallest = dist
        if dist > largest:
            largest = dist
    average = total/10
    print(average)
    print(smallest)
    print(largest)
    '''

    ###   Sc-Coreset  ######
    #if coresetSize == 350:
    #    continue
    fileName = '/home/khalid/Datasets/neurons_coreset_' + str(coresetSize) + '_1.txt'
    Y = np.loadtxt(fileName, delimiter=',')
    dist = printError(X,Y)
    print(dist)

'''
#  Plot Hausdorff Distance
x = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
y_geoSketch =       [0.87, 0.83, 0.81, 0.79, 0.78, 0.77, 0.74, 0.74, 0.73, 0.73]
y_geoSketch_lower = [0.84, 0.81, 0.79, 0.77, 0.76, 0.76, 0.73, 0.72, 0.71, 0.71]
y_geoSketch_upper = [0.89, 0.85, 0.84, 0.80, 0.80, 0.79, 0.76, 0.77, 0.74, 0.75]
y_uniform =         [0.79, 0.75, 0.73, 0.71, 0.70, 0.69, 0.68, 0.67, 0.66, 0.66]
y_uniform_lower =   [0.77, 0.72, 0.71, 0.69, 0.68, 0.67, 0.67, 0.65, 0.65, 0.65]
y_uniform_upper =   [0.81, 0.78, 0.77, 0.72, 0.71, 0.72, 0.71, 0.70, 0.68, 0.67]
y_scCoreset =       [0.78, 0.76, 0.74, 0.70, 0.69, 0.68, 0.67, 0.65, 0.62, 0.61]
y_scCoreset_lower = [0.77, 0.75, 0.73, 0.69, 0.68, 0.67, 0.66, 0.64, 0.61, 0.60]
y_scCoreset_upper = [0.79, 0.77, 0.75, 0.71, 0.70, 0.69, 0.69, 0.66, 0.63, 0.62]

fig, ax = plt.subplots()

#ax.fill_between(x, y_geoSketch_lower, y_geoSketch_upper, alpha=0.2)
ax.plot(x, y_geoSketch, '-', label = 'Geo-Sketch')
#ax.fill_between(x, y_uniform_lower, y_uniform_upper, alpha=0.2)
ax.plot(x, y_uniform, '-', label = 'Uniform Sampling')
#ax.fill_between(x, y_scCoreset_lower, y_scCoreset_upper, alpha=0.2)
ax.plot(x, y_scCoreset, '-', label = 'Sc-Coreset')

ax.yaxis.grid(False)
plt.xticks(fontsize='x-large')
plt.yticks(fontsize='x-large')
ax.legend(fontsize='x-large')
ax.legend(loc='center left', fontsize='x-large', bbox_to_anchor=(1, 0.5))

plt.savefig('/Users/mac/Desktop/pbmc68k_withoutStd_hausdorff.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/Users/mac/Desktop/pbmc68k_withoutStd_hausdorff.png', dpi=300, bbox_inches='tight')
'''
