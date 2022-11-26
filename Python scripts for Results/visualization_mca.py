#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 09:25:56 2021

@author: mac
"""

import scanpy as sc
import numpy as np
import pandas as pd
from geosketch import gs
from sklearn import preprocessing
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import accuracy_score

sc.settings.verbosity = 0            # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_header()
sc.settings.set_figure_params(dpi=300)

def measurePerformance(y_test, y_pred, method):
    conf_matrix = confusion_matrix(y_test, y_pred)
    '''
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.matshow(conf_matrix, cmap=plt.cm.Oranges, alpha=0.3)
    for i in range(conf_matrix.shape[0]):
        for j in range(conf_matrix.shape[1]):
            ax.text(x=j, y=i,s=conf_matrix[i, j], va='center', ha='center', size='xx-large')
 
    plt.xlabel('Predicted', fontsize=18)
    plt.ylabel('Actuals', fontsize=18)
    plt.title('Confusion Matrix', fontsize=18)
    plt.show()
    filepng = 'neurons_confusion_' + method + '.png'
    filepdf = 'neurons_confusion_' + method + '.pdf'
    fig.savefig(filepng, dpi=300)
    fig.savefig(filepdf, dpi=300)
    '''

    microavg = precision_score(y_test, y_pred, average='micro')
    macroavg = precision_score(y_test, y_pred, average='macro')
    print(microavg, macroavg)
    
    print('Precision')
    microavg = precision_score(y_test, y_pred, average='micro')
    macroavg = precision_score(y_test, y_pred, average='macro')
    print(microavg, macroavg)
    print('Recall')
    microavg = recall_score(y_test, y_pred, average='micro')
    macroavg = recall_score(y_test, y_pred, average='macro')
    print(microavg, macroavg)
    print('F-Score')
    microavg = f1_score(y_test, y_pred, average='micro')
    macroavg = f1_score(y_test, y_pred, average='macro')
    print(microavg, macroavg)
    
    bal_acc = balanced_accuracy_score(y_test, y_pred)
    print(bal_acc)
    
    acc = accuracy_score(y_test, y_pred)
    print(acc)

    #apr = average_precision_score(y_test, y_pred, average="micro")
    #print(apr)
    
    #roc_auc = roc_auc_score(y_test, y_pred, average = 'macro', multi_class='ovr')
    #print(roc_auc)

def transferLabel(X_train, y_train):
    # Create KNN classifier
    knn = KNeighborsClassifier(n_neighbors = 5)# Fit the classifier to the data
    knn.fit(X_train, y_train)
    Y_pred = knn.predict(X_full)
    return Y_pred

def printFrequencies(Y):
    (unique, counts) = np.unique(Y, return_counts=True)
    frequencies = np.asarray((unique, counts)).T
    print(frequencies)
    
def saveClusters(labels, saveImage):
    adata = sc.read_text(fileName, delimiter = ',')
    sc.pp.neighbors(adata, n_pcs=50, n_neighbors=10) # n_neighbors=10 for zeisel and default for pbmc3k
    sc.tl.leiden(adata, resolution=0.8) # 0.1 for zeisel and default for pbmc3k
    sc.tl.umap(adata, min_dist=0.1) # 0.1 for zeisel and 0.3 for pbmc3k
    #sc.tl.tsne(adata, perplexity=30)
    
    adata.obs['leiden'] = labels
    adata.obs['leiden'] = adata.obs['leiden'].astype('category')
    sc.pl.umap(adata, color='leiden', s=10, legend_fontsize = 'large', save=saveImage+'.png')
    sc.pl.umap(adata, color='leiden', s=10, legend_fontsize = 'large', save=saveImage+'.pdf')
    #sc.pl.tsne(adata, color='leiden', s=10, legend_fontsize = 'large', save=saveImage)

X_full = np.loadtxt('/home/khalid/Datasets/mca_fbpca75.txt', delimiter=',')
#Y_full = pd.read_csv('/home/khalid/Datasets/.txt', delimiter=',', header=None, usecols=range(1,2))
#le = preprocessing.LabelEncoder()
#Y_full = le.fit_transform(Y_full.iloc[:,0])
Y_full = np.loadtxt('/home/khalid/Datasets/mca_label.txt', delimiter=',').astype(int)
printFrequencies(Y_full)

coresetSize = 5000

fileName = "/home/khalid/Datasets/sample.txt"

###   Uniform   ###
print('Uniform')
for i in range(10):
    idx = np.random.randint(X_full.shape[0], size=coresetSize)
    X_uniform = X_full[idx]
    Y_uniform = Y_full[idx].astype(int)
    Y_pred = transferLabel(X_uniform, Y_uniform) 
    #printFrequencies(Y_uniform)
    #printFrequencies(Y_pred)
    measurePerformance(Y_full, Y_pred, 'uniform')
    #saveClusters(Y_pred, saveImage='_neurons_uniform')

###   GeoSketch   ###
print('GeoSketch')
for i in range(3):
    sketch_idx = gs(X_full, coresetSize, replace=False)
    X_sketch = X_full[sketch_idx]
    Y_sketch = Y_full[sketch_idx].astype(int)
    Y_pred = transferLabel(X_sketch, Y_sketch) 
    #printFrequencies(Y_sketch)
    #printFrequencies(Y_pred)
    measurePerformance(Y_full, Y_pred, 'geo-Sketch')
    #saveClusters(Y_pred, saveImage='_neurons_geoSketch')

###   Sc-Coreset  ######
print('SC-Coreset')
coreset_idx = np.loadtxt('/home/khalid/Datasets/mca_coreset_5000_idx.txt', delimiter=',').astype(int)
coreset_idx -= 1
X_coreset = X_full[coreset_idx]
Y_coreset = Y_full[coreset_idx].astype(int)
Y_pred = transferLabel(X_coreset, Y_coreset) 
#printFrequencies(Y_coreset)
printFrequencies(Y_pred)
measurePerformance(Y_full, Y_pred, 'sc-Coreset')
#saveClusters(Y_pred, saveImage='_pbmc68k_scCoreset')
