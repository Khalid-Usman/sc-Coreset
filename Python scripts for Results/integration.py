#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 07:34:42 2021

@author: mac
"""

import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import datetime

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white')

#adata_ref = sc.read_text('/Users/mac/pbmc3k/pbmc3k_processed.txt', delimiter=',')
#adata_ref_genes = pd.read_csv('/Users/mac/pbmc3k/pbmc3k_processed_genes.txt', header=None)
#adata_ref_genes = list(adata_ref_genes.iloc[:,0])
#adata_ref.var_names = adata_ref_genes
#adata_ref.var_names.unique()
#print(adata_ref)

adata_ref = sc.datasets.pbmc3k_processed()  # this is an earlier version of the dataset from the pbmc3k tutorial

adata = sc.read_text('/Users/mac/pbmc68k/pbmc68k_HVG.txt', delimiter=' ')
adata_genes = pd.read_csv('/Users/mac/pbmc68k/pbmc68k_genes.txt', header=None, usecols=range(1,2))
adata_genes = list(adata_genes.iloc[:1000,0])
adata.var_names = adata_genes
adata.var_names.unique()
print(adata)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.louvain(adata)
sc.pp.subsample(adata, n_obs=5000)

labels = pd.read_csv('/Users/mac/pbmc68k/pbmc68k_label.txt', header=None)
labels = labels.iloc[:,0]
adata.obs['bulk_labels'] = labels

X_train = np.loadtxt('/Users/mac/pbmc68k/pbmc68k_coreset_5000.txt', delimiter=',')
idx = np.loadtxt('/Users/mac/pbmc68k/pbmc68k_coreset_idx_5000.txt')
adata.obsm['X_pca'] = X_train
adata.obs['louvain'] = labels[idx]
adata.obs['louvain'] = adata.obs['louvain'].astype('category')

#adata = sc.datasets.pbmc68k_reduced()

start = datetime.datetime.now()

var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref, color='louvain')

sc.tl.ingest(adata, adata_ref, obs='louvain')
adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix colors
sc.pl.umap(adata, color=['louvain', 'bulk_labels'], wspace=0.5)

adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
adata_concat.obs.louvain = adata_concat.obs.louvain.astype('category')
adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories, inplace=True)  # fix category ordering
adata_concat.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix category colors
sc.pl.umap(adata_concat, color=['batch', 'louvain'])

end = datetime.datetime.now()
print((end-start).total_seconds())

# Using BKNN
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')  # running bbknn 1.3.6
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['batch', 'louvain'])