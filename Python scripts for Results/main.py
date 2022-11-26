import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import coresets
import algorithms

# Scanpy
matplotlib.use('Agg')
sc.settings.verbosity = 2             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True
#sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300)

#var_names = pd.read_csv('C:/Users/Khalid-IIIS/Desktop/pbmc68k_analysis/pbmc68k_all_genes.txt', sep=',', usecols=range(1,2), header=None)

var_names = pd.read_csv('/home/khalid/pbmc68k_analysis/pbmc68k_all_genes.txt', sep=',', usecols=range(1,2), header=None)

#############################    Compute Coreset    ########################################################

cores = 1
adata = sc.read_text('/home/khalid/pbmc68k_analysis/pbmc68k_full.txt', delimiter = ",")
#sc.pp.subsample(adata, fraction=0.15)
adata.var_names_make_unique()
n = adata.X.shape[0]
nn = int(n/cores)
coreset_size = int(5000/cores)
genes = list(adata.var_names)

X = pd.DataFrame()
for i in range(cores):
    data = adata.X[i * nn:(i + 1) * nn, :]
    km_coreset_gen = coresets.KMeansCoreset(data)
    C, w = km_coreset_gen.generate_coreset(int(coreset_size))
    S = C * w[:, np.newaxis]
    X = X.append(pd.DataFrame(data=S, columns=genes), ignore_index=True)
print(X.shape)
X.to_csv('/home/khalid/pbmc68k_analysis/pbmc68k_full_temp.txt', sep=',', index=False)

#adata = sc.read_text('C:/Users/Khalid-IIIS/Desktop/pbmc68k_analysis/pbmc68k_coreset_10k_combine.txt', delimiter = ",")
adata = sc.read_text('/home/khalid/pbmc68k_analysis/pbmc68k_full_temp.txt', delimiter = ",")

adata.var_names = var_names.values.transpose().ravel()
adata.var_names_make_unique('.')
print(adata)

#adata.X = adata.X.astype(np.int)
adata.raw = adata

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='cell_ranger') # for pbmc68k
#adata = adata[:, adata.var['highly_variable']]
#print(adata)

sc.tl.pca(adata, svd_solver='arpack', n_comps=50)  # 50 for pbmc68k
#sc.pl.pca(adata, color=['louvain'])
#sc.pl.pca_variance_ratio(adata, log=True)

# Finding Neighboring
pbmc68k_Neighbors = 10
pbmc68k_Distance = 0.3    # 0.3 for pbmc68k
pbmc68k_Resolution = 0.5   # 0.8 for pbmc68k

sc.pp.neighbors(adata, n_neighbors=pbmc68k_Neighbors, n_pcs=50) # use_rep='X' , n_pcs=10,  metric='correlation'
sc.tl.umap(adata, min_dist=pbmc68k_Distance)

# Louvain Clustering
sc.tl.louvain(adata, resolution=pbmc68k_Resolution)

#sc.tl.leiden(adata, resolution=pbmc68k_Resolution)
#sc.pl.umap(adata, color=['leiden'], s = 5)

# Find Marker Genes
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')

#sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', weights=None)
#sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False)

marker_genes = ['GZMK', 'CD8A', 'TNFRSF18', 'SIGLEC7', 'GNLY', 'LGALS3', 'CCR10',
                'CD4', 'CLEC4C', 'PF4', 'PTCRA', 'CD8B', 'ID3', 'CD79A']

# Plotting
sc.pl.umap(adata, color=['louvain'],  save='_LWCS_Solo_10K_Louvain.png')
sc.pl.rank_genes_groups(adata, save='_LWCS_Solo_10K_RankingGenes.pdf')
#sc.pl.dotplot(adata, marker_genes, groupby='louvain', save='_130K.png')
#sc.pl.stacked_violin(adata, marker_genes, groupby='louvain', save='_130K.png')
sc.pl.heatmap(adata, marker_genes, groupby='louvain', figsize=(5, 8), use_raw=True, vmin=-3, vmax=3, cmap='bwr', var_group_rotation=0, dendrogram=True, save='_LWCS_Solo_10K_Heatmap.pdf')

