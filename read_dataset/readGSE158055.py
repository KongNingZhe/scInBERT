import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from scanpy._utils import _empty
import os

path = "/home/zhangc/data/GSE158055/"
from scanpy._utils import _empty
adata = sc.read('/home/zhangc/data/GSE158055/GSE158055_covid19_matrix.mtx.gz',
                cache=True,
                cache_compression = _empty)
genes = pd.read_csv('/home/zhangc/data/GSE158055/GSE158055_covid19_features.tsv.gz', header=None)
cells = pd.read_csv('/home/zhangc/data/GSE158055/GSE158055_covid19_barcodes.tsv.gz', header=None)
adata.var_names = anndata.utils.make_index_unique(pd.Index(genes[0].values))
adata.obs_names = cells[0].values


def top_n_idx_sparse(matrix, n):
    '''Return index of top n values in each row of a sparse matrix'''
    top_n_idx = []
    for le, ri in zip(matrix.indptr[:-1], matrix.indptr[1:]):
        n_row_pick = min(n, ri - le)
        top_n_idx.append(matrix.indices[le + np.argpartition(matrix.data[le:ri], -n_row_pick)[-n_row_pick:]])
    return top_n_idx

with open(datapath+'filelist.txt') as f:
    lines = f.readlines()
    for i in range(2,len(lines)):
        datalist.append(lines[i].split('\t')[1])

out_filename = "GSE168215_train.txt"
out_f = open(out_filename,'w')

for data in datalist:
    adata = sc.read_10x_h5(datapath + data)
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=100)   
    sc.pp.filter_genes(adata, min_cells=3)

    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.var["mt"] = adata.var_names.str.match("^GRCh38______MT-|^GRCh38______mt-")
    adata.var["RP"] = adata.var_names.str.match("^GRCh38______RPS|^GRCh38______RPL|^GRCh38______Rps|^GRCh38______Rpl")
    adata.obs['percent_mito'] = np.sum(adata[:, adata.var["mt"]].X, axis=1) / np.array(adata.obs['n_counts']).reshape(-1,1)

    adata = adata[adata.obs['n_genes'] < 4000, :]
    adata = adata[adata.obs['percent_mito'] < 0.2, :]

    gene_list = np.array([l.split('______')[-1] for l in adata.var_names])

    top_n = top_n_idx_sparse(adata.X,126)
    for i in range(adata.n_obs):
        cell_gene = list(gene_list[top_n[i][np.argsort(-adata.X[i,top_n[i]].toarray())].reshape(-1)])
        s_cell_gene = '\t'.join(cell_gene) +'\n'
        out_f.write(s_cell_gene)

out_f.close()



