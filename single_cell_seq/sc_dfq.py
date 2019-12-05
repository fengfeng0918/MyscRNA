"""
Description: some test.
Date: 2019年12月2日13:21:21
Author: dfq
"""


import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from gprofiler import gprofiler

import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri

# 导入R包
"""
from rpy2.robjects.packages import importr
scran = importr("scran")
RColorBrewer = importr("RColorBrewer")
slingshot = importr("slingshot")
monocle = importr("monocle")
gam = importr("gam")
clusterExperiment = importr("clusterExperiment")
ggplot2 = importr("ggplot2")
plyr = importr("plyr")
MAST = importr("MAST")

# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
# %load_ext rpy2.ipython

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200, dpi_save=300)
# sc.logging.print_versions()
"""


# Set up data loading
def get_raw_data():
    # Data files  设置文件的识别类型
    sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
    sample_id_strings = ['3', '4', '5', '6', '7', '8']
    # file_base = '/home/fdong/data/sc_data/GSM283657'
    file_base = 'E:\sc_data\GSE92332\GSE92332_RAW\GSM283657\GSM283657'
    exp_string = '_Regional_'
    data_file_end = '_matrix.mtx'
    barcode_file_end = '_barcodes.tsv'
    gene_file_end = '_genes.tsv'
    # cc_genes_file = '~/data/Macosko_cell_cycle_genes.txt'
    cc_genes_file = 'E:\\sc_data\\Macosko_cell_cycle_genes.txt'

    # First data set load & annotation
    #Parse Filenames
    sample = sample_strings.pop(0)
    sample_id = sample_id_strings.pop(0)
    data_file = file_base+sample_id+exp_string+sample+data_file_end
    barcode_file = file_base+sample_id+exp_string+sample+barcode_file_end
    gene_file = file_base+sample_id+exp_string+sample+gene_file_end


    #Load data
    adata = sc.read(data_file, cache=True)    # 主要数据存储对象
    adata = adata.transpose()
    adata.X = adata.X.toarray()


    barcodes = pd.read_csv(barcode_file, header=None, sep='\t')
    genes = pd.read_csv(gene_file, header=None, sep='\t')

    #Annotate data
    barcodes.rename(columns={0:'barcode'}, inplace=True)
    barcodes.set_index('barcode', inplace=True)
    adata.obs = barcodes
    adata.obs['sample'] = [sample]*adata.n_obs
    adata.obs['region'] = [sample.split("_")[0]]*adata.n_obs
    adata.obs['donor'] = [sample.split("_")[1]]*adata.n_obs

    genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
    genes.set_index('gene_symbol', inplace=True)
    adata.var = genes

    for i in range(len(sample_strings)):
        # Parse Filenames
        sample = sample_strings[i]
        sample_id = sample_id_strings[i]
        data_file = file_base + sample_id + exp_string + sample + data_file_end
        barcode_file = file_base + sample_id + exp_string + sample + barcode_file_end
        gene_file = file_base + sample_id + exp_string + sample + gene_file_end

        # Load data
        adata_tmp = sc.read(data_file, cache=True)
        adata_tmp = adata_tmp.transpose()
        adata_tmp.X = adata_tmp.X.toarray()

        barcodes_tmp = pd.read_csv(barcode_file, header=None, sep='\t')
        genes_tmp = pd.read_csv(gene_file, header=None, sep='\t')

        # Annotate data
        barcodes_tmp.rename(columns={0: 'barcode'}, inplace=True)
        barcodes_tmp.set_index('barcode', inplace=True)
        adata_tmp.obs = barcodes_tmp
        adata_tmp.obs['sample'] = [sample] * adata_tmp.n_obs
        adata_tmp.obs['region'] = [sample.split("_")[0]] * adata_tmp.n_obs
        adata_tmp.obs['donor'] = [sample.split("_")[1]] * adata_tmp.n_obs

        genes_tmp.rename(columns={0: 'gene_id', 1: 'gene_symbol'}, inplace=True)
        genes_tmp.set_index('gene_symbol', inplace=True)
        adata_tmp.var = genes_tmp
        adata_tmp.var_names_make_unique()

        # Concatenate to main adata object
        adata = adata.concatenate(adata_tmp, batch_key='sample_id')
        # adata.var['gene_id'] = adata.var['gene_id-1']
        # adata.var.drop(columns=['gene_id-1', 'gene_id-0'], inplace=True)
        adata.obs.drop(columns=['sample_id'], inplace=True)
        adata.obs_names = [c.split("-")[0] for c in adata.obs_names]
        adata.obs_names_make_unique(join='_')


    #Assign variable names and gene id columns
    adata.var_names = [g.split("_")[1] for g in adata.var_names]
    adata.var['gene_id'] = [g.split("_")[1] for g in adata.var['gene_id']]

    return adata

def pre_qc(adata):
    # Quality control - calculate QC covariates
    adata.obs['n_counts'] = adata.X.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (adata.X > 0).sum(1)

    mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts']  # 计算每个细胞中线粒体基因比例
    return adata


def sc_qc_violin_plot(keys, group='sample', **kwds):
    if keys in adata.obs.columns:
        sc.pl.violin(adata, keys, groupby=group, cut=0, save="_{0}.png".format(keys), show=False, **kwds)
        sc.pl.violin(adata, keys, groupby=group, cut=0, save="_{0}.pdf".format(keys), show=False, **kwds)
    else:
        print("adata.obs.columns没有该{0}数据列".format(keys))


def sc_qc_scatter_plot(adata,x,y,color=None,save="00",**kwds):
    if (x in adata.obs.columns) and (y in adata.obs.columns):
        sc.pl.scatter(adata, x, y, color=color, save="_{0}.png".format(save), show=False, **kwds)
        sc.pl.scatter(adata, x, y, color=color, save="_{0}.pdf".format(save), show=False, **kwds)
    else:
        print("adata.obs.columns没有该{0}or{1}数据列".format(x,y))


def sc_top_gene_plot():
    aa = pd.DataFrame(adata.X)
    aa.index = adata.obs.index
    aa.columns = adata.var_names

    sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
    cc = aa[adata.obs["sample"] == "Duo_M2"].T
    gg = cc.loc[cc.sum(1).sort_values()[-20:].index,:]  # 总数量最多
    plt.boxplot(gg, vert=False, showfliers=False, labels=gg.index)

    dd = aa[adata.obs["sample"] == "Duo_M2"].T
    dd = dd/dd.sum(0)
    hh = dd.loc[(dd.T.describe().loc["50%",:]).sort_values()[-20:].index, :]  # 基因counts比例中位数最高
    fig5, ax = plt.subplots()
    ax.boxplot(hh, vert=False,showfliers=True,labels=hh.index,flierprops={"markersize":2,"markerfacecolor":"red","markeredgecolor":"red"})
    ax.set_title("Duo_M2")
    ax.set_xlabel("percent of total counts")
    plt.show()


    # all
    gg = aa.T.loc[aa.T.sum(1).sort_values()[-20:].index, :]
    plt.boxplot(gg, vert=False,showfliers=False,labels=gg.index)
    plt.show()



adata = get_raw_data()
adata = pre_qc(adata)
sc_qc_violin_plot("n_counts", group='sample', **{"size":1, "log":True, "show":False})
sc_qc_violin_plot("mt_frac", group='sample', **{"size":1, "show":False})
sc_qc_scatter_plot(adata, 'n_counts', 'n_genes', color='mt_frac',save=1)
sc_qc_scatter_plot(adata[adata.obs['n_counts']<10000], 'n_counts', 'n_genes', color='mt_frac',save=2)


# Annotate the data sets
print(adata.obs['region'].value_counts())
print('')
print(adata.obs['donor'].value_counts())
print('')
print(adata.obs['sample'].value_counts())




# Quality control - plot QC metrics
#Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='sample', size=2, log=True, cut=0,save="n_counts.png")
t2 = sc.pl.violin(adata, 'mt_frac', groupby='sample',save="mt_frac.png")

#Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac',save="_P1.png")
p2 = sc.pl.scatter(adata[adata.obs['n_counts']<10000], 'n_counts', 'n_genes', color='mt_frac')


#Thresholding decision: counts
p3 = sb.distplot(adata.obs['n_counts'], kde=False)

p4 = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000], kde=False, bins=60)

p5 = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']>10000], kde=False, bins=60)




#Thresholding decision: genes
p6 = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)

p7 = sb.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=60)

# Filter cells according to identified QC thresholds:
print('Total number of cells: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_counts = 1500)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, max_counts = 40000)
print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

adata = adata[adata.obs['mt_frac'] < 0.2]
print('Number of cells after MT filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_genes = 700)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))


#Filter genes:
print('Total number of genes: {:d}'.format(adata.n_vars))

# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(adata, min_cells=20)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))



#Perform a clustering for scran normalization in clusters
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)