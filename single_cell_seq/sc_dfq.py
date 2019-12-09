"""
Description: some test.
Date: 2019年12月2日13:21:21
Author: dfq
"""

import os
import scanpy as sc
import numpy as np
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt

# import scipy as sp
# from matplotlib import rcParams
# from matplotlib import colors
# from gprofiler import gprofiler
#
# import rpy2.rinterface_lib.callbacks
# import logging
# from rpy2.robjects import pandas2ri
# import anndata2ri

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
    file_base = '/home/fdong/data/sc_data/GSM283657'
    # file_base = 'E:\sc_data\GSE92332\GSE92332_RAW\GSM283657\GSM283657'
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

def sc_qc_violin_plot(adata, keys, group='sample', **kwds):
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

def sc_qc_dist_plot(data,keys,*arg,**kwargs):
    fig,ax = plt.subplots()
    sb.distplot(data,*arg,**kwargs)
    fig.savefig("figures/{}.png".format(keys))
    fig.savefig("figures/{}.pdf".format(keys))
    print("figures/{} has been ploted!".format(keys))
    plt.close()

def sc_top_gene_plot(adata,sample_strings):

    def _top_plot(data,name):
        """
        ploting.
        :param data: pd.DataFrame
        :param name: save filename
        :return:
        """
        fig, ax = plt.subplots()
        ax.boxplot(data, vert=False, showfliers=True, labels=hh.index,
                   flierprops={"markersize": 2, "markerfacecolor": "red", "markeredgecolor": "red"})
        ax.set_title(name)
        ax.set_xlabel("percent(%) of total counts")
        fig.savefig("figures/{0}_top20.png".format(name))
        fig.savefig("figures/{0}_top20.pdf".format(name))
        print("{0}_top20  has been saved!".format(name))
        plt.close()

    # sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
    # 生成所有数据集
    aa = pd.DataFrame(adata.X)
    aa.index = adata.obs.index
    aa.columns = adata.var_names

    # select top20 genes for each sample
    for sample in sample_strings:
        # 生成样本的数据集 如：Duo_M2
        dd = aa[adata.obs["sample"] == sample].T
        dd = dd/dd.sum(0)*100
        hh = dd.loc[dd.T.quantile([0.5]).T.sort_values(by=0.5)[-20:].index, :]  # 基因counts比例中位数最高
        _top_plot(hh, sample)

    # select top20 genes for all sample
    hh = aa.T.loc[aa.quantile([0.5]).T.sort_values(by=0.5)[-20:].index, :]
    _top_plot(hh, "all")

def sc_qc_filter(adata, **kwargs):
    before_cells = adata.shape[0]
    print("Total number of cells: [{0}]".format(before_cells))
    # Filter cells according to identified QC thresholds:
    for key,value in kwargs.items():
        if key == "min_counts":
            adata = adata[adata.obs['n_counts'] > value]
        elif key == "min_genes":
            adata = adata[adata.obs['n_genes'] > value]
        elif key == "max_counts":
            adata = adata[adata.obs['n_counts'] < value]
        elif key == "max_genes":
            adata = adata[adata.obs['n_genes'] < value]
        elif key =="mt_frac":
            adata = adata[adata.obs['mt_frac'] < value]
        else:
            print("Have nothing to do with {}={}.".format(key,value))
        after_cells = adata.shape[0]
        gap = before_cells - after_cells
        before_cells = after_cells

        print("Number of cells after {0} filter:[{1}]\tand [{2}] cells filtered out"
              .format(key,after_cells,gap))
    # sc.pp.filter_cells(adata, max_counts = 40000)
    # adata = adata[adata.obs['n_counts'] < 40000]

    sc.pp.filter_genes(adata, min_cells=20)
    return adata

def get_size_factors(adata_pp):

    # 保存数据 供R分析
    adata_pp.obs.to_csv("obs.txt",sep="\t",index_label="cell")
    aa = pd.DataFrame(adata_pp.X)
    # aa.index = adata_pp.obs.index
    aa.columns = adata_pp.var_names
    aa.to_csv("X.txt",sep="\t")
    os.system("Rscript Normalization.R")

    size_factors= pd.read_csv("size_factors.txt",header=0)
    return size_factors

def to_visualization(adata):
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata)

    sc.tl.tsne(adata, n_jobs=12)  # Note n_jobs works for MulticoreTSNE, but not regular implementation)
    sc.tl.umap(adata)
    sc.tl.diffmap(adata)
    sc.tl.draw_graph(adata)

    sc.pl.pca_scatter(adata, color='n_counts',show=False, save="pca.png")
    sc.pl.tsne(adata, color='n_counts',show=False, save="tsne.png")
    sc.pl.umap(adata, color='n_counts',show=False, save="umap.png")
    sc.pl.diffmap(adata, color='n_counts', components=['1,2', '1,3'],show=False, save="diffmap.png")
    sc.pl.draw_graph(adata, color='n_counts',show=False, save="draw_graph.png")

def cell_cycle_scoring(adata,cc_genes_file):
    # Score cell cycle and visualize the effect:
    cc_genes = pd.read_table(cc_genes_file, delimiter='\t')
    s_genes = cc_genes['S'].dropna()
    g2m_genes = cc_genes['G2.M'].dropna()

    s_genes_mm = [gene.lower().capitalize() for gene in s_genes]
    g2m_genes_mm = [gene.lower().capitalize() for gene in g2m_genes]

    s_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
    g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)

    sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False,show=False,save="S_score_G2M_score.png")
    sc.pl.umap(adata, color='phase', use_raw=False,show=False,save="cell_cycle_phase.png")


def main():
    sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
    cc_genes_file = '~/data/Macosko_cell_cycle_genes.txt'
    # cc_genes_file = 'E:\\sc_data\\Macosko_cell_cycle_genes.txt'

    # 获得数据
    adata = get_raw_data()
    adata = pre_qc(adata)
    #
    # # qc 画图
    # ## 查看在样本中的分布
    # sc_qc_violin_plot("n_counts", group='sample', **{"size": 1, "log": True, "show": False})
    # sc_qc_violin_plot("mt_frac", group='sample', **{"size": 1, "show": False})
    #
    # ## 查看数据之间的关系
    # sc_qc_scatter_plot(adata, 'n_counts', 'n_genes', color='mt_frac', save=1)
    # sc_qc_scatter_plot(adata[adata.obs['n_counts'] < 10000], 'n_counts', 'n_genes', color='mt_frac', save=2)
    #
    # ## 质控，确定删选的阈值
    # sc_qc_dist_plot(adata.obs['n_counts'], "n_counts", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_counts'][adata.obs['n_counts'] < 4000], "n_counts_less_4000", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_counts'][adata.obs['n_counts'] > 10000], "n_counts_more_10000", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_genes'], "n_genes", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_genes'][adata.obs['n_genes'] < 1000], "n_genes_less_1000", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['mt_frac'], "mt_frac", kde=False, bins=60)
    # # 以上得到n_counts范围1500~40000；n_genes范围为>700;mt_frac范围为<0.2
    #
    # filter cells
    adata = sc_qc_filter(adata, min_counts=1500, max_counts=40000, min_genes=700, mt_frac=0.2)
    #
    # Filter genes:
    print('Total number of genes: {:d}'.format(adata.n_vars))
    # Min 20 cells - filters out 0 count genes
    sc.pp.filter_genes(adata, min_cells=20)
    print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
    #
    #
    # # select top genes
    # sc_top_gene_plot(adata, sample_strings)


    # normalization

    # TODO 解释以下代码意义
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)
    # Normalization 得到size_factors
    size_factors = get_size_factors(adata_pp)
    del adata_pp
    size_factors.index = adata.obs.index
    adata.obs["size_factors"] = size_factors

    sc_qc_scatter_plot(adata, 'size_factors', 'n_counts',save="size_factors_n_counts")
    sc_qc_scatter_plot(adata, 'size_factors', 'n_genes',save="size_factors_n_genes")
    sc_qc_dist_plot(adata.obs['size_factors'],'size_factors', bins=50, kde=False)

    # Keep the count data in a counts layer
    adata.layers["counts"] = adata.X.copy()
    # Normalize adata
    adata.X /= adata.obs['size_factors'].values[:, None]
    sc.pp.log1p(adata)

    # Batch Correction
    # sc.pp.combat(adata, key='sample')

    # Highly Variable Genes
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
    print('\n', 'Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
    sc.pl.highly_variable_genes(adata,show=False, save="highly_variable_genes.png")




if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    main()






